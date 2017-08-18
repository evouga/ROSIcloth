#include "ROSIHook.h"
#include "EigenQP.h"

void ROSIHook::initializeNails()
{
    int nnails = requestedNails;
    unsigned int seed = (unsigned int)std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator(seed);
    std::uniform_real_distribution<double> distribution(-1.0, 1.0);
    nails.clear();
    for (int i = 0; i < nnails; i++)
    {
        Eigen::Vector2d pos(distribution(generator), distribution(generator));
        nails.push_back(pos);
    }

    nailsP.resize(2*nnails, 3);
    nailsE.resize(nnails, 2);
    nailsC.resize(nnails, 3);
    for (int i = 0; i < nnails; i++)
    {
        nailsP(2*i, 0) = nails[i][0];
        nailsP(2*i, 1) = nails[i][1];
        nailsP(2*i, 2) = 0;
        nailsP(2*i+1, 0) = nails[i][0];
        nailsP(2*i+1, 1) = nails[i][1];
        nailsP(2*i+1, 2) = -1;

        nailsE(i, 0) = 2 * i;
        nailsE(i, 1) = 2 * i + 1;

        nailsC(i, 0) = 1;
        nailsC(i, 1) = 0;
        nailsC(i, 2) = 0;
    }
}

void ROSIHook::initSimulation()
{
    Q.resize(4, 3);
    Q << -1, -1, 1,
        1, -1, 1,
        -1, 1, 1,
        1, 1, 1;
    V.resize(4, 3);
    V.setZero();
    F.resize(2, 3);
    F << 0, 1, 2,
        2, 1, 3;

    Qorig = Q;

    initializeNails();
    dt = 1e-6;
    g = -9.8;
    k = 100;
    bounced = false;

    Minv.resize(3 * Q.rows(), 3 * Q.rows());
    std::vector<Eigen::Triplet<double> > Minvvals;
    int nfaces = F.rows();
    for (int i = 0; i < nfaces; i++)
    {
        Eigen::Vector3i face = F.row(i);
        Eigen::Vector3d e1 = Q.row(face[1]) - Q.row(face[0]);
        Eigen::Vector3d e2 = Q.row(face[2]) - Q.row(face[0]);
        double area = 0.5*e1.cross(e2).norm();
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                Minvvals.push_back(Eigen::Triplet<double>(3 * face[j] + k, 3 * face[j] + k, area));
            }
        }
    }
    Minv.setFromTriplets(Minvvals.begin(), Minvvals.end());
}

bool ROSIHook::computeBarycentricWeights(const Eigen::Vector2d &nailpos, int &faceidx, double &a, double &b, double &c)
{
    int nfaces = F.rows();
    for (int i = 0; i < nfaces; i++)
    {
        Eigen::Vector3i face = F.row(i);
        Eigen::Vector3d e1 = Q.row(face[1]) - Q.row(face[0]);
        Eigen::Vector3d e2 = Q.row(face[2]) - Q.row(face[0]);
        Eigen::Matrix2d edges;
        edges << e1[0], e2[0],
            e1[1], e2[1];
        Eigen::Vector2d rhs = nailpos - Eigen::Vector2d(Q(face[0], 0), Q(face[0], 1));
        Eigen::Vector2d bary = edges.inverse()*rhs;
        if (bary[0] >= 0.0 && bary[1] >= 0.0 && bary[0] + bary[1] <= 1.0)
        {
            faceidx = i;
            a = 1.0 - bary[0] - bary[1];
            b = bary[0];
            c = bary[1];
            return true;
        }
    }
    return false;
}

Eigen::VectorXd ROSIHook::nailConstraintNormal(int nail)
{
    assert(nail >= 0 && nail < (int)nails.size());
    int faceidx;
    double w[3];
    bool found = computeBarycentricWeights(nails[nail], faceidx, w[0], w[1], w[2]);
    assert(found);
    Eigen::VectorXd result(3 * Q.rows());
    result.setZero();
    for (int i = 0; i < 3; i++)
    {
        result[3 * F(faceidx, i) + 2] = w[i];
    }
    double norm = sqrt(result.transpose()*(Minv*result));
    result /= norm;
    return result;
}

void ROSIHook::doGaussSeidel()
{
    bool done = false;
    int nnails = nails.size();
    Eigen::VectorXd Vflat(3 * V.rows());
    for (int i = 0; i < V.rows(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Vflat[3 * i + j] = V(i, j);
        }
    }
    while (!done)
    {
        done = true;
        for (int i = 0; i < nnails; i++)
        {
            Eigen::VectorXd normal = nailConstraintNormal(i);
            double relvel = normal.dot(Vflat);
            if (relvel < 0)
            {
                done = false;
                Vflat -= 2.0*relvel*Minv*normal;
            }
        }
    }
    for (int i = 0; i < V.rows(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            V(i,j) = Vflat[3 * i + j];
        }
    }
}

void ROSIHook::doROSI()
{
    int nnails = nails.size();
    Eigen::VectorXd Vflat(3 * V.rows());
    for (int i = 0; i < V.rows(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Vflat[3 * i + j] = V(i, j);
        }
    }

    bool done = false;
    while (!done)
    {
        std::vector<Eigen::VectorXd> violated;
        for (int i = 0; i < nnails; i++)
        {
            Eigen::VectorXd normal = nailConstraintNormal(i);
            double relvel = normal.dot(Vflat);
            if (relvel < 0)
            {
                violated.push_back(normal);
            }
        }

        if (violated.empty())
        {
            done = true;
        }
        else
        {
            Eigen::MatrixXd NV(Vflat.size(), violated.size());
            for (int i = 0; i < (int)violated.size(); i++)
            {
                for (int j = 0; j < (int)Vflat.size(); j++)
                {
                    NV(j, i) = violated[i][j];
                }
            }

            Eigen::MatrixXd halfObj = Minv * NV;
            Eigen::MatrixXd G = halfObj.transpose() * halfObj;
            double reg = 1e-6;
            for (int i = 0; i < (int)violated.size(); i++)
            {
                G(i, i) += reg;
            }
            Eigen::VectorXd g0 = 2.0*Vflat.transpose()*Minv*NV;
            Eigen::MatrixXd CE(violated.size(), 0);
            Eigen::VectorXd ce0(0);
            Eigen::MatrixXd CI(violated.size(), violated.size());
            CI.setIdentity();
            Eigen::VectorXd ci0(violated.size());
            ci0.setZero();
            Eigen::VectorXd lambda(violated.size());
            lambda.setZero();
            QP::solve_quadprog(G, g0, CE, ce0, CI, ci0, lambda);
            std::cout << "lambda: " << lambda.transpose() << std::endl;
            Eigen::VectorXd update = Minv * (NV*lambda);
            std::cout << "update: " << update.transpose() << std::endl;
            Vflat += update;
        }
    }
    for (int i = 0; i < V.rows(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            V(i, j) = Vflat[3 * i + j];
        }
    }
}

bool ROSIHook::simulateOneStep()
{
    Q += dt * V;
    int nfaces = F.rows();
    for (int i = 0; i < nfaces; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            int idx1 = F(i, j);
            int idx2 = F(i, (j + 1) % 3);
            Eigen::Vector3d edge = Q.row(idx1) - Q.row(idx2);
            Eigen::Vector3d edge0 = Qorig.row(idx1) - Qorig.row(idx2);
            double len = edge.norm();
            double len0 = edge0.norm();
            V.row(idx1) -= dt * k * (len - len0) * edge / len;
            V.row(idx2) += dt * k * (len - len0) * edge / len;
        }
    }
    for (int i = 0; i < V.rows(); i++)
    {
        V(i, 2) += dt * g;
    }

    if (Q(0, 2) < 0 && !bounced)
    {
        bounced = true;
        if (contactMethod == CM_GS)
            doGaussSeidel();
        else if(contactMethod == CM_ROSI)
            doROSI();
    }

    return false;
}

