#ifndef ROSIHOOK_H
#define ROSIHOOK_H

#include<Eigen/StdVector>
#include "PhysicsHook.h"
#include <random>

class ROSIHook : public PhysicsHook
{
public:
    enum ContactMethod {
        CM_GS,
        CM_ROSI
    };

    ROSIHook(int numnails, ContactMethod method) : PhysicsHook(), requestedNails(numnails), contactMethod(method) {}

    void initializeNails();
    virtual void initSimulation();

    virtual void updateRenderGeometry()
    {
        renderQ = Q;
        renderF = F;        
    }

    bool computeBarycentricWeights(const Eigen::Vector2d &nailpos, int &face, double &a, double &b, double &c);
    Eigen::VectorXd nailConstraintNormal(int nail);

    void doROSI();
    void doGaussSeidel();

    virtual bool simulateOneStep();

    virtual void renderRenderGeometry(igl::viewer::Viewer &viewer)
    {
        viewer.data.set_mesh(renderQ, renderF);
        viewer.data.set_edges(nailsP, nailsE, nailsC);
    }

    void setNumNails(int nnails) { requestedNails = nnails; }
    void setContactMethod(ContactMethod method) { contactMethod = method; }

private:
    std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d> > nails;

    bool bounced;
    double g;
    double dt;
    double k;
    Eigen::MatrixXd Q;
    Eigen::MatrixXd Qorig;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::SparseMatrix<double> Minv;

    Eigen::MatrixXd renderQ;
    Eigen::MatrixXi renderF;
    Eigen::MatrixXd nailsP;
    Eigen::MatrixXi nailsE;
    Eigen::MatrixXd nailsC;

    int requestedNails;
    ContactMethod contactMethod;
};

#endif