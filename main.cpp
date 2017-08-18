// shut up warnings inside libigl
#ifdef _MSC_VER
#pragma warning(push, 0)   
#endif
#include <igl/viewer/Viewer.h>
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#include <thread>
#include "PhysicsHook.h"
#include "ROSIHook.h"

static PhysicsHook *hook = NULL;

static int numnails;

static ROSIHook::ContactMethod method;

void toggleSimulation()
{
    if (!hook)
        return;

    if (hook->isPaused())
        hook->run();
    else
        hook->pause();
}

void resetSimulation()
{
    if (!hook)
        return;

    ((ROSIHook *)hook)->setNumNails(numnails);
    ((ROSIHook *)hook)->setContactMethod(method);
    hook->reset();
}

bool drawCallback(igl::viewer::Viewer &viewer)
{
    if (!hook)
        return false;

    hook->render(viewer);
    return false;
}

bool keyCallback(igl::viewer::Viewer& viewer, unsigned int key, int modifiers)
{
    if (key == ' ')
    {
        toggleSimulation();
        return true;
    }
    return false;
}

bool initGUI(igl::viewer::Viewer &viewer)
{
    viewer.ngui->addButton("Run/Pause Sim", toggleSimulation);
    viewer.ngui->addButton("Reset Sim", resetSimulation);
    viewer.ngui->addVariable("Nails", numnails);
    viewer.ngui->addVariable("Method", method, true)->setItems({ "Gauss-Seidel", "ROSI" });
    viewer.screen->performLayout();
    return false;
}

int main(int argc, char *argv[])
{
  igl::viewer::Viewer viewer;

  numnails = 4;
  method = ROSIHook::ContactMethod::CM_GS;
  hook = new ROSIHook(numnails, method);
  hook->reset();
  
  viewer.data.set_face_based(true);
  viewer.core.is_animating = true;
  viewer.callback_key_pressed = keyCallback;
  viewer.callback_pre_draw = drawCallback;
  viewer.callback_init = initGUI;
  viewer.launch();
}
