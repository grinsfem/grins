#include "grins_solver.h"
#include "multiphysics_sys.h"

int main(int argc, char* argv[])
{

  std::string application_options; // TODO: ADD OPTIONS LATER

  GRINS::Solver<GRINS::MultiphysicsSystem> grins_solver(application_options);

  return 0;
}
