#include "grins_solver.h"
#include "low_mach_num_navier_stokes_sys.h"

int main(int argc, char* argv[])
{

  std::string application_options; // TODO: ADD OPTIONS LATER

  GRINS::Solver<GRINS::LowMachNumberNavierStokesSystem> grins_solver(application_options);

  return 0;
}
