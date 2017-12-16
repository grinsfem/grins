//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
// Copyright (C) 2010-2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-


#include "grins_config.h"

#include "grins/reacting_low_mach_navier_stokes_base.h"
#include "reacting_low_mach_navier_stokes.C"
#include "reacting_low_mach_navier_stokes_stab_base.C"
#include "reacting_low_mach_navier_stokes_spgsm_stab.C"

#ifdef GRINS_HAVE_CANTERA

#include "grins/cantera_mixture.h"
#include "grins/cantera_evaluator.h"

template class GRINS::ReactingLowMachNavierStokesBase<GRINS::CanteraMixture>;
template class GRINS::ReactingLowMachNavierStokes<GRINS::CanteraMixture,GRINS::CanteraEvaluator>;
template class GRINS::ReactingLowMachNavierStokesStabilizationBase<GRINS::CanteraMixture,GRINS::CanteraEvaluator>;
template class GRINS::ReactingLowMachNavierStokesSPGSMStabilization<GRINS::CanteraMixture,GRINS::CanteraEvaluator>;

#endif // GRINS_HAVE_CANTERA



#ifdef GRINS_HAVE_ANTIOCH

#include "grins/reacting_low_mach_navier_stokes_macro.h"

INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_CONSTANT_MIXTURE(ReactingLowMachNavierStokesBase);

INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_CONSTANT_MIXTURE_AND_CONSTANT_EVALUATOR(ReactingLowMachNavierStokes);
INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_CONSTANT_MIXTURE_AND_CONSTANT_EVALUATOR(ReactingLowMachNavierStokesStabilizationBase);
INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_CONSTANT_MIXTURE_AND_CONSTANT_EVALUATOR(ReactingLowMachNavierStokesSPGSMStabilization);

INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_MIXTURE_ONLY(ReactingLowMachNavierStokesBase);
INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_MIXTURE_AND_EVALUATOR(ReactingLowMachNavierStokes);
INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_MIXTURE_AND_EVALUATOR(ReactingLowMachNavierStokesStabilizationBase);
INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_MIXTURE_AND_EVALUATOR(ReactingLowMachNavierStokesSPGSMStabilization);

#endif //GRINS_HAVE_ANTIOCH
