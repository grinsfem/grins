//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2019 Paul T. Bauman, Roy H. Stogner
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

#include "catalytic_wall_base.C"
#include "gas_recombination_catalytic_wall.C"
#include "gas_solid_catalytic_wall.C"

#ifdef GRINS_HAVE_CANTERA

#include "grins/cantera_mixture.h"

template class GRINS::CatalyticWallBase<GRINS::CanteraMixture>;
template class GRINS::GasRecombinationCatalyticWall<GRINS::CanteraMixture>;
template class GRINS::GasSolidCatalyticWall<GRINS::CanteraMixture>;

#endif // GRINS_HAVE_CANTERA


#ifdef GRINS_HAVE_ANTIOCH

#include "grins/antioch_chemistry.h"

template class GRINS::CatalyticWallBase<GRINS::AntiochChemistry>;
template class GRINS::GasRecombinationCatalyticWall<GRINS::AntiochChemistry>;
template class GRINS::GasSolidCatalyticWall<GRINS::AntiochChemistry>;

#endif //GRINS_HAVE_ANTIOCH
