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

#ifdef GRINS_HAVE_ANTIOCH

// GRINS
#include "grins/antioch_mixture_averaged_transport_mixture.h"

// Antioch
#include "antioch_config.h"
#include "antioch/vector_utils_decl.h"
#include "antioch/vector_utils.h"

// This class
#include "antioch_mixture_averaged_transport_mixture.C"

#include "grins/antioch_instantiation_macro.h"

INSTANTIATE_ANTIOCH_TRANSPORT(AntiochMixtureAveragedTransportMixture);

#ifdef ANTIOCH_HAVE_GSL
INSTANTIATE_ANTIOCH_KINETICS_THEORY_TRANSPORT(AntiochMixtureAveragedTransportMixture);
#endif // ANTIOCH_HAVE_GSL


#endif //GRINS_HAVE_ANTIOCH
