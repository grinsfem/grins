//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
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
#include "grins/antioch_wilke_transport_mixture.h"

// Antioch
#include "antioch/vector_utils_decl.h"
#include "antioch/vector_utils.h"

// This class
#include "antioch_wilke_transport_mixture.C"

template class GRINS::AntiochWilkeTransportMixture<Antioch::StatMechThermodynamics<libMesh::Real>,
                                                   Antioch::MixtureViscosity<Antioch::SutherlandViscosity<libMesh::Real> >,
                                                   Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >,
                                                   Antioch::ConstantLewisDiffusivity<libMesh::Real> >;

template class GRINS::AntiochWilkeTransportMixture<Antioch::StatMechThermodynamics<libMesh::Real>,
                                                   Antioch::MixtureViscosity<Antioch::BlottnerViscosity<libMesh::Real> >,
                                                   Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >,
                                                   Antioch::ConstantLewisDiffusivity<libMesh::Real> >;

#endif //GRINS_HAVE_ANTIOCH
