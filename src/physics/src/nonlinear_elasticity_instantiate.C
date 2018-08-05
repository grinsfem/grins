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

#include "elastic_membrane_base.C"
#include "elastic_membrane.C"
#include "elastic_cable_base.C"
#include "elastic_cable.C"
#include "elastic_cable_rayleigh_damping.C"
#include "elastic_membrane_rayleigh_damping.C"
#include "elastic_membrane_pressure.C"

#include "grins/hookes_law.h"
#include "grins/hookes_law_1d.h"
#include "grins/incompressible_plane_stress_hyperelasticity.h"
#include "grins/mooney_rivlin.h"

#include "grins/constant_pressure.h"
#include "grins/parsed_pressure.h"

template class GRINS::ElasticMembraneBase<GRINS::HookesLaw>;
template class GRINS::ElasticMembraneBase<GRINS::IncompressiblePlaneStressHyperelasticity<GRINS::MooneyRivlin> >;
template class GRINS::ElasticMembrane<GRINS::HookesLaw>;
template class GRINS::ElasticMembrane<GRINS::IncompressiblePlaneStressHyperelasticity<GRINS::MooneyRivlin> >;
template class GRINS::ElasticMembraneRayleighDamping<GRINS::HookesLaw>;
template class GRINS::ElasticMembraneRayleighDamping<GRINS::IncompressiblePlaneStressHyperelasticity<GRINS::MooneyRivlin> >;

template class GRINS::ElasticCableBase<GRINS::HookesLaw1D>;
template class GRINS::ElasticCable<GRINS::HookesLaw1D>;
template class GRINS::ElasticCableRayleighDamping<GRINS::HookesLaw1D>;

template class GRINS::ElasticMembranePressure<GRINS::ConstantPressure>;
template class GRINS::ElasticMembranePressure<GRINS::ParsedPressure>;
