//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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

#include "elastic_membrane.C"
#include "elastic_cable.C"

#include "grins/hookes_law.h"
#include "grins/hookes_law_1d.h"
#include "grins/incompressible_plane_stress_hyperelasticity.h"
#include "grins/mooney_rivlin.h"

template class GRINS::ElasticMembrane<GRINS::HookesLaw>;
template class GRINS::ElasticMembrane<GRINS::IncompressiblePlaneStressHyperelasticity<GRINS::MooneyRivlin> >;

template class GRINS::ElasticCable<GRINS::HookesLaw1D>;
