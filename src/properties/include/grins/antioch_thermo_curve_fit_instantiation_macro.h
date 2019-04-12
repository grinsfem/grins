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

#ifndef GRINS_ANTIOCH_THERMO_CURVE_FIT_INSTANTIATION_MACRO_H
#define GRINS_ANTIOCH_THERMO_CURVE_FIT_INSTANTIATION_MACRO_H

#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

#define INSTANTIATE_ANTIOCH_THERMO_CURVE_FIT(class_name)                \
  template class GRINS::class_name<Antioch::CEACurveFit<libMesh::Real> >; \
  template class GRINS::class_name<Antioch::NASA7CurveFit<libMesh::Real> >; \
  template class GRINS::class_name<Antioch::NASA9CurveFit<libMesh::Real> >

#define INSTANTIATE_ANTIOCH_THERMO_CURVE_FIT_WITH_CONDUCTIVITY(class_name,conductivity) \
  template class GRINS::class_name<Antioch::CEACurveFit<libMesh::Real>,GRINS::conductivity>; \
  template class GRINS::class_name<Antioch::NASA7CurveFit<libMesh::Real>,GRINS::conductivity>; \
  template class GRINS::class_name<Antioch::NASA9CurveFit<libMesh::Real>,GRINS::conductivity>

#endif // GRINS_HAVE_ANTIOCH

#endif // GRINS_ANTIOCH_THERMO_CURVE_FIT_INSTANTIATION_MACRO_H
