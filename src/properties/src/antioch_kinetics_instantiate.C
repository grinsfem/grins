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
#include "grins/antioch_kinetics.h"

// Antioch
#include "antioch/cea_curve_fit.h"
#include "antioch/nasa9_curve_fit.h"
#include "antioch/nasa7_curve_fit.h"

// This class
#include "antioch_kinetics.C"

template class GRINS::AntiochKinetics<Antioch::NASA7CurveFit<libMesh::Real> >;
template class GRINS::AntiochKinetics<Antioch::NASA9CurveFit<libMesh::Real> >;
template class GRINS::AntiochKinetics<Antioch::CEACurveFit<libMesh::Real> >;

#endif //GRINS_HAVE_ANTIOCH
