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
#include "grins/antioch_evaluator.h"

// Antioch
#include "antioch/nasa_evaluator.h"
#include "antioch/nasa9_curve_fit.h"
#include "antioch/nasa7_curve_fit.h"
#include "antioch/cea_curve_fit.h"
#include "antioch/stat_mech_thermo.h"
#include "antioch/ideal_gas_thermo.h"

// This class
#include "antioch_evaluator.C"

template class GRINS::AntiochEvaluator<Antioch::NASA7CurveFit<libMesh::Real>,Antioch::StatMechThermodynamics<libMesh::Real> >;
template class GRINS::AntiochEvaluator<Antioch::NASA9CurveFit<libMesh::Real>,Antioch::StatMechThermodynamics<libMesh::Real> >;
template class GRINS::AntiochEvaluator<Antioch::CEACurveFit<libMesh::Real>,Antioch::StatMechThermodynamics<libMesh::Real> >;

template class GRINS::AntiochEvaluator<Antioch::NASA7CurveFit<libMesh::Real>, Antioch::IdealGasThermo<Antioch::NASA7CurveFit<libMesh::Real>,libMesh::Real> >;
template class GRINS::AntiochEvaluator<Antioch::NASA9CurveFit<libMesh::Real>, Antioch::IdealGasThermo<Antioch::NASA9CurveFit<libMesh::Real>,libMesh::Real> >;
template class GRINS::AntiochEvaluator<Antioch::CEACurveFit<libMesh::Real>, Antioch::IdealGasThermo<Antioch::CEACurveFit<libMesh::Real>,libMesh::Real> >;

#endif //GRINS_HAVE_ANTIOCH
