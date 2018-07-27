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

#ifndef GRINS_ANTIOCH_INSTANTIATION_MACRO_H
#define GRINS_ANTIOCH_INSTANTIATION_MACRO_H

#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

#define INSTANTIATE_ANTIOCH_TRANSPORT_RAW(class_name,curve_fit)         \
  template class GRINS::class_name<curve_fit,                           \
                                   Antioch::StatMechThermodynamics<libMesh::Real>, \
                                   Antioch::SutherlandViscosity<libMesh::Real>, \
                                   Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >, \
                                   Antioch::ConstantLewisDiffusivity<libMesh::Real> >; \
  template class GRINS::class_name<curve_fit,                           \
                                   Antioch::StatMechThermodynamics<libMesh::Real>, \
                                   Antioch::BlottnerViscosity<libMesh::Real>, \
                                   Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >, \
                                   Antioch::ConstantLewisDiffusivity<libMesh::Real> >; \
  template class GRINS::class_name<curve_fit,                           \
                                   Antioch::IdealGasThermo<curve_fit,libMesh::Real>, \
                                   Antioch::SutherlandViscosity<libMesh::Real>, \
                                   Antioch::EuckenThermalConductivity<Antioch::IdealGasThermo<curve_fit,libMesh::Real> >, \
                                   Antioch::ConstantLewisDiffusivity<libMesh::Real> >; \
  template class GRINS::class_name<curve_fit,                           \
                                   Antioch::IdealGasThermo<curve_fit,libMesh::Real>, \
                                   Antioch::BlottnerViscosity<libMesh::Real>, \
                                   Antioch::EuckenThermalConductivity<Antioch::IdealGasThermo<curve_fit,libMesh::Real> >, \
                                   Antioch::ConstantLewisDiffusivity<libMesh::Real> >

#define INSTANTIATE_ANTIOCH_TRANSPORT(class_name)                       \
  INSTANTIATE_ANTIOCH_TRANSPORT_RAW(class_name,Antioch::CEACurveFit<libMesh::Real>); \
  INSTANTIATE_ANTIOCH_TRANSPORT_RAW(class_name,Antioch::NASA7CurveFit<libMesh::Real>); \
  INSTANTIATE_ANTIOCH_TRANSPORT_RAW(class_name,Antioch::NASA9CurveFit<libMesh::Real>)

#ifdef ANTIOCH_HAVE_GSL
#define INSTANTIATE_ANTIOCH_KINETICS_THEORY_TRANSPORT_RAW(class_name,curve_fit) \
  template class GRINS::class_name<curve_fit,                           \
                                   Antioch::StatMechThermodynamics<libMesh::Real>, \
                                   Antioch::KineticsTheoryViscosity<libMesh::Real,Antioch::GSLSpliner>, \
                                   Antioch::KineticsTheoryThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real>,libMesh::Real>, \
                                   Antioch::MolecularBinaryDiffusion<libMesh::Real,Antioch::GSLSpliner> >; \
  template class GRINS::class_name<curve_fit,                           \
                                   Antioch::IdealGasThermo<curve_fit,libMesh::Real>, \
                                   Antioch::KineticsTheoryViscosity<libMesh::Real,Antioch::GSLSpliner>, \
                                   Antioch::KineticsTheoryThermalConductivity<Antioch::IdealGasThermo<curve_fit,libMesh::Real>, libMesh::Real>, \
                                   Antioch::MolecularBinaryDiffusion<libMesh::Real,Antioch::GSLSpliner> >

#define INSTANTIATE_ANTIOCH_KINETICS_THEORY_TRANSPORT(class_name)       \
  INSTANTIATE_ANTIOCH_KINETICS_THEORY_TRANSPORT_RAW(class_name,Antioch::CEACurveFit<libMesh::Real>); \
  INSTANTIATE_ANTIOCH_KINETICS_THEORY_TRANSPORT_RAW(class_name,Antioch::NASA7CurveFit<libMesh::Real>); \
  INSTANTIATE_ANTIOCH_KINETICS_THEORY_TRANSPORT_RAW(class_name,Antioch::NASA9CurveFit<libMesh::Real>)

#endif // ANTIOCH_HAVE_GSL

#endif // GRINS_HAVE_ANTIOCH

#endif // GRINS_ANTIOCH_INSTANTIATION_MACRO_H
