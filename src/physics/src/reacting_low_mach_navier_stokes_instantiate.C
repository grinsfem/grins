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


#include "grins_config.h"

#include "reacting_low_mach_navier_stokes_base.C"
#include "reacting_low_mach_navier_stokes.C"

#ifdef GRINS_HAVE_CANTERA

#include "grins/cantera_mixture.h"
#include "grins/cantera_evaluator.h"

template class GRINS::ReactingLowMachNavierStokesBase<GRINS::CanteraMixture,GRINS::CanteraEvaluator>;
template class GRINS::ReactingLowMachNavierStokes<GRINS::CanteraMixture,GRINS::CanteraEvaluator>;

#endif // GRINS_HAVE_CANTERA



#ifdef GRINS_HAVE_ANTIOCH

#include "grins/antioch_mixture_averaged_transport_mixture.h"
#include "grins/antioch_mixture_averaged_transport_evaluator.h"

#include "grins/antioch_constant_transport_mixture.h"
#include "grins/antioch_constant_transport_evaluator.h"

/* -------------------- ReactingLowMachNavierStokes -------------------- */
template class GRINS::ReactingLowMachNavierStokes<GRINS::AntiochMixtureAveragedTransportMixture<Antioch::StatMechThermodynamics<libMesh::Real>,
                                                                                                Antioch::SutherlandViscosity<libMesh::Real>,
                                                                                                Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >,
                                                                                                Antioch::ConstantLewisDiffusivity<libMesh::Real> >,
                                                  GRINS::AntiochMixtureAveragedTransportEvaluator<Antioch::StatMechThermodynamics<libMesh::Real>,
                                                                                                  Antioch::SutherlandViscosity<libMesh::Real>,
                                                                                                  Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >,
                                                                                                  Antioch::ConstantLewisDiffusivity<libMesh::Real> > >;

template class GRINS::ReactingLowMachNavierStokes<GRINS::AntiochMixtureAveragedTransportMixture<Antioch::StatMechThermodynamics<libMesh::Real>,
                                                                                                Antioch::BlottnerViscosity<libMesh::Real>,
                                                                                                Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >,
                                                                                                Antioch::ConstantLewisDiffusivity<libMesh::Real> >,
                                                  GRINS::AntiochMixtureAveragedTransportEvaluator<Antioch::StatMechThermodynamics<libMesh::Real>,
                                                                                                  Antioch::BlottnerViscosity<libMesh::Real>,
                                                                                                  Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >,
                                                                                                  Antioch::ConstantLewisDiffusivity<libMesh::Real> > >;

template class GRINS::ReactingLowMachNavierStokes<GRINS::AntiochMixtureAveragedTransportMixture<Antioch::StatMechThermodynamics<libMesh::Real>,
                                                                                                Antioch::KineticsTheoryViscosity<libMesh::Real,Antioch::GSLSpliner>,
                                                                                                Antioch::KineticsTheoryThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real>,libMesh::Real>,
                                                                                                Antioch::MolecularBinaryDiffusion<libMesh::Real,Antioch::GSLSpliner> >,
                                                  GRINS::AntiochMixtureAveragedTransportEvaluator<Antioch::StatMechThermodynamics<libMesh::Real>,
                                                                                                  Antioch::KineticsTheoryViscosity<libMesh::Real,Antioch::GSLSpliner>,
                                                                                                  Antioch::KineticsTheoryThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real>,libMesh::Real>,
                                                                                                  Antioch::MolecularBinaryDiffusion<libMesh::Real,Antioch::GSLSpliner> > >;

template class GRINS::ReactingLowMachNavierStokesBase<GRINS::AntiochMixtureAveragedTransportMixture<Antioch::StatMechThermodynamics<libMesh::Real>,
                                                                                          Antioch::SutherlandViscosity<libMesh::Real>,
                                                                                          Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >,
                                                                                          Antioch::ConstantLewisDiffusivity<libMesh::Real> >,
                                                      GRINS::AntiochMixtureAveragedTransportEvaluator<Antioch::StatMechThermodynamics<libMesh::Real>,
                                                                                            Antioch::SutherlandViscosity<libMesh::Real>,
                                                                                            Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >,
                                                                                            Antioch::ConstantLewisDiffusivity<libMesh::Real> > >;

template class GRINS::ReactingLowMachNavierStokesBase<GRINS::AntiochMixtureAveragedTransportMixture<Antioch::StatMechThermodynamics<libMesh::Real>,
                                                                                                    Antioch::BlottnerViscosity<libMesh::Real>,
                                                                                                    Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >,
                                                                                                    Antioch::ConstantLewisDiffusivity<libMesh::Real> >,
                                                      GRINS::AntiochMixtureAveragedTransportEvaluator<Antioch::StatMechThermodynamics<libMesh::Real>,
                                                                                                      Antioch::BlottnerViscosity<libMesh::Real>,
                                                                                                      Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >,
                                                                                                      Antioch::ConstantLewisDiffusivity<libMesh::Real> > >;

template class GRINS::ReactingLowMachNavierStokesBase<GRINS::AntiochMixtureAveragedTransportMixture<Antioch::StatMechThermodynamics<libMesh::Real>,
                                                                                                    Antioch::KineticsTheoryViscosity<libMesh::Real,Antioch::GSLSpliner>,
                                                                                                    Antioch::KineticsTheoryThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real>,libMesh::Real>,
                                                                                                    Antioch::MolecularBinaryDiffusion<libMesh::Real,Antioch::GSLSpliner> >,
                                                      GRINS::AntiochMixtureAveragedTransportEvaluator<Antioch::StatMechThermodynamics<libMesh::Real>,
                                                                                                      Antioch::KineticsTheoryViscosity<libMesh::Real,Antioch::GSLSpliner>,
                                                                                                      Antioch::KineticsTheoryThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real>,libMesh::Real>,
                                                                                                      Antioch::MolecularBinaryDiffusion<libMesh::Real,Antioch::GSLSpliner> > >;


template class GRINS::ReactingLowMachNavierStokes<GRINS::AntiochConstantTransportMixture<GRINS::ConstantConductivity>,
                                                  GRINS::AntiochConstantTransportEvaluator<Antioch::StatMechThermodynamics<libMesh::Real>, GRINS::ConstantConductivity> >;

template class GRINS::ReactingLowMachNavierStokes<GRINS::AntiochConstantTransportMixture<GRINS::ConstantPrandtlConductivity>,
                                                  GRINS::AntiochConstantTransportEvaluator<Antioch::StatMechThermodynamics<libMesh::Real>, GRINS::ConstantPrandtlConductivity> >;

template class GRINS::ReactingLowMachNavierStokes<GRINS::AntiochConstantTransportMixture<GRINS::ConstantConductivity>,
                                                  GRINS::AntiochConstantTransportEvaluator<Antioch::CEAEvaluator<libMesh::Real>, GRINS::ConstantConductivity> >;

template class GRINS::ReactingLowMachNavierStokes<GRINS::AntiochConstantTransportMixture<GRINS::ConstantPrandtlConductivity>,
                                                  GRINS::AntiochConstantTransportEvaluator<Antioch::CEAEvaluator<libMesh::Real>, GRINS::ConstantPrandtlConductivity> >;

template class GRINS::ReactingLowMachNavierStokesBase<GRINS::AntiochConstantTransportMixture<GRINS::ConstantConductivity>,
                                                      GRINS::AntiochConstantTransportEvaluator<Antioch::StatMechThermodynamics<libMesh::Real>, GRINS::ConstantConductivity> >;

template class GRINS::ReactingLowMachNavierStokesBase<GRINS::AntiochConstantTransportMixture<GRINS::ConstantPrandtlConductivity>,
                                                      GRINS::AntiochConstantTransportEvaluator<Antioch::StatMechThermodynamics<libMesh::Real>, GRINS::ConstantPrandtlConductivity> >;

template class GRINS::ReactingLowMachNavierStokesBase<GRINS::AntiochConstantTransportMixture<GRINS::ConstantConductivity>,
                                                      GRINS::AntiochConstantTransportEvaluator<Antioch::CEAEvaluator<libMesh::Real>, GRINS::ConstantConductivity> >;

template class GRINS::ReactingLowMachNavierStokesBase<GRINS::AntiochConstantTransportMixture<GRINS::ConstantPrandtlConductivity>,
                                                      GRINS::AntiochConstantTransportEvaluator<Antioch::CEAEvaluator<libMesh::Real>, GRINS::ConstantPrandtlConductivity> >;

#endif //GRINS_HAVE_ANTIOCH
