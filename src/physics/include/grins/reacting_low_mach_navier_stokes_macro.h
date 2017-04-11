#ifndef GRINS_REACTING_LOW_MACH_NAVIER_STOKES_MACRO_H
#define GRINS_REACTING_LOW_MACH_NAVIER_STOKES_MACRO_H

#include "grins/antioch_mixture_averaged_transport_mixture.h"
#include "grins/antioch_mixture_averaged_transport_evaluator.h"

#include "grins/antioch_constant_transport_mixture.h"
#include "grins/antioch_constant_transport_evaluator.h"

#define INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_MIXTURE_AND_EVALUATOR(class_name) \
  /*
   *MixtureAveraged: StatMech thermo, Sutherland viscosity, Eucken conductivity, constant Lewis diffusivity
   */ \
template class GRINS::class_name<GRINS::AntiochMixtureAveragedTransportMixture<Antioch::StatMechThermodynamics<libMesh::Real>, \
                                                                               Antioch::SutherlandViscosity<libMesh::Real>, \
                                                                               Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >, \
                                                                               Antioch::ConstantLewisDiffusivity<libMesh::Real> >, \
                                 GRINS::AntiochMixtureAveragedTransportEvaluator<Antioch::StatMechThermodynamics<libMesh::Real>, \
                                                                                 Antioch::SutherlandViscosity<libMesh::Real>, \
                                                                                 Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >, \
                                                                                 Antioch::ConstantLewisDiffusivity<libMesh::Real> > >; \
  /*
   *MixtureAveraged: StatMech thermo, Blottner viscosity, Eucken conductivity, constant Lewis diffusivity
   */ \
  template class GRINS::class_name<GRINS::AntiochMixtureAveragedTransportMixture<Antioch::StatMechThermodynamics<libMesh::Real>, \
                                                                                 Antioch::BlottnerViscosity<libMesh::Real>, \
                                                                                 Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >, \
                                                                                 Antioch::ConstantLewisDiffusivity<libMesh::Real> >, \
                                   GRINS::AntiochMixtureAveragedTransportEvaluator<Antioch::StatMechThermodynamics<libMesh::Real>, \
                                                                                   Antioch::BlottnerViscosity<libMesh::Real>, \
                                                                                   Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >, \
                                                                                   Antioch::ConstantLewisDiffusivity<libMesh::Real> > >; \
  /*
   *MixtureAveraged: StatMech thermo, Kinetic theory viscosity, Kinetic theory  conductivity, Molecular binary diffusivity
   */ \
  template class GRINS::class_name<GRINS::AntiochMixtureAveragedTransportMixture<Antioch::StatMechThermodynamics<libMesh::Real>, \
                                                                                 Antioch::KineticsTheoryViscosity<libMesh::Real,Antioch::GSLSpliner>, \
                                                                                 Antioch::KineticsTheoryThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real>,libMesh::Real>, \
                                                                                 Antioch::MolecularBinaryDiffusion<libMesh::Real,Antioch::GSLSpliner> >, \
                                   GRINS::AntiochMixtureAveragedTransportEvaluator<Antioch::StatMechThermodynamics<libMesh::Real>, \
                                                                                   Antioch::KineticsTheoryViscosity<libMesh::Real,Antioch::GSLSpliner>, \
                                                                                   Antioch::KineticsTheoryThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real>,libMesh::Real>, \
                                                                                   Antioch::MolecularBinaryDiffusion<libMesh::Real,Antioch::GSLSpliner> > >; \
  /*
   *MixtureAveraged: IdealGasMicroThermo/CEA thermo, Sutherland viscosity, Eucken conductivity, constant Lewis diffusivity
   */ \
  template class GRINS::class_name<GRINS::AntiochMixtureAveragedTransportMixture<Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real, Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real >, \
                                                                                 Antioch::SutherlandViscosity<libMesh::Real>, \
                                                                                 Antioch::EuckenThermalConductivity<Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real, Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real > >, \
                                                                                 Antioch::ConstantLewisDiffusivity<libMesh::Real> >, \
                                   GRINS::AntiochMixtureAveragedTransportEvaluator<Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real, Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real >, \
                                                                                   Antioch::SutherlandViscosity<libMesh::Real>, \
                                                                                   Antioch::EuckenThermalConductivity<Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real, Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real > >, \
                                                                                   Antioch::ConstantLewisDiffusivity<libMesh::Real> > >; \
  /*
   *MixtureAveraged: IdealGasMicroThermo/CEA thermo, Blottner viscosity, Eucken conductivity, constant Lewis diffusivity
   */ \
  template class GRINS::class_name<GRINS::AntiochMixtureAveragedTransportMixture<Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real, Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real >, \
                                                                                 Antioch::BlottnerViscosity<libMesh::Real>, \
                                                                                 Antioch::EuckenThermalConductivity<Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real, Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real > >, \
                                                                                 Antioch::ConstantLewisDiffusivity<libMesh::Real> >, \
                                   GRINS::AntiochMixtureAveragedTransportEvaluator<Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real, Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real >, \
                                                                                   Antioch::BlottnerViscosity<libMesh::Real>, \
                                                                                   Antioch::EuckenThermalConductivity<Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real, Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real > >, \
                                                                                   Antioch::ConstantLewisDiffusivity<libMesh::Real> > >; \
  /*
   *MixtureAveraged: IdealGasMicroThermo/CEA thermo, Kinetic theory viscosity, Kinetic theory  conductivity, Molecular binary diffusivity
   */ \
  template class GRINS::class_name<GRINS::AntiochMixtureAveragedTransportMixture<Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real, Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real >, \
                                                                                 Antioch::KineticsTheoryViscosity<libMesh::Real,Antioch::GSLSpliner>, \
                                                                                 Antioch::KineticsTheoryThermalConductivity<Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real, Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real >,libMesh::Real>, \
                                                                                 Antioch::MolecularBinaryDiffusion<libMesh::Real,Antioch::GSLSpliner> >, \
                                   GRINS::AntiochMixtureAveragedTransportEvaluator<Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real, Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real >, \
                                                                                   Antioch::KineticsTheoryViscosity<libMesh::Real,Antioch::GSLSpliner>, \
                                                                                   Antioch::KineticsTheoryThermalConductivity<Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real, Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real >,libMesh::Real>, \
                                                                                   Antioch::MolecularBinaryDiffusion<libMesh::Real,Antioch::GSLSpliner> > >; \
  /*
   *Constant viscosity, constant Lewis diffusivity, constant conductivity, StatMech thermo
   */ \
  template class GRINS::class_name<GRINS::AntiochConstantTransportMixture<GRINS::ConstantConductivity>, \
                                   GRINS::AntiochConstantTransportEvaluator<Antioch::StatMechThermodynamics<libMesh::Real>, GRINS::ConstantConductivity> >; \
  /*
   *Constant viscosity, constant Lewis diffusivity, constant Prandtl conductivity, StatMech thermo
   */ \
  template class GRINS::class_name<GRINS::AntiochConstantTransportMixture<GRINS::ConstantPrandtlConductivity>,\
                                   GRINS::AntiochConstantTransportEvaluator<Antioch::StatMechThermodynamics<libMesh::Real>, GRINS::ConstantPrandtlConductivity> >; \
  /*
   *Constant viscosity, constant Lewis diffusivity, constant conductivity, CEA thermo
   */ \
  template class GRINS::class_name<GRINS::AntiochConstantTransportMixture<GRINS::ConstantConductivity>,\
                                   GRINS::AntiochConstantTransportEvaluator<Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real,Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real>, GRINS::ConstantConductivity> >; \
  /*
   *Constant viscosity, constant Lewis diffusivity, constant Prandtl conductivity, CEA thermo
   */ \
  template class GRINS::class_name<GRINS::AntiochConstantTransportMixture<GRINS::ConstantPrandtlConductivity>, \
                                   GRINS::AntiochConstantTransportEvaluator<Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real,Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real>, GRINS::ConstantPrandtlConductivity> >


#define INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_MIXTURE_ONLY(class_name) \
  /*
   *MixtureAveraged: StatMech thermo, Sutherland viscosity, Eucken conductivity, constant Lewis diffusivity
   */ \
  template class GRINS::class_name<GRINS::AntiochMixtureAveragedTransportMixture<Antioch::StatMechThermodynamics<libMesh::Real>, \
                                                                                 Antioch::SutherlandViscosity<libMesh::Real>, \
                                                                                 Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >, \
                                                                                 Antioch::ConstantLewisDiffusivity<libMesh::Real> > >;\
  /*
   *MixtureAveraged: StatMech thermo, Blottner viscosity, Eucken conductivity, constant Lewis diffusivity
   */ \
  template class GRINS::class_name<GRINS::AntiochMixtureAveragedTransportMixture<Antioch::StatMechThermodynamics<libMesh::Real>, \
                                                                                 Antioch::BlottnerViscosity<libMesh::Real>, \
                                                                                 Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >, \
                                                                                 Antioch::ConstantLewisDiffusivity<libMesh::Real> > >; \
  /*
   *MixtureAveraged: StatMech thermo, Kinetic theory viscosity, Kinetic theory  conductivity, Molecular binary diffusivity
   */ \
  template class GRINS::class_name<GRINS::AntiochMixtureAveragedTransportMixture<Antioch::StatMechThermodynamics<libMesh::Real>, \
                                                                                 Antioch::KineticsTheoryViscosity<libMesh::Real,Antioch::GSLSpliner>, \
                                                                                 Antioch::KineticsTheoryThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real>,libMesh::Real>, \
                                                                                 Antioch::MolecularBinaryDiffusion<libMesh::Real,Antioch::GSLSpliner> > >; \
  /*
   *MixtureAveraged: IdealGasMicroThermo/CEA thermo, Sutherland viscosity, Eucken conductivity, constant Lewis diffusivity
   */ \
  template class GRINS::class_name<GRINS::AntiochMixtureAveragedTransportMixture<Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real, Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real >, \
                                                                                 Antioch::SutherlandViscosity<libMesh::Real>, \
                                                                                 Antioch::EuckenThermalConductivity<Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real, Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real > >, \
                                                                                 Antioch::ConstantLewisDiffusivity<libMesh::Real> > >;\
  /*
   *MixtureAveraged: IdealGasMicroThermo/CEA thermo, Blottner viscosity, Eucken conductivity, constant Lewis diffusivity
   */ \
  template class GRINS::class_name<GRINS::AntiochMixtureAveragedTransportMixture<Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real, Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real >, \
                                                                                 Antioch::BlottnerViscosity<libMesh::Real>, \
                                                                                 Antioch::EuckenThermalConductivity<Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real, Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real > >, \
                                                                                 Antioch::ConstantLewisDiffusivity<libMesh::Real> > >;\
  /*
   *MixtureAveraged: StatMech thermo, Kinetic theory viscosity, Kinetic theory  conductivity, Molecular binary diffusivity
   */ \
  template class GRINS::class_name<GRINS::AntiochMixtureAveragedTransportMixture<Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real, Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real >, \
                                                                                 Antioch::KineticsTheoryViscosity<libMesh::Real,Antioch::GSLSpliner>, \
                                                                                 Antioch::KineticsTheoryThermalConductivity<Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real, Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real >,libMesh::Real>, \
                                                                                 Antioch::MolecularBinaryDiffusion<libMesh::Real,Antioch::GSLSpliner> > >; \
  /*
   *Constant transport: constant conductivity
   */ \
  template class GRINS::class_name<GRINS::AntiochConstantTransportMixture<GRINS::ConstantConductivity> >; \
  /*
   *Constant transport: constant Prandtl conductivity
   */ \
  template class GRINS::class_name<GRINS::AntiochConstantTransportMixture<GRINS::ConstantPrandtlConductivity> >

#endif // GRINS_REACTING_LOW_MACH_NAVIER_STOKES_MACRO_H
