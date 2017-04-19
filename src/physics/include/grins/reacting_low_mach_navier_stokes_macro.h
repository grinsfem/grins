#ifndef GRINS_REACTING_LOW_MACH_NAVIER_STOKES_MACRO_H
#define GRINS_REACTING_LOW_MACH_NAVIER_STOKES_MACRO_H

#include "grins/antioch_mixture_averaged_transport_mixture.h"
#include "grins/antioch_mixture_averaged_transport_evaluator.h"

#include "grins/antioch_constant_transport_mixture.h"
#include "grins/antioch_constant_transport_evaluator.h"

namespace GRINSPrivate
{
  typedef  Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real,Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real> CEAIdealGasThermo;
}

#define INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_CONSTANT_MIXTURE_RAW(class_name,curve_fit,conductivity) \
  template class GRINS::class_name<GRINS::AntiochConstantTransportMixture<curve_fit,GRINS::conductivity> >

#define INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_CONSTANT_MIXTURE_CURVEFIT_RAW(class_name,conductivity) \
  template class GRINS::class_name<GRINS::AntiochConstantTransportMixture<Antioch::CEACurveFit<libMesh::Real>,GRINS::conductivity> >

#define INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_CONSTANT_MIXTURE(class_name) \
  INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_CONSTANT_MIXTURE_CURVEFIT_RAW(class_name,ConstantConductivity); \
  INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_CONSTANT_MIXTURE_CURVEFIT_RAW(class_name,ConstantPrandtlConductivity)


#define INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_CONSTANT_MIXTURE_AND_CONSTANT_EVALUATOR_RAW(class_name,curve_fit,conductivity,thermo) \
  template class GRINS::class_name<GRINS::AntiochConstantTransportMixture<curve_fit,GRINS::conductivity>, \
                                   GRINS::AntiochConstantTransportEvaluator<curve_fit,thermo,GRINS::conductivity> >

#define INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_CONSTANT_MIXTURE_AND_CONSTANT_EVALUATOR_CURVEFIT_THERMO_RAW(class_name,conductivity) \
  INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_CONSTANT_MIXTURE_AND_CONSTANT_EVALUATOR_RAW(class_name,Antioch::CEACurveFit<libMesh::Real>,conductivity,Antioch::StatMechThermodynamics<libMesh::Real>);\
  INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_CONSTANT_MIXTURE_AND_CONSTANT_EVALUATOR_RAW(class_name,Antioch::CEACurveFit<libMesh::Real>,conductivity,GRINSPrivate::CEAIdealGasThermo)

#define INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_CONSTANT_MIXTURE_AND_CONSTANT_EVALUATOR(class_name) \
  INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_CONSTANT_MIXTURE_AND_CONSTANT_EVALUATOR_CURVEFIT_THERMO_RAW(class_name,ConstantConductivity);\
  INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_CONSTANT_MIXTURE_AND_CONSTANT_EVALUATOR_CURVEFIT_THERMO_RAW(class_name,ConstantPrandtlConductivity)



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
                                                                                   Antioch::MolecularBinaryDiffusion<libMesh::Real,Antioch::GSLSpliner> > >


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
                                                                                 Antioch::MolecularBinaryDiffusion<libMesh::Real,Antioch::GSLSpliner> > >

#endif // GRINS_REACTING_LOW_MACH_NAVIER_STOKES_MACRO_H
