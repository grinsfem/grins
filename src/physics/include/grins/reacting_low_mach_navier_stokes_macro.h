#ifndef GRINS_REACTING_LOW_MACH_NAVIER_STOKES_MACRO_H
#define GRINS_REACTING_LOW_MACH_NAVIER_STOKES_MACRO_H

#include "grins/antioch_mixture_averaged_transport_mixture.h"
#include "grins/antioch_mixture_averaged_transport_evaluator.h"

#include "grins/antioch_constant_transport_mixture.h"
#include "grins/antioch_constant_transport_evaluator.h"

namespace GRINSPrivate
{
  // Need typedefs for these because the commas in the template arguments screw up the C preprocessor
  // when putting the full types in the argument list
  typedef Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real,Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real> CEAIdealGasThermo;
  typedef Antioch::KineticsTheoryViscosity<libMesh::Real,Antioch::GSLSpliner> KineticsViscosity;
  typedef Antioch::MolecularBinaryDiffusion<libMesh::Real,Antioch::GSLSpliner> BinaryDiffusion;
  typedef Antioch::KineticsTheoryThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real>,libMesh::Real> KineticsConductivityStatMech;
  typedef Antioch::KineticsTheoryThermalConductivity<CEAIdealGasThermo,libMesh::Real> KineticsConductivityCEA;
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
  INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_CONSTANT_MIXTURE_AND_CONSTANT_EVALUATOR_RAW(class_name,Antioch::CEACurveFit<libMesh::Real>,conductivity,Antioch::StatMechThermodynamics<libMesh::Real>); \
  INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_CONSTANT_MIXTURE_AND_CONSTANT_EVALUATOR_RAW(class_name,Antioch::CEACurveFit<libMesh::Real>,conductivity,GRINSPrivate::CEAIdealGasThermo)

#define INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_CONSTANT_MIXTURE_AND_CONSTANT_EVALUATOR(class_name) \
  INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_CONSTANT_MIXTURE_AND_CONSTANT_EVALUATOR_CURVEFIT_THERMO_RAW(class_name,ConstantConductivity); \
  INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_CONSTANT_MIXTURE_AND_CONSTANT_EVALUATOR_CURVEFIT_THERMO_RAW(class_name,ConstantPrandtlConductivity)



#define INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_MIXTUREAVERAGED_MIXTURE_ONLY_RAW(class_name,curve_fit,thermo,viscosity,conductivity,diffusivity) \
  template class GRINS::class_name<GRINS::AntiochMixtureAveragedTransportMixture<curve_fit,thermo,viscosity,conductivity,diffusivity> >

#define INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_MIXTUREAVERAGED_MIXTURE_ONLY_CURVEFIT_THERMO_CONDUCTIVITY_CONSTLEWIS_RAW(class_name,viscosity) \
  INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_MIXTUREAVERAGED_MIXTURE_ONLY_RAW(class_name, \
                                                                          Antioch::CEACurveFit<libMesh::Real>, \
                                                                          Antioch::StatMechThermodynamics<libMesh::Real>, \
                                                                          viscosity, \
                                                                          Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >, \
                                                                          Antioch::ConstantLewisDiffusivity<libMesh::Real>); \
  INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_MIXTUREAVERAGED_MIXTURE_ONLY_RAW(class_name, \
                                                                          Antioch::CEACurveFit<libMesh::Real>, \
                                                                          GRINSPrivate::CEAIdealGasThermo, \
                                                                          viscosity, \
                                                                          Antioch::EuckenThermalConductivity<GRINSPrivate::CEAIdealGasThermo>, \
                                                                          Antioch::ConstantLewisDiffusivity<libMesh::Real>)

#define INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_MIXTUREAVERAGED_MIXTURE_ONLY_KINETICS_THEORY_RAW(class_name) \
  INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_MIXTUREAVERAGED_MIXTURE_ONLY_RAW(class_name, \
                                                                          Antioch::CEACurveFit<libMesh::Real>, \
                                                                          Antioch::StatMechThermodynamics<libMesh::Real>, \
                                                                          GRINSPrivate::KineticsViscosity, \
                                                                          GRINSPrivate::KineticsConductivityStatMech, \
                                                                          GRINSPrivate::BinaryDiffusion); \
  INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_MIXTUREAVERAGED_MIXTURE_ONLY_RAW(class_name, \
                                                                          Antioch::CEACurveFit<libMesh::Real>, \
                                                                          GRINSPrivate::CEAIdealGasThermo, \
                                                                          GRINSPrivate::KineticsViscosity, \
                                                                          GRINSPrivate::KineticsConductivityCEA, \
                                                                          GRINSPrivate::BinaryDiffusion)

#define INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_MIXTURE_ONLY(class_name) \
  INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_MIXTUREAVERAGED_MIXTURE_ONLY_CURVEFIT_THERMO_CONDUCTIVITY_CONSTLEWIS_RAW(class_name,Antioch::SutherlandViscosity<libMesh::Real>); \
  INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_MIXTUREAVERAGED_MIXTURE_ONLY_CURVEFIT_THERMO_CONDUCTIVITY_CONSTLEWIS_RAW(class_name,Antioch::BlottnerViscosity<libMesh::Real>); \
  INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_MIXTUREAVERAGED_MIXTURE_ONLY_KINETICS_THEORY_RAW(class_name)



#define INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_MIXTUREAVERAGED_MIXTURE_AND_EVALUATOR_RAW(class_name,curve_fit,thermo,viscosity,conductivity,diffusivity) \
  template class GRINS::class_name<GRINS::AntiochMixtureAveragedTransportMixture<curve_fit,thermo,viscosity,conductivity,diffusivity>, \
                                   GRINS::AntiochMixtureAveragedTransportEvaluator<curve_fit,thermo,viscosity,conductivity,diffusivity> >

#define INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_MIXTUREAVERAGED_MIXTURE_AND_EVALUATOR_CURVEFIT_THERMO_CONDUCTIVITY_CONSTLEWIS_RAW(class_name,viscosity) \
  INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_MIXTUREAVERAGED_MIXTURE_AND_EVALUATOR_RAW(class_name, \
                                                                                   Antioch::CEACurveFit<libMesh::Real>, \
                                                                                   Antioch::StatMechThermodynamics<libMesh::Real>, \
                                                                                   viscosity, \
                                                                                   Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >, \
                                                                                   Antioch::ConstantLewisDiffusivity<libMesh::Real>); \
  INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_MIXTUREAVERAGED_MIXTURE_AND_EVALUATOR_RAW(class_name, \
                                                                                   Antioch::CEACurveFit<libMesh::Real>, \
                                                                                   GRINSPrivate::CEAIdealGasThermo, \
                                                                                   viscosity, \
                                                                                   Antioch::EuckenThermalConductivity<GRINSPrivate::CEAIdealGasThermo>, \
                                                                                   Antioch::ConstantLewisDiffusivity<libMesh::Real>)

#define INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_MIXTUREAVERAGED_MIXTURE_AND_EVALUATOR_KINETICS_THEORY_RAW(class_name) \
  INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_MIXTUREAVERAGED_MIXTURE_AND_EVALUATOR_RAW(class_name, \
                                                                                   Antioch::CEACurveFit<libMesh::Real>, \
                                                                                   Antioch::StatMechThermodynamics<libMesh::Real>, \
                                                                                   GRINSPrivate::KineticsViscosity, \
                                                                                   GRINSPrivate::KineticsConductivityStatMech, \
                                                                                   GRINSPrivate::BinaryDiffusion); \
  INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_MIXTUREAVERAGED_MIXTURE_AND_EVALUATOR_RAW(class_name, \
                                                                                   Antioch::CEACurveFit<libMesh::Real>, \
                                                                                   GRINSPrivate::CEAIdealGasThermo, \
                                                                                   GRINSPrivate::KineticsViscosity, \
                                                                                   GRINSPrivate::KineticsConductivityCEA, \
                                                                                   GRINSPrivate::BinaryDiffusion)

#define INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_MIXTURE_AND_EVALUATOR(class_name) \
  INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_MIXTUREAVERAGED_MIXTURE_AND_EVALUATOR_CURVEFIT_THERMO_CONDUCTIVITY_CONSTLEWIS_RAW(class_name,Antioch::SutherlandViscosity<libMesh::Real>); \
  INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_MIXTUREAVERAGED_MIXTURE_AND_EVALUATOR_CURVEFIT_THERMO_CONDUCTIVITY_CONSTLEWIS_RAW(class_name,Antioch::BlottnerViscosity<libMesh::Real>); \
  INSTANTIATE_REACTING_LOW_MACH_SUBCLASS_MIXTUREAVERAGED_MIXTURE_AND_EVALUATOR_KINETICS_THEORY_RAW(class_name)

#endif // GRINS_REACTING_LOW_MACH_NAVIER_STOKES_MACRO_H
