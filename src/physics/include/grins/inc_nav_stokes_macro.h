#ifndef GRINS_INC_NAV_STOKES_MACRO_H
#define GRINS_INC_NAV_STOKES_MACRO_H

#include "grins/constant_viscosity.h"
#include "grins/parsed_viscosity.h"
#include "grins/spalart_allmaras_viscosity.h"

#define INSTANTIATE_INC_NS_SUBCLASS(class_name)                         \
  template class GRINS::class_name<GRINS::ConstantViscosity>;           \
  template class GRINS::class_name<GRINS::ParsedViscosity>;             \
  template class GRINS::class_name<GRINS::SpalartAllmarasViscosity<GRINS::ConstantViscosity> >

#endif // GRINS_INC_NAV_STOKES_MACRO_H
