#include "grins/constant_conductivity.h"
#include "grins/parsed_conductivity.h"

#ifndef GRINS_HEAT_TRANSFER_MACRO_H
#define GRINS_HEAT_TRANSFER_MACRO_H

#define INSTANTIATE_HEAT_TRANSFER_SUBCLASS(class_name)                  \
  template class GRINS::class_name<GRINS::ConstantConductivity>;        \
  template class GRINS::class_name<GRINS::ParsedConductivity>

#endif // GRINS_HEAT_TRANSFER_MACRO_H
