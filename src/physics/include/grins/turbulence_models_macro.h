#ifndef GRINS_TURBULENCE_MODELS_MACRO_H
#define GRINS_TURBULENCE_MODELS_MACRO_H

#define INSTANTIATE_TURBULENCE_MODELS_SUBCLASS(class_name) \
template class GRINS::class_name<GRINS::ConstantViscosity>; \
template class GRINS::class_name<GRINS::ParsedViscosity>

#endif // GRINS_TURBULENCE_MODELS_MACRO_H
