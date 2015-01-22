#ifndef GRINS_INC_NAV_STOKES_MACRO_H
#define GRINS_INC_NAV_STOKES_MACRO_H

#define INSTANTIATE_INC_NS_SUBCLASS(class_name) \
template class GRINS::class_name<GRINS::ConstantViscosity>; \
template class GRINS::class_name<GRINS::ParsedViscosity>

#endif // GRINS_INC_NAV_STOKES_MACRO_H
