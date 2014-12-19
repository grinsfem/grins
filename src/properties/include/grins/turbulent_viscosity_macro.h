#ifndef GRINS_TURBULENT_VISCOSITY_MACRO_H
#define GRINS_TURBULENT_VISCOSITY_MACRO_H

#define INSTANTIATE_TURBULENT_VISCOSITY_SUBCLASS(class_name) \
template class GRINS::class_name<GRINS::ConstantViscosity>; \
 template class GRINS::class_name<GRINS::ParsedViscosity>; \
template class GRINS::class_name<GRINS::SpalartAllmarasViscosity<GRINS::ConstantViscosity> >

#endif // GRINS_TURBULENT_VISCOSITY_MACRO_H
