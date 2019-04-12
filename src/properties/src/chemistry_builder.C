//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2019 Paul T. Bauman, Roy H. Stogner
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

// This class
#include "grins/chemistry_builder.h"

#include "grins_config.h"

// GRINS
#include "grins/antioch_mixture_builder_base.h"

#ifdef GRINS_HAVE_CANTERA
#include "grins/cantera_mixture.h"
#endif

#ifdef GRINS_HAVE_ANTIOCH
#include "grins/antioch_chemistry.h"
#endif

namespace GRINS
{
#ifdef GRINS_HAVE_CANTERA
  template<>
  void ChemistryBuilder::build_chemistry(const GetPot & input,const std::string & material,
                                         std::unique_ptr<CanteraMixture> & chem_ptr )
  {
    chem_ptr.reset( new CanteraMixture(input,material) );
  }
#endif // GRINS_HAVE_CANTERA

#ifdef GRINS_HAVE_ANTIOCH
  template<>
  void ChemistryBuilder::build_chemistry(const GetPot & input,const std::string & material,
                                         std::unique_ptr<AntiochChemistry> & chem_ptr )
  {
    AntiochMixtureBuilderBase builder;

    std::unique_ptr<Antioch::ChemicalMixture<libMesh::Real> > antioch_chem_mix
      = builder.build_chem_mix(input,material);

    chem_ptr.reset( new AntiochChemistry(antioch_chem_mix) );
  }
#endif // GRINS_HAVE_ANTIOCH

} // end namespace GRINS
