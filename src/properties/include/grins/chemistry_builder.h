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


#ifndef GRINS_CHEMISTRY_BUILDER_H
#define GRINS_CHEMISTRY_BUILDER_H

// libMesh
#include "libmesh/auto_ptr.h" // std::unique_ptr
#include "libmesh/getpot.h"

// C++
#include <string>

namespace GRINS
{
  class ChemistryBuilder
  {
  public:
    ChemistryBuilder(){}
    ~ChemistryBuilder(){}

    template<typename ChemistryType>
    void build_chemistry( const GetPot & input,const std::string & material,
                          std::unique_ptr<ChemistryType> & chem_ptr );
  };

} // end namespace GRINS

#endif // GRINS_CHEMISTRY_BUILDER_H
