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

#ifndef GRINS_CONSTRAINED_POINTS_H
#define GRINS_CONSTRAINED_POINTS_H

// libMesh
#include "libmesh/system.h"

// libMesh forward declarations
class GetPot;

namespace GRINS
{
  //! Manages construction of Constraint application object.
  class ConstrainedPoints : public libMesh::System::Constraint
  {
  public:

    ConstrainedPoints ( const GetPot& input,
                        libMesh::System& system );

    virtual ~ConstrainedPoints(){}

    virtual void constrain ();

  private:

    libMesh::System & _sys;

    struct ConstrainingPoint : public libMesh::Point {
      libMesh::Number coeff;
      unsigned short var;
    };

    struct ConstrainedPoint : public libMesh::Point {
      std::string name;
      std::vector<ConstrainingPoint> constrainers;
      libMesh::Number rhs;
      unsigned short var;
      bool forbid_overwrite;
    };

    std::vector<ConstrainedPoint> _constrained_pts;
  };
} // end namespace GRINS

#endif // GRINS_CONSTRAINED_POINTS_H
