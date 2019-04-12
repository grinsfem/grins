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

#ifndef GRINS_BOUNDARY_CONDITION_NAMES_H
#define GRINS_BOUNDARY_CONDITION_NAMES_H

namespace GRINS
{
  //! Class to contain names of sections and types of boundary conditions
  class BoundaryConditionNames
  {
  public:

    //! Outer section name for boundary conditions section in input file
    static std::string bc_section()
    { return "BoundaryConditions"; }

    //! Variable for list boundary ids to correspond with bc_id_name_map
    static std::string bc_ids_var()
    { return BoundaryConditionNames::bc_section()+"/bc_ids"; }

    //! Names of boundaries to correspond with bc_ids
    /*! These will be the subsections of bc_section() that are parsed
      for the boundary condition types and values for each variable */
    static std::string bc_id_name_map_var()
    { return BoundaryConditionNames::bc_section()+"/bc_id_name_map"; }

    //! "Standard" input variable for Neumann flux in input file
    static std::string bc_flux_var()
    { return "normal_flux"; }

    //! Input variable for tractions in input file
    /*! Effectively synonomous with normal_flux, just syntax sugar. */
    static std::string traction_var()
    { return "traction"; }

    //! Default boundary name prefix if bc_id_name_map/bc_ids are not used
    /*! If the user opts to not provide bc_id_name_map/bc_ids, then we'll
      parse for this prefix appened with the boundary id from the mesh. */
    static std::string bc_name_prefix_default()
    { return "Boundary"; }

    static std::string axisymmetric()
    { return "axisymmetric"; }

    static std::string periodic()
    { return "periodic"; }

  };
} // end namespace GRINS

#endif // GRINS_BOUNDARY_CONDITION_NAMES_H
