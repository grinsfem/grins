//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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
#include "grins/constrained_points.h"

// GRINS
#include "grins/builder_helper.h"
#include "grins/variables_parsing.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/system.h"

namespace GRINS
{
  ConstrainedPoints::ConstrainedPoints( const GetPot& input,
                                        libMesh::System& system )
    : _sys(system)
  {
    std::vector<std::string> var_sections;

    BuilderHelper::parse_var_sections_vector(input, var_sections);

    for (std::vector<std::string>::const_iterator
           it = var_sections.begin(); it != var_sections.end(); ++it)
      {
        const std::string & section_name = *it;
        const std::string constraint_section =
          VariablesParsing::variables_section() + '/' +
          section_name + '/' + "ConstrainedPoints";

        if (input.have_section (constraint_section))
          {
            std::vector<std::string> point_sections =
              input.get_subsection_names(constraint_section);

            for (std::vector<std::string>::const_iterator
                   pt_it = point_sections.begin();
                   pt_it != point_sections.end(); ++pt_it)
              {
                ConstrainedPoint pt;
                
                pt.name = *pt_it;

                const std::string point_section =
                  constraint_section + '/' + pt.name;

                const unsigned int constrained_loc_dim =
                  input.vector_variable_size(point_section+"/constraint_location");
                libmesh_assert_less_equal(constrained_loc_dim, LIBMESH_DIM);
                for (unsigned int d=0; d != constrained_loc_dim; ++d)
                  pt(d) = input(point_section+"/constraint_location", 0.0, d );

                const unsigned int n_constraining_points =
                  input.vector_variable_size(point_section+"/constraining_points_x");

#ifndef NDEBUG
                const unsigned int n_constraining_points_y =
                  input.vector_variable_size(point_section+"/constraining_points_y");
                if (n_constraining_points_y)
                  libmesh_assert_equal_to(n_constraining_points_x, n_constraining_points_y);

                const unsigned int n_constraining_points_z =
                  input.vector_variable_size(point_section+"/constraining_points_z");
                if (n_constraining_points_z)
                  libmesh_assert_equal_to(n_constraining_points_x, n_constraining_points_z);

                const unsigned int n_constraining_points_var =
                  input.vector_variable_size(point_section+"/constraining_points_var");
                libmesh_assert_equal_to(n_constraining_points_x, n_constraining_points_var);

                const unsigned int n_constraining_points_coeff =
                  input.vector_variable_size(point_section+"/constraining_points_coeff");
                if (n_constraining_points_coeff)
                  libmesh_assert_equal_to(n_constraining_points_x, n_constraining_points_coeff);
#endif // !NDEBUG

                for (unsigned int n=0; n != n_constraining_points; ++n)
                  {
                    ConstrainingPoint constraining_pt;
                    constraining_pt(0) = input(point_section+"/constraining_points_x", 0.0, n);
                    if (LIBMESH_DIM > 1)
                      constraining_pt(1) = input(point_section+"/constraining_points_y", 0.0, n);
                    if (LIBMESH_DIM > 2)
                      constraining_pt(2) = input(point_section+"/constraining_points_z", 0.0, n);

                    constraining_pt.coeff =
                      input(point_section+"/constraining_points_coeff", libMesh::Number(0), n);

                    const std::string var_name
                      (input(point_section+"/constraining_points_var", std::string(), n));

                    constraining_pt.var = _sys.variable_number(var_name);

                    pt.constrainers.push_back(constraining_pt);
                  }

                _constrained_pts.push_back(pt);
              }
          }
      }
  }

  void ConstrainedPoints::constrain ()
  {
  }
} // end namespace GRINS
