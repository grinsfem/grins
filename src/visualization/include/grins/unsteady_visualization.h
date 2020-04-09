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


#ifndef GRINS_UNSTEADY_VISUALIZATION_H
#define GRINS_UNSTEADY_VISUALIZATION_H

// GRINS
#include "grins/visualization.h"

namespace GRINS
{
  class UnsteadyVisualization : public Visualization
  {
  public:

    UnsteadyVisualization(const GetPot& input,
                          const libMesh::Parallel::Communicator &comm );
    ~UnsteadyVisualization();

    virtual void output_residual( std::shared_ptr<libMesh::EquationSystems> equation_system,
                                  MultiphysicsSystem* system,
                                  const unsigned int time_step,
                                  const libMesh::Real time);

    virtual void output_residual_sensitivities
    (std::shared_ptr<libMesh::EquationSystems> equation_system,
     MultiphysicsSystem* system,
     const libMesh::ParameterVector & params,
     const unsigned int time_step,
     const libMesh::Real time);

    virtual void output_adjoint( std::shared_ptr<libMesh::EquationSystems> equation_system,
                                 MultiphysicsSystem* system,
                                 const unsigned int time_step,
                                 const libMesh::Real time );

    virtual void output_solution_sensitivities
    (std::shared_ptr<libMesh::EquationSystems> equation_system,
     MultiphysicsSystem* system,
     const libMesh::ParameterVector & params,
     const unsigned int time_step,
     const libMesh::Real time);
  };

} // end namespace GRINS
#endif // GRINS_UNSTEADY_VISUALIZATION_H
