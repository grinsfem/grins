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


#ifndef GRINS_RAYFIRE_MESH_H
#define GRINS_RAYFIRE_MESH_H

// libMesh
#include "libmesh/quadrature.h"
#include "libmesh/mesh.h"

// GRINS
#include "grins/qoi_base.h"
#include "grins/variable_name_defaults.h"

namespace GRINS
{
    
  //! RayfireMesh
  /*!
  This class creates a 1D mesh of EDGE2 elements to allow
  for integration of functions along a line across the domain
  */
  class RayfireMesh
  {
  public:
    
    //! 2D Constructor
    /*! 
    The init() function must be called after the constructor
    @param origin Origin point (x,y) of the rayfire on mesh boundary
    @param theta Spherical polar angle (in radians)
    */
    RayfireMesh(libMesh::Point& origin, libMesh::Real theta);
    
    //! 3D Constructor
    /*! 
    The init() function must be called after the constructor
    @param origin Origin point (x,y,z) of the rayfire on mesh boundary
    @param theta  Spherical polar angle (in radians)
    @param phi    Spherical azimuthal angle (in radians)
    */
    RayfireMesh(libMesh::Point& origin, libMesh::Real theta, libMesh::Real phi);
    
    //! Initialization
    /*! 
    This function performs the rayfire and assembles the 1D mesh 
    @param mesh_base Reference to the main mesh
    */
    void init(const libMesh::MeshBase& mesh_base);
    
    /*! 
    This function takes in an elem_id on the main mesh and returns an elem from the 1D rayfire mesh 
    @param elem_id The ID of the elem on the main mesh
    */
    const libMesh::Elem* translate(const libMesh::dof_id_type elem_id);


  private:
    //! Dimension of the main mesh 
    const unsigned int _dim;
    
    //! Origin point
    libMesh::Point& _origin;
    
    //! Rayfire Spherical polar angle (in radians)
    libMesh::Real   _theta;
    
    //! Rayfire Spherical azimuthal angle (in radians)
    libMesh::Real   _phi;
    
    //! Internal 1D mesh of EDGE2 elements
    SharedPtr<libMesh::Mesh> _mesh;
    
    //! Map of main mesh elem_id to rayfire mesh elems
    std::map<libMesh::dof_id_type,libMesh::Elem*> _elem_id_map;
    
    
    //! Calculate the intersection point
    /*!
    @param[out] end_point The intersection point in (x,y) coordinates
    @return Elem* if intersection is found
    @return NULL  if no intersection point found (i.e. start_point is on a boundary)
    */
    const libMesh::Elem* _get_next_elem(const libMesh::Elem* cur_elem, libMesh::Point* start_point, libMesh::Point* end_point);
    
    //! Ensure the calculated intersection point is on the edge_elem and is not the start_point
    bool _check_valid_point(libMesh::Point& intersection_point, libMesh::Point& start_point, libMesh::Elem& edge_elem, libMesh::Point* next_point);
    
    //! Knowing the intersection_point, get the appropraite next elem along the path
    const libMesh::Elem* _get_correct_neighbor(libMesh::Point& end_point, const libMesh::Elem* cur_elem, unsigned int side);
    
    //! Ensure the supplied origin is on a boundary of the mesh
    void _check_origin_on_boundary(const libMesh::Elem* start_elem);
    
    //! Iterative solver for calculating the intersection point of the rayfire on the side of an elem
    /*!
    @param[out] intersection_point The calculated intersection point (only set if convergence is achieved)
    @return Whether or not the solver converged before hitting the iteration limit
    */
    bool _newton_solver(libMesh::Point& initial_point, const libMesh::Elem* edge_elem, libMesh::Point* intersection_point);
  };
      
}
#endif //GRINS_RAYFIRE_MESH_H
