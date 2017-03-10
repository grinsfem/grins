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
    This class performs a rayfire across a given mesh.
    The user supplies an origin point on the mesh boundary,
    and spherical angle(s) theta (and phi for 3D) to describe the direction of the rayfire.
    As is standard in spherical coordinates, theta is the angle from the x-axis in the xy-plane,
    with counterclockwise being positive. -2&pi; &le; theta &le; +2&pi;.
    Phi is the angle from the xy-plane, with phi>0 being in the positive half of the z-plane.
    -&pi; &le; phi &le; +&pi;

    Starting at the given origin point, this class will walk in a straight line
    across the mesh in the prescribed direction. On each element along the line,
    the class will calculate the point where it enters and leaves that element.
    Knowing the entering and exiting points, an EDGE2 element
    is created with those points as nodes, and is added to an internal 1D mesh.
    This repeats until a main mesh boundary is reached.

    The intent of this class is to use the 1D mesh to integrate quantities along the
    rayfire path. Using the map_to_rayfire_elem() function, a separate integration class
    can access the 1D elements, prescribe a desired quadrature rule, and perform
    integration directly along the line, rather than using the entire 2D/3D elements
    of the main mesh.

    Refinement and coarsening of the rayfire mesh are supported through the reinit() function.
  */
  class RayfireMesh
  {
  public:

    //! 2D Constructor
    /*!
      This is restricted to 2D meshes. The init() function must
      be called to perform the actual rayfire.
      @param origin Origin point (x,y) of the rayfire on mesh boundary
      @param theta Spherical polar angle (in radians)
    */
    RayfireMesh(libMesh::Point& origin, libMesh::Real theta);

    //! 3D Constructor
    /*!
      This is restricted to 3D meshes. The init() function must
      be called to perform the actual rayfire.
      @param origin Origin point (x,y,z) of the rayfire on mesh boundary
      @param theta  Spherical polar angle (in radians)
      @param phi    Spherical azimuthal angle (in radians)
    */
    RayfireMesh(libMesh::Point& origin, libMesh::Real theta, libMesh::Real phi);

    //! Initialization
    /*!
      This function performs the rayfire and assembles the 1D mesh.
      Call this immediately after the constructor.
      Must be called before any refinements of the main mesh.
      @param mesh_base Reference to the main mesh
    */
    void init(const libMesh::MeshBase& mesh_base);

    /*!
      This function takes in an elem_id on the main mesh and returns an elem from the 1D rayfire mesh.
      @param elem_id The ID of the elem on the main mesh
      @return Elem* to 1D elem from the rayfire mesh if the elem_id falls along the rayfire path
      @return NULL if the given elem_id does not correspond to an elem through which the rayfire passes
    */
    const libMesh::Elem* map_to_rayfire_elem(const libMesh::dof_id_type elem_id);

    /*!
      Returns a vector of elem IDs that are currently along the rayfire
    */
    void elem_ids_in_rayfire(std::vector<libMesh::dof_id_type>& id_vector);

    /*!
      Checks for refined and coarsened main mesh elements along the rayfire path.
      They are then passed to refine() and coarsen(), respectively, to update the rayfire mesh.

      Only 1 level of refinement and/or coarsening can be done between reinit() calls.
      Note that this is not limited to doing either refinement or coarsening between reinits.
      Both can be done on different elements, as long as they are only 1 level in either direction.
      @param mesh: reference to main mesh, needed to get Elem* from the stored elem_id's
    */
    void reinit(const libMesh::MeshBase& mesh_base);


  private:
    //! Dimension of the main mesh
    const unsigned int _dim;

    //! Origin point
    libMesh::Point _origin;

    //! Rayfire Spherical polar angle (in radians)
    libMesh::Real   _theta;

    //! Rayfire Spherical azimuthal angle (in radians)
    libMesh::Real   _phi;

    //! Internal 1D mesh of EDGE2 elements
    SharedPtr<libMesh::Mesh> _mesh;

    //! Map of main mesh elem_id to rayfire mesh elem_id
    std::map<libMesh::dof_id_type,libMesh::dof_id_type> _elem_id_map;

    //! User should never call the default constructor
    RayfireMesh();

    //! Private function to get a rayfire elem from main_mesh elem ID
    /*!
      Does not return a const pointer, and is used within this class to simplify
      rayfire elem access.

      Also used by map_to_rayfire_elem(), which adds a const to prevent modification
      outside this class.
    */
    libMesh::Elem* get_rayfire_elem(const libMesh::dof_id_type elem_id);

    //! Calculate the intersection point
    /*!
      @param[out] end_point The intersection point in (x,y) coordinates
      @return Elem* if intersection is found
      @return NULL  if no intersection point found (i.e. start_point is on a boundary)
    */
    const libMesh::Elem* get_next_elem(const libMesh::Elem* cur_elem, libMesh::Point& start_point, libMesh::Point& end_point, bool same_parent = false);

    //! Ensure the calculated intersection point is on the edge_elem and is not the start_point
    bool check_valid_point(libMesh::Point& intersection_point, libMesh::Point& start_point, libMesh::Elem& edge_elem, libMesh::Point& next_point);

    //! Knowing the end_point, get the appropraite next elem along the path
    const libMesh::Elem* get_correct_neighbor(libMesh::Point& end_point, const libMesh::Elem* cur_elem, unsigned int side, bool same_parent);

    //! Ensure the supplied origin is on a boundary of the mesh
    void check_origin_on_boundary(const libMesh::Elem* start_elem);

    //! Get the correct starting elem corresponding to _origin
    /*!
      @return NULL Origin is not on the mesh
      @return Elem* First elem on the rayfire
    */
    const libMesh::Elem* get_start_elem(const libMesh::MeshBase& mesh_base);

    //! Walks a short distance along the rayfire and checks if elem contains that point
    bool rayfire_in_elem(const libMesh::Point& end_point, const libMesh::Elem* elem);

    //! Iterative solver for calculating the intersection point of the rayfire on the side of an elem
    /*!
      The rayfire line can be represented in point-slope form:

      y-y0 = m(x-x0)

      where m is the slope and (x0,y0) is the known initial_point. To find the intersection point,
      we will use a parametric representation of the edge utilizing the shape functions as such:

      X(&xi;) = sum( x<SUB>j</SUB> &phi;<SUB>j</SUB>(&xi;) )

      Y(&xi;) = sum( y<SUB>j</SUB> &phi;<SUB>j</SUB>(&xi;) )

      where x<SUB>j</SUB> and y<SUB>j</SUB> are the node x- and y-coordinates, respectively,
      and &phi;<SUB>j</SUB> is the value of shape function j at reference coordinate &xi;

      The intersection point is then where the x and y coordinates from the point-slope equation
      equal the X and Y coordinates from the parametric equations for the given edge, respectively.
      We can then form a residual by substituting the parametric coordinates into the point slope
      form and bringing all terms to one side of the equation as:

      f = 0 = m(X(&xi;)-x0) - (Y(&xi;)-y0)

      The residual is then a function of a single parameter, namely the reference coordinate &xi;.
      The derivative is rather trivial and can be expressed as

      df = m*dX - dY

      where

      dX = sum( x<SUB>j</SUB> d&phi;<SUB>j</SUB>/d&xi; )

      dY = sum( y<SUB>j</SUB> d&phi;<SUB>j</SUB>/d&xi; )

      @param[out] intersection_point The calculated intersection point (only set if convergence is achieved)
      @return Whether or not the solver converged before hitting the iteration limit
    */
    bool newton_solve_intersection(libMesh::Point& initial_point, const libMesh::Elem* edge_elem, libMesh::Point& intersection_point);

    //! Refinement of a rayfire element whose main mesh counterpart was refined
    void refine(const libMesh::Elem* main_elem, libMesh::Elem* rayfire_elem);

    //! Coarsening of a rayfire element whose main mesh counterpart was coarsened
    void coarsen(const libMesh::Elem* rayfire_elem);

  };

}
#endif //GRINS_RAYFIRE_MESH_H
