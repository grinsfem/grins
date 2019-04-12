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
    As is standard in spherical coordinates, theta is the azimuthal angle from
    the positive x-axis in the xy-plane, with counterclockwise being positive;
    \f$-2\pi \leq \theta \leq +2\pi\f$.
    Phi is the polar angle from the positive z-axis, with \f$\phi = \pi/2\f$ being the xy-plane;
    \f$ 0 \leq \phi \leq \pi \f$.

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
      @param theta Spherical azimuthal angle (in radians)
    */
    RayfireMesh(libMesh::Point & origin, libMesh::Real theta);

    //! 3D Constructor
    /*!
      This is restricted to 3D meshes. The init() function must
      be called to perform the actual rayfire.
      @param origin Origin point (x,y,z) of the rayfire on mesh boundary
      @param theta  Spherical azimuthal angle (in radians)
      @param phi    Spherical polar angle (in radians)
    */
    RayfireMesh(libMesh::Point & origin, libMesh::Real theta, libMesh::Real phi);

    //! Input File Constructor
    /*!
      Creates rayfire from parameters specified in input file section: QoI/'qoi_string'/Rayfire/
    */
    RayfireMesh(const GetPot & input, const std::string & qoi_string);

    //! Enable rayfire output after init()
    void enable_rayfire_output(const std::string & filename)
    {
      _output_filename = filename;
    }

    //! Initialization
    /*!
      This function performs the rayfire and assembles the 1D mesh.
      Call this immediately after the constructor.
      Must be called before any refinements of the main mesh.
      @param mesh_base Reference to the main mesh
    */
    void init(const libMesh::MeshBase & mesh_base);

    /*!
      This function takes in an elem_id on the main mesh and returns an elem from the 1D rayfire mesh.
      @param elem_id The ID of the elem on the main mesh
      @return Elem* to 1D elem from the rayfire mesh if the elem_id falls along the rayfire path
      @return NULL if the given elem_id does not correspond to an elem through which the rayfire passes
    */
    const libMesh::Elem * map_to_rayfire_elem(const libMesh::dof_id_type elem_id);

    /*!
      Returns a vector of elem IDs that are currently along the rayfire
    */
    void elem_ids_in_rayfire(std::vector<libMesh::dof_id_type> & id_vector) const;

    /*!
      Checks for refined and coarsened main mesh elements along the rayfire path.
      They are then passed to refine() and coarsen(), respectively, to update the rayfire mesh.

      Only 1 level of refinement and/or coarsening can be done between reinit() calls.
      Note that this is not limited to doing either refinement or coarsening between reinits.
      Both can be done on different elements, as long as they are only 1 level in either direction.
      @param mesh: reference to main mesh, needed to get Elem* from the stored elem_id's
    */
    void reinit(const libMesh::MeshBase & mesh_base);


  private:
    //! Dimension of the main mesh
    const unsigned int _dim;

    //! Origin point
    libMesh::Point _origin;

    //! Spherical azimuthal angle (in radians)
    libMesh::Real _theta;

    //! Spherical polar angle (in radians)
    libMesh::Real _phi;

    //! Internal 1D mesh of EDGE2 elements
    std::unique_ptr<libMesh::Mesh> _mesh;

    //! Map of main mesh elem_id to rayfire mesh elem_id
    std::map<libMesh::dof_id_type,libMesh::dof_id_type> _elem_id_map;

    //! Filename to output rayfire after init()
    //! Defaults to empty string (no output)
    std::string _output_filename;

    //! User should never call the default constructor
    RayfireMesh();

    void validate_rayfire_angles();

    //! Private function to get a rayfire elem from main_mesh elem ID
    /*!
      Does not return a const pointer, and is used within this class to simplify
      rayfire elem access.

      Also used by map_to_rayfire_elem(), which adds a const to prevent modification
      outside this class.
    */
    libMesh::Elem * get_rayfire_elem(const libMesh::dof_id_type elem_id);

    //! Calculate the intersection point
    /*!
      @param[out] end_point The intersection point in (x,y) coordinates
      @return Elem* if intersection is found
      @return NULL  if no intersection point found (i.e. start_point is on a boundary)
    */
    const libMesh::Elem * get_next_elem(const libMesh::Elem * cur_elem, libMesh::Point & start_point, libMesh::Point & end_point, bool same_parent = false);

    //! Ensure the calculated intersection point is on the edge_elem and is not the start_point
    bool check_valid_point(libMesh::Point & intersection_point, libMesh::Point & start_point, const libMesh::Elem & edge_elem, libMesh::Point & next_point);

    //! Knowing the end_point, get the appropraite next elem along the path
    const libMesh::Elem * get_correct_neighbor(libMesh::Point & start_point, libMesh::Point & end_point, const libMesh::Elem * cur_elem, unsigned int side, bool same_parent);

    //! Ensure the supplied origin is on a boundary of the mesh
    void check_origin_on_boundary(const libMesh::Elem* start_elem);

    //! Get the correct starting elem corresponding to _origin
    /*!
      @return NULL Origin is not on the mesh
      @return Elem* First elem on the rayfire
    */
    const libMesh::Elem * get_start_elem(const libMesh::MeshBase& mesh_base);

    //! Ensures the rayfire doesn't wander into the middle of an elem
    bool validate_edge(const libMesh::Point & start_point, const libMesh::Point & end_point, const libMesh::Elem * side_elem, const libMesh::Elem * neighbor);

    //! Walks a short distance along the rayfire and checks if elem contains that point
    bool rayfire_in_elem(const libMesh::Point & end_point, const libMesh::Elem* elem);

    //! Find the intersection point of the rayfire and a side of cur_elem. Returns libMesh::invalid_uint if no intersection is found
    unsigned int calculate_intersection_point(libMesh::Point & initial_point, const libMesh::Elem * cur_elem, libMesh::Point & intersection_point);

    //! Edges of 2D FIRST order elements can be represented in point-slope form, so there is an analytical solution for the intersection point
    unsigned int intersection_2D_first_order(libMesh::Point & initial_point, const libMesh::Elem * cur_elem, libMesh::Point & intersection_point, unsigned int initial_side = libMesh::invalid_uint);

    //! Edges of 2D SECOND order elements can be nonlinear, so use a standard Newton iterative solver to find the intersection point
    /*!
      Knowing the starting point \f$(x_r,y_r)\f$ of the rayfire on the current element,
      the rayfire line can be parameterized as

      \f$ x = x_r + L*cos(\theta) \f$

      \f$ y = y_r + L*sin(\theta) \f$

      We can also represent a point \f$(X,Y)\f$ on a given edge of the element as a function of
      the reference coordinate \f$\xi\f$ along that edge

      \f$ X = \sum_i x_i \phi_i(\xi) \f$

      \f$ Y = \sum_i y_i \phi_i(\xi) \f$

      At the intersection point of the rayfire line and the edge,
      \f$x=X\f$ and \f$y=Y\f$ so we can form a residual vector

      \f$ F = F(L,\xi) = \begin{bmatrix}
                            x_r + L*cos(\theta) - \sum_i x_i \phi_i(\xi) \\
                            y_r + L*sin(\theta) - \sum_i y_i \phi_i(\xi)
                          \end{bmatrix}
            = 0
      \f$

      with Jacobian \f$J\f$

      \f$ J = \begin{bmatrix}
                cos(\theta) & - \sum_i x_i \frac{d \phi_i}{d\xi} \\
                sin(\theta) & - \sum_i y_i \frac{d \phi_i}{d\xi}
              \end{bmatrix}
      \f$

      Then we can use a Newton iteration
      \f$ \begin{bmatrix}
                \Delta L \\
                \Delta \xi
              \end{bmatrix}
          = \delta = - J^{-1} F\f$
      and converge when \f$|| \delta || <\f$ libMesh::TOLERANCE
    */
    unsigned int intersection_2D_second_order(libMesh::Point & initial_point, const libMesh::Elem * cur_elem, libMesh::Point & intersection_point);

    //! Faces of 3D elements (of FIRST or SECOND) can be nonlinear, so use a standard Newton iterative solver to find the intersection point
    /*!
      Knowing the starting point \f$(x_r,y_r,z_r)\f$ of the rayfire on the current element,
      the rayfire line can be parameterized as

      \f$ x = x_r + L*cos(\theta)*sin(\phi) \f$

      \f$ y = y_r + L*sin(\theta)*sin(phi) \f$
      
      \f$ z = z_r + L*cos(\phi) \f$

      We can also represent a point \f$(X,Y,Z)\f$ on a given face of the element as a function of
      the reference coordinates \f$\xi,\eta\f$ on that face

      \f$ X = \sum_i x_i \phi_i(\xi,\eta) \f$

      \f$ Y = \sum_i y_i \phi_i(\xi,\eta) \f$
      
      \f$ Z = \sum_i z_i \phi_i(\xi,\eta) \f$

      At the intersection point of the rayfire line and the edge,
      \f$x=X\f$, \f$y=Y\f$, and \f$z=Z\f$ so we can form a residual vector

      \f$ F = F(L,\xi,\eta) = \begin{bmatrix}
                                x_r + L*cos(\theta)*sin(\phi) - \sum_i x_i \phi_i(\xi) \\
                                y_r + L*sin(\theta)*sin(\phi) - \sum_i y_i \phi_i(\xi) \\
                                z_r + L*cos(\phi) - \sum_i z_i \phi_i(\xi)
                              \end{bmatrix}
            = 0
      \f$

      with Jacobian \f$J\f$

      \f$ J = \begin{bmatrix}
                cos(\theta)*sin(\phi) & - \sum_i x_i \frac{d \phi_i}{d\xi} & - \sum_i x_i \frac{d \phi_i}{d\eta} \\
                sin(\theta)*sin(\phi) & - \sum_i y_i \frac{d \phi_i}{d\xi} & - \sum_i y_i \frac{d \phi_i}{d\eta} \\
                cos(\phi)             & - \sum_i z_i \frac{d \phi_i}{d\xi} & - \sum_i z_i \frac{d \phi_i}{d\eta}
              \end{bmatrix}
      \f$

      Then we can use a Newton iteration
      \f$ \begin{bmatrix}
                \Delta L \\
                \Delta \xi \\
                \Delta \eta
              \end{bmatrix}
          = \delta = - J^{-1} F\f$
      and converge when \f$|| \delta || <\f$ libMesh::TOLERANCE
    */
    unsigned int intersection_3D(libMesh::Point & initial_point, const libMesh::Elem * cur_elem, libMesh::Point & intersection_point);

    //! Create a FIRST order version of a SECOND order element
    libMesh::Elem * copy_to_first_order(const libMesh::Elem * elem);

    //! Refinement of a rayfire element whose main mesh counterpart was refined
    void refine(const libMesh::Elem * main_elem, libMesh::Elem * rayfire_elem);

    //! Coarsening of a rayfire element whose main mesh counterpart was coarsened
    void coarsen(const libMesh::Elem * rayfire_elem);

    //! Helper function that solves a 3x3 system
    //! @return false if the matrix A is singular (i.e. there is no solution)
    bool system_solve_3x3(libMesh::DenseMatrix<libMesh::Real> & A, libMesh::DenseVector<libMesh::Real> & b, libMesh::DenseVector<libMesh::Real> & x);

  };

}
#endif //GRINS_RAYFIRE_MESH_H
