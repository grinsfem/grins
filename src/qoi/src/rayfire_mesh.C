//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
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
#include "grins/rayfire_mesh.h"

// GRINS
#include "grins/math_constants.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/elem.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/fe.h"
#include "libmesh/namebased_io.h"

namespace GRINS
{
  RayfireMesh::RayfireMesh(libMesh::Point& origin, libMesh::Real theta, libMesh::Real phi) :
    _dim(3),
    _origin(origin),
    _theta(theta),
    _phi(phi),
    _output_filename("")
  {
    libmesh_not_implemented();
  }


  RayfireMesh::RayfireMesh(libMesh::Point& origin, libMesh::Real theta) :
    _dim(2),
    _origin(origin),
    _theta(theta),
    _phi(-7.0), // bounds on angles are +/- 2pi
    _output_filename("")
  {
    if (std::abs(_theta) > 2.0*Constants::pi)
      libmesh_error_msg("Please supply a theta value between -2*pi and 2*pi");
  }

  RayfireMesh::RayfireMesh(const GetPot & input, const std::string & qoi_string) :
    _dim(2),
    _phi(-7.0),
    _output_filename(input("QoI/"+qoi_string+"/Rayfire/output_filename",""))
  {
    unsigned int rayfire_dim = input.vector_variable_size("QoI/"+qoi_string+"/Rayfire/origin");

    if (rayfire_dim != 2)
      libmesh_error_msg("ERROR: Only 2D Rayfires are currently supported");

    if (!input.have_variable("QoI/"+qoi_string+"/Rayfire/origin"))
      libmesh_error_msg("ERROR: No origin specified for Rayfire");

    if (input.vector_variable_size("QoI/"+qoi_string+"/Rayfire/origin") != 2)
      libmesh_error_msg("ERROR: Please specify a 2D point (x,y) for the rayfire origin");

    libMesh::Point origin;
    _origin(0) = input("QoI/"+qoi_string+"/Rayfire/origin", 0.0, 0);
    _origin(1) = input("QoI/"+qoi_string+"/Rayfire/origin", 0.0, 1);

    if (input.have_variable("QoI/"+qoi_string+"/Rayfire/theta"))
      _theta = input("QoI/"+qoi_string+"/Rayfire/theta", -7.0);
    else
      libmesh_error_msg("ERROR: Spherical polar angle theta must be given for Rayfire");

    if (std::abs(_theta) > 2.0*Constants::pi)
      libmesh_error_msg("Please supply a theta value between -2*pi and 2*pi");

    if (input.have_variable("QoI/"+qoi_string+"/Rayfire/phi"))
      libmesh_error_msg("ERROR: cannot specify spherical azimuthal angle phi for Rayfire, only 2D is currently supported");
  }

  void RayfireMesh::init(const libMesh::MeshBase& mesh_base)
  {
    // check if rayfire has already been initialized
    if (_elem_id_map.size() == 0)
    {
      // consistency check
      if(mesh_base.mesh_dimension() != _dim)
        {
          std::stringstream ss;
          ss <<"The supplied mesh object is " <<mesh_base.mesh_dimension()
             <<"D, but the RayfireMesh object was created with the "
             <<_dim <<"D constructor";

          libmesh_error_msg(ss.str());
        }

      _mesh.reset( new libMesh::Mesh(mesh_base.comm(),(unsigned char)1) );

      libMesh::Point start_point(_origin);

      // get first element
      const libMesh::Elem* start_elem = this->get_start_elem(mesh_base);

      if (!start_elem)
        libmesh_error_msg("Origin is not on mesh");

      // ensure the origin is on a boundary element
      // AND on the boundary of said element
      this->check_origin_on_boundary(start_elem);

      // add the origin point to the point list
      libMesh::Node * start_node = _mesh->add_point(start_point);
      libMesh::Node * end_node = NULL;

      libMesh::Point end_point;

      const libMesh::Elem* next_elem;
      const libMesh::Elem* prev_elem = start_elem;

      do
        {
          // calculate the end point and
          // get the next elem in the rayfire
          next_elem = this->get_next_elem(prev_elem,start_point,end_point);

#ifndef NDEBUG
          // make sure we are only picking up active elements
          if (next_elem)
            libmesh_assert( next_elem->active() );
#endif

          // add end point as node on the rayfire mesh
          end_node = _mesh->add_point(end_point);
          libMesh::Elem* elem = _mesh->add_elem(new libMesh::Edge2);
          elem->set_node(0) = start_node;
          elem->set_node(1) = end_node;

          // warn if rayfire elem is shorter than TOLERANCE
          if ( (start_point-end_point).norm() < libMesh::TOLERANCE)
            {
              std::stringstream ss;
              ss  <<"********\n"
                  <<"WARNING\n"
                  <<"Detected rayfire element shorter than TOLERANCE\n"
                  <<"Element ID: " <<prev_elem->id() <<", rayfire element length: " <<(start_point-end_point).norm();

              libmesh_warning(ss.str());
            }

          // add new rayfire elem to the map
          _elem_id_map[prev_elem->id()] = elem->id();

          start_point = end_point;
          start_node = end_node;
          prev_elem = next_elem;
        } while(next_elem);

      _mesh->prepare_for_use();

      if ( !(_output_filename.empty()) )
      {
        libMesh::NameBasedIO io(*(_mesh.get()));
        io.write(_output_filename);
      }

    }

  }


  const libMesh::Elem* RayfireMesh::map_to_rayfire_elem(const libMesh::dof_id_type elem_id)
  {
    return this->get_rayfire_elem(elem_id);
  }


  void RayfireMesh::elem_ids_in_rayfire(std::vector<libMesh::dof_id_type>& id_vector) const
  {
    std::map<libMesh::dof_id_type,libMesh::dof_id_type>::const_iterator it = _elem_id_map.begin();
    for(; it != _elem_id_map.end(); it++)
      {
        if (_mesh->elem_ptr(it->second)->active())
          id_vector.push_back(it->first);
      }
  }


  void RayfireMesh::reinit(const libMesh::MeshBase& mesh_base)
  {
    // we don't want to reinit() multiple times
    // at the same AMR level
    bool do_reinit = true;

    // store the elems to be refined until later
    // so we don't mess with the _elem_id_map while we
    // iterate over it
    std::vector<std::pair<const libMesh::Elem*,libMesh::Elem*> > elems_to_refine;

    // same with elems to coarsen
    // store the main_mesh elem
    std::vector<const libMesh::Elem*> elems_to_coarsen;

    // iterate over all main elems along the rayfire and look for
    // refinement: INACTIVE parent with JUST_REFINED children
    // coarsening: JUST_COARSENED parent
    std::map<libMesh::dof_id_type,libMesh::dof_id_type>::iterator it = _elem_id_map.begin();
    for(; it != _elem_id_map.end(); it++)
      {
        const libMesh::Elem * main_elem = mesh_base.elem_ptr(it->first);
        libmesh_assert(main_elem);
        
        libMesh::Elem * rayfire_elem = _mesh->elem_ptr(it->second);
        libmesh_assert(rayfire_elem);

        if (main_elem->parent())
          {
            if (main_elem->parent()->refinement_flag() == libMesh::Elem::RefinementState::JUST_COARSENED)
              {
                // if the rayfire_elem is INACTIVE, then we have already done a reinit() at this AMR level
                if (rayfire_elem->refinement_flag() == libMesh::Elem::RefinementState::INACTIVE)
                  {
                    do_reinit = false;
                    break;
                  }
                else
                  elems_to_coarsen.push_back(main_elem);
              }
          }

        libMesh::Elem::RefinementState state = main_elem->refinement_flag();

        if (state == libMesh::Elem::RefinementState::INACTIVE)
          {
            if (main_elem->has_children())
              if (main_elem->child_ptr(0)->refinement_flag() == libMesh::Elem::RefinementState::JUST_REFINED)
                {
                  // if the rayfire_elem is INACTIVE, then we have already done a reinit() at this AMR level
                  if (rayfire_elem->refinement_flag() == libMesh::Elem::RefinementState::INACTIVE)
                    {
                      do_reinit = false;
                      break;
                    }
                  else
                    elems_to_refine.push_back(std::pair<const libMesh::Elem *, libMesh::Elem *>(main_elem,rayfire_elem));
                }
          }
      }

    if (do_reinit)
      {
        // refine the elements that need it
        for (unsigned int i=0; i<elems_to_refine.size(); i++)
          this->refine(elems_to_refine[i].first, elems_to_refine[i].second);

        // coarsen the elements that need it
        for (unsigned int i=0; i<elems_to_coarsen.size(); i++)
          this->coarsen(elems_to_coarsen[i]);
      }

  }


  // private functions

  void RayfireMesh::check_origin_on_boundary(const libMesh::Elem* start_elem)
  {
    libmesh_assert(start_elem);

    // first, make sure the elem is on a boundary
    if ( !(start_elem->on_boundary()) )
      libmesh_error_msg("The supplied origin point is not on a boundary element");

    // second, check all boundary sides of the elem
    // to see if one of them conatins the origin point
    bool valid = false;

    for (unsigned int s=0; s<start_elem->n_sides(); s++)
      {
        // neighbor_ptr() returns NULL on boundary elems
        if ( start_elem->neighbor_ptr(s) )
          continue;

        // we found a boundary elem, so make an edge and see if it contains the origin
        std::unique_ptr<const libMesh::Elem> edge_elem = start_elem->build_edge_ptr(s);
        valid |= edge_elem->contains_point(_origin);
      }

    if (!valid)
      libmesh_error_msg("The supplied origin point is not on the boundary of the starting element");
  }


  const libMesh::Elem* RayfireMesh::get_start_elem(const libMesh::MeshBase& mesh_base)
  {
    const libMesh::Elem* start_elem = NULL;

    std::unique_ptr<libMesh::PointLocatorBase> locator = mesh_base.sub_point_locator();
    const libMesh::Elem* elem = (*locator)(_origin);

    // elem would be NULL if origin is not on mesh
    if (elem)
      {
        if (this->rayfire_in_elem(_origin,elem))
          start_elem = elem;
        else
          {
            for (unsigned int i=0; i<elem->n_neighbors(); i++)
              {
                const libMesh::Elem* neighbor_elem = elem->neighbor_ptr(i);
                if (!neighbor_elem)
                  continue;

                if (this->rayfire_in_elem(_origin,neighbor_elem))
                  {
                    start_elem = neighbor_elem;
                    break;
                  }
              }
          }
      }

    return start_elem; // might be NULL if _origin is not on mesh
  }


  libMesh::Elem* RayfireMesh::get_rayfire_elem(const libMesh::dof_id_type elem_id)
  {
    // return value; set if valid rayfire elem is found
    libMesh::Elem* retval = NULL;

    std::map<libMesh::dof_id_type,libMesh::dof_id_type>::iterator it;
    it = _elem_id_map.find(elem_id);
    if (it != _elem_id_map.end())
      if (_mesh->elem_ptr(it->second)->refinement_flag() != libMesh::Elem::RefinementState::INACTIVE)
        retval = _mesh->elem_ptr(it->second);

    return retval;
  }


  const libMesh::Elem* RayfireMesh::get_next_elem(const libMesh::Elem* cur_elem, libMesh::Point& start_point, libMesh::Point& next_point, bool same_parent)
  {
    libmesh_assert(cur_elem);

    unsigned int intersection_side = this->calculate_intersection_point(start_point,cur_elem,next_point);

    if (intersection_side == libMesh::invalid_uint)
      {
        std::stringstream ss;
        ss  <<"ERROR: Could not find next element along rayfire" <<std::endl
            <<"Current element ID: " <<cur_elem->id() <<std::endl
            <<"Current point: " <<start_point <<std::endl;

        libmesh_error_msg(ss.str());
      }

    // will return NULL if intersection_side is a mesh boundary
    return this->get_correct_neighbor(start_point,next_point,cur_elem,intersection_side,same_parent);
  }


  bool RayfireMesh::check_valid_point(libMesh::Point& intersection_point, libMesh::Point& start_point, const libMesh::Elem& edge_elem, libMesh::Point& next_point)
  {
    bool is_not_start = !(intersection_point.absolute_fuzzy_equals(start_point));
    bool is_on_edge = edge_elem.contains_point(intersection_point);

    bool is_valid = is_not_start && is_on_edge;

    if ( is_valid )
      {
        next_point(0) = intersection_point(0);
        next_point(1) = intersection_point(1);
      }

    return is_valid;
  }


  bool RayfireMesh::rayfire_in_elem(const libMesh::Point& end_point, const libMesh::Elem* elem)
  {
    libmesh_assert(elem);

    // move a little bit along the rayfire and see if we are still in the elem
    // need to move more than TOLERANCE to avoid false positive on contains_point()
    libMesh::Real L = 2*libMesh::TOLERANCE;

    // If the elem is too small, need to shorten L so we stay within the elem
    if ( elem->hmin() < libMesh::TOLERANCE )
      L = elem->hmin() * 0.1;

    // parametric representation of rayfire line
    libMesh::Real x = end_point(0) + L*std::cos(_theta);
    libMesh::Real y = end_point(1) + L*std::sin(_theta);

    return elem->contains_point(libMesh::Point(x,y));
  }


  bool RayfireMesh::validate_edge(const libMesh::Point & start_point, const libMesh::Point & end_point, const libMesh::Elem * side_elem, const libMesh::Elem * neighbor)
  {
    bool is_valid = true;

    const libMesh::Node * node0 = side_elem->node_ptr(0);
    const libMesh::Node * node1 = side_elem->node_ptr(1);

    unsigned int side = libMesh::invalid_uint;

    for (unsigned int s=0; s<neighbor->n_sides(); s++)
      {
        std::unique_ptr<const libMesh::Elem> side_elem = neighbor->build_side_ptr(s);

        if ( (side_elem->contains_point(*node0)) && (side_elem->contains_point(*node1)) )
          {
            side = s;
            break;
          }
      }

    if ( side != libMesh::invalid_uint )
      {
        std::unique_ptr<const libMesh::Elem> edge_elem = neighbor->build_side_ptr(side);

        bool start_point_on_edge = edge_elem->contains_point(start_point);
        bool end_point_on_edge = edge_elem->contains_point(end_point);

        if ( end_point_on_edge && start_point_on_edge )
          {
            bool l_to_r = ( start_point.absolute_fuzzy_equals(*(edge_elem->node_ptr(0))) ) && ( end_point.absolute_fuzzy_equals(*(edge_elem->node_ptr(1))) );
            bool r_to_l = ( end_point.absolute_fuzzy_equals(*(edge_elem->node_ptr(0))) )   && ( start_point.absolute_fuzzy_equals(*(edge_elem->node_ptr(1))) );

            is_valid &= ( l_to_r || r_to_l );
          }
      }

    return is_valid;
  }


  const libMesh::Elem* RayfireMesh::get_correct_neighbor(libMesh::Point & start_point, libMesh::Point & end_point, const libMesh::Elem * cur_elem, unsigned int side, bool same_parent)
  {
    libmesh_assert(cur_elem);
    libmesh_assert(cur_elem->active());

    const libMesh::Elem * neighbor = cur_elem->neighbor_ptr(side);

    // a valid neighbor is an ACTIVE elem or null
    bool is_valid = true;
    if (neighbor)
      is_valid = neighbor->active();

    if (is_valid)
      {
        // check if the intersection point is a vertex
        bool is_vertex = false;
        const libMesh::Node * vertex = NULL;
        for(unsigned int n=0; n<cur_elem->n_nodes(); n++)
          {
            if ((cur_elem->node_ptr(n))->absolute_fuzzy_equals(end_point))
              {
                is_vertex = true;
                vertex = cur_elem->node_ptr(n);
                break;
              }
          }

        if (is_vertex)
          {
            // rayfire goes through vertex

            // check elem neighbors first
            for (unsigned int s=0; s<cur_elem->n_sides(); s++)
              {
                std::unique_ptr<const libMesh::Elem> side_elem = cur_elem->build_side_ptr(s);
                if (side_elem->contains_point(end_point))
                  {
                    const libMesh::Elem * neighbor = cur_elem->neighbor_ptr(s);
                    if (!neighbor)
                      continue;

                    if (same_parent)
                      {
                        if (neighbor->parent())
                          {
                            if (neighbor->parent()->id() != cur_elem->parent()->id())
                              continue;
                            else
                              if (this->rayfire_in_elem(end_point,neighbor))
                                if (this->validate_edge(start_point,end_point,side_elem.get(),neighbor))
                                  return neighbor;
                          }
                        else
                          continue;
                      }
                    else
                      if (this->rayfire_in_elem(end_point,neighbor))
                        if (this->validate_edge(start_point,end_point,side_elem.get(),neighbor))
                          return neighbor;
                  }
              }

            // next elem is not a neighbor,
            // so get all elems that share this vertex
            std::set<const libMesh::Elem*> elem_set;
            cur_elem->find_point_neighbors(*vertex,elem_set);
            std::set<const libMesh::Elem *>::const_iterator       it  = elem_set.begin();
            const std::set<const libMesh::Elem *>::const_iterator end = elem_set.end();

            // iterate over each elem
            for (; it != end; ++it)
              {
                const libMesh::Elem* elem = *it;

                if (elem == cur_elem) // skip the current elem
                  continue;

                if (this->rayfire_in_elem(end_point,elem))
                  return elem;
              }

          }
        else
          {
            // check if side is a boundary
            if( !(cur_elem->neighbor_ptr(side)) )
              return NULL;

            // not a vertex or boundary, so just get the elem on that side
            return cur_elem->neighbor_ptr(side);
          }
      } // if (is_valid)
    else
      {
        for (unsigned int c=0; c<neighbor->n_children(); c++)
          {
            const libMesh::Elem * child = neighbor->child_ptr(c);

            if (child->contains_point(end_point))
              if (this->rayfire_in_elem(end_point,child))
                return child;
          }
        // could not find a valid child elem, so we return NULL below
      }

    return NULL;
  }


  unsigned int RayfireMesh::calculate_intersection_point(libMesh::Point & initial_point, const libMesh::Elem * cur_elem, libMesh::Point & intersection_point)
  {
    libmesh_assert(cur_elem);

    unsigned int intersect_side = libMesh::invalid_uint;

    switch(cur_elem->dim())
      {
        case 2:
          switch(cur_elem->default_order())
            {
              case libMesh::Order::FIRST:
                intersect_side = this->intersection_2D_first_order(initial_point,cur_elem,intersection_point);
                break;
              case libMesh::Order::SECOND:
                intersect_side = this->intersection_2D_second_order(initial_point,cur_elem,intersection_point);
                break;
              default:
                libmesh_error_msg("Unknown/unsupported 2D element order");
                break;
            }
          break;
        default:
          libmesh_error_msg("Unknown/unsupported element dim");
      }

    return intersect_side;
  }


  unsigned int RayfireMesh::intersection_2D_first_order(libMesh::Point & initial_point, const libMesh::Elem * cur_elem, libMesh::Point & intersection_point)
  {
    libmesh_assert(cur_elem);
    libmesh_assert_equal_to(cur_elem->dim(),2);
    libmesh_assert(cur_elem->default_order() == libMesh::Order::FIRST);

    // return value
    unsigned int intersect_side = libMesh::invalid_uint;

    // precompute repeated term
    libMesh::Real tan_theta = std::tan(_theta);

    // cos(theta)==0 when tan(theta)==Inf
    // Truncating pi can result in tan(theta) being very large but not Inf,
    // so we can check cos(theta)<TOLERANCE to identify vertical rayfires
    libMesh::Real cos_theta = std::cos(_theta);

    // rayfire starting points
    libMesh::Real xr = initial_point(0);
    libMesh::Real yr = initial_point(1);

    for (unsigned int s=0; s<cur_elem->n_sides(); ++s)
      {
        std::unique_ptr<const libMesh::Elem> edge_elem = cur_elem->build_edge_ptr(s);

        // using the default tol can cause a false positive when start_point is near a node,
        // causing this loop to skip over an otherwise valid edge to check
        if (edge_elem->contains_point(initial_point,libMesh::TOLERANCE*0.1))
          continue;

        // since cur_elem is first order, the sides are always linear
        // and can be represented in point-slope form
        libMesh::Real x0 = (edge_elem->node_ref(0))(0);
        libMesh::Real y0 = (edge_elem->node_ref(0))(1);
        libMesh::Real x1 = (edge_elem->node_ref(1))(0);
        libMesh::Real y1 = (edge_elem->node_ref(1))(1);

        // slope of edge
        libMesh::Real m = (y1-y0)/(x1-x0);

        // intersection point
        libMesh::Real x_hat;
        libMesh::Real y_hat;

        if ( std::abs(cos_theta) < libMesh::TOLERANCE )
          {
            // rayfire is vertical, so we know x_hat == xr
            x_hat = xr;
            y_hat = y0 + m*(x_hat-x0);
          }
        else
          {
            if ( std::isinf(m) )
              {
                // edge is vertical, so we already know the x coordinate of the intersection
                x_hat = x0;
              }
            else
              {
                if ( std::abs(tan_theta - m) < libMesh::TOLERANCE )
                  continue; // rayfire is parallel to this edge

                x_hat = ( xr*tan_theta - yr + y0 - m*x0 )/( tan_theta - m );
              }

            y_hat = yr + (x_hat-xr)*tan_theta;

          }

        libMesh::Point intersect(x_hat,y_hat);

        if (edge_elem->contains_point(intersect,libMesh::TOLERANCE*0.1))
          {
            intersect_side = s;

            intersection_point(0) = intersect(0);
            intersection_point(1) = intersect(1);
            break;
          }

      }

    return intersect_side;
  }


  unsigned int RayfireMesh::intersection_2D_second_order(libMesh::Point & initial_point, const libMesh::Elem * cur_elem, libMesh::Point & intersection_point)
  {
    libmesh_assert(cur_elem);
    libmesh_assert_equal_to(cur_elem->dim(),2);
    libmesh_assert(cur_elem->default_order() == libMesh::Order::SECOND);

    // return value
    unsigned int intersect_side = libMesh::invalid_uint;

    // based on limit from libMesh::FEInterface::inverse_map()
    unsigned int iter_max = 20;

    // precompute repeated terms
    libMesh::Real sin_theta = std::sin(_theta);
    libMesh::Real cos_theta = std::cos(_theta);

    // rayfire starting points
    libMesh::Real xr = initial_point(0);
    libMesh::Real yr = initial_point(1);

    for (unsigned int s=0; s<cur_elem->n_sides(); ++s)
      {
        std::unique_ptr<const libMesh::Elem> edge_elem = cur_elem->build_edge_ptr(s);

        // using the default tol can cause a false positive when start_point is near a node,
        // causing this loop to skip over an otherwise valid edge to check
        if (edge_elem->contains_point(initial_point,libMesh::TOLERANCE*0.1))
          continue;

        // the number of shape functions needed for the edge_elem
        unsigned int n_sf = libMesh::FE<1,libMesh::LAGRANGE>::n_shape_functions(edge_elem->type(),edge_elem->default_order());

        // shape functions and derivatives w.r.t reference coordinate
        std::vector<libMesh::Real> phi(n_sf);
        std::vector<libMesh::Real> dphi(n_sf);

        // Initial xi guess is center of the edge_elem
        libMesh::Real xi = 0.0;

        // Initial L guess
        libMesh::Real L = std::abs(xr-(edge_elem->node_ref(2))(0));

        // Newton iteration
        for(unsigned int it=0; it<iter_max; it++)
          {
            // Get the shape function and derivative values at the reference coordinate
            // phi.size() == dphi.size()
            for(unsigned int i=0; i<phi.size(); i++)
              {
                phi[i] = libMesh::FE<1,libMesh::LAGRANGE>::shape(edge_elem->type(),
                                                                 edge_elem->default_order(),
                                                                 i,
                                                                 xi);

                dphi[i] = libMesh::FE<1,libMesh::LAGRANGE>::shape_deriv(edge_elem->type(),
                                                                        edge_elem->default_order(),
                                                                        i,
                                                                        0, // const unsigned int libmesh_dbg_varj
                                                                        xi);
              } // for i

            libMesh::Real X = 0.0;
            libMesh::Real Y = 0.0;
            libMesh::Real dX_dxi = 0.0;
            libMesh::Real dY_dxi = 0.0;

            for(unsigned int i=0; i<phi.size(); i++)
              {
                X  += (*(edge_elem->node_ptr(i)))(0) * phi[i];
                Y  += (*(edge_elem->node_ptr(i)))(1) * phi[i];
                dX_dxi += (*(edge_elem->node_ptr(i)))(0) * dphi[i];
                dY_dxi += (*(edge_elem->node_ptr(i)))(1) * dphi[i];
              }

          libMesh::DenseVector<libMesh::Real> F(2);
          F(0) = xr + L*cos_theta - X;
          F(1) = yr + L*sin_theta - Y;

          // Jacobian entries
          libMesh::Real a = cos_theta;
          libMesh::Real b = -dX_dxi;
          libMesh::Real c = sin_theta;
          libMesh::Real d = -dY_dxi;

          libMesh::Real J_det = a*d - b*c;

          // inverse of the Jacobian
          libMesh::DenseMatrix<libMesh::Real> J_inv(2,2);
          J_inv(0,0) =  d;
          J_inv(0,1) = -b;
          J_inv(1,0) = -c;
          J_inv(1,1) =  a;
          J_inv *= 1.0/J_det;

          // delta will be the newton step
          libMesh::DenseVector<libMesh::Real> delta(2);
          J_inv.vector_mult(delta,F);

          // check for convergence
          if ( delta.l2_norm() < libMesh::TOLERANCE )
            {
              libMesh::Point intersect(X,Y);

              // newton solver converged, now make sure it converged to a point on the edge_elem
              if (edge_elem->contains_point(intersect,libMesh::TOLERANCE*0.1))
              {
                intersect_side = s;

                intersection_point(0) = intersect(0);
                intersection_point(1) = intersect(1);
              }
              break; // break out of 'for it'
            }
          else
            {
              L  -= delta(0);
              xi -= delta(1);
            }
        } //for it

        if (intersect_side != libMesh::invalid_uint)
          break; // break out of 'for s'

      } // for s

    return intersect_side;
  }


  void RayfireMesh::refine(const libMesh::Elem* main_elem, libMesh::Elem* rayfire_elem)
  {
    libmesh_assert(main_elem);
    libmesh_assert(rayfire_elem);

    libmesh_assert_equal_to(main_elem->refinement_flag(),libMesh::Elem::RefinementState::INACTIVE);

    // these nodes cannot change
    libMesh::Node* start_node = rayfire_elem->node_ptr(0);
    libMesh::Node* end_node   = rayfire_elem->node_ptr(1);

    // set the rayfire_elem as INACTIVE
    rayfire_elem->set_refinement_flag(libMesh::Elem::RefinementState::INACTIVE);

    // find which child elem we start with
    libMesh::dof_id_type start_child = libMesh::invalid_uint;
    for(unsigned int i=0; i<main_elem->n_children(); i++)
      {
        if ( (main_elem->child_ptr(i))->contains_point(*start_node) )
          {
            // check if the start point is a vertex
            bool is_vertex = false;
            for(unsigned int n=0; n<main_elem->child_ptr(i)->n_nodes(); n++)
              {
                if ((main_elem->child_ptr(i)->node_ptr(n))->absolute_fuzzy_equals(*start_node))
                  {
                    is_vertex = true;
                    break;
                  }
              }

            if (is_vertex)
              {
                if ( this->rayfire_in_elem(*start_node, main_elem->child_ptr(i)) )
                  {
                    start_child = i;
                    break;
                  }
              }
            else
              {
                start_child = i;
                break;
              }
          }
      }

    // make sure we found a starting child elem
    libmesh_assert_not_equal_to(start_child,libMesh::invalid_uint);

    // we found the starting element, so perform the rayfire
    // until we reach the stored end_node
    libMesh::Point start_point = *start_node;
    libMesh::Point end_point;

    const libMesh::Elem* next_elem;
    const libMesh::Elem* prev_elem = main_elem->child_ptr(start_child);

    // if prev_elem is INACTIVE, then more than one refinement
    // has taken place between reinit() calls and will
    // break this
    libmesh_assert_equal_to( prev_elem->refinement_flag(), libMesh::Elem::RefinementState::JUST_REFINED );

    libMesh::Node * prev_node = start_node;
    libMesh::Node * new_node = NULL;

    // iterate until we reach the stored end_node
    do
      {
        next_elem = this->get_next_elem(prev_elem,start_point,end_point,true);

#ifndef NDEBUG
        // make sure we are only picking up active elements
        if (next_elem)
          libmesh_assert( next_elem->active() );
#endif

        // add end point as node on the rayfire mesh
        new_node = _mesh->add_point(end_point);
        libMesh::Elem* elem = _mesh->add_elem(new libMesh::Edge2);

        elem->set_node(0) = prev_node;
        elem->set_node(1) = new_node;

        libmesh_assert_less( (*(elem->node_ptr(0))-_origin).norm(),  (*(elem->node_ptr(1))-_origin).norm());

        // warn if rayfire elem is shorter than TOLERANCE
        if ( (start_point-end_point).norm() < libMesh::TOLERANCE)
          {
            std::stringstream ss;
            ss  <<"********\n"
                <<"WARNING\n"
                <<"Detected rayfire element shorter than TOLERANCE on refinement\n"
                <<"Element ID: " <<prev_elem->id() <<", rayfire element length: " <<(start_point-end_point).norm();

            libmesh_warning(ss.str());
          }

        // set rayfire_elem as the parent of this new elem
        // in case it gets coarsened
        elem->set_parent(rayfire_elem);

        // add new rayfire elem to the map
        _elem_id_map[prev_elem->id()] = elem->id();
        start_point = end_point;
        prev_elem = next_elem;
        prev_node = new_node;

      } while(!(new_node->absolute_fuzzy_equals(*end_node)));

  }


  void RayfireMesh::coarsen(const libMesh::Elem* child_elem)
  {
    libmesh_assert(child_elem);

    if (this->get_rayfire_elem(child_elem->id()))
      {
        const libMesh::Elem* parent_elem = child_elem->parent();
        libmesh_assert(parent_elem);

        const libMesh::Node* start_node = NULL;
        const libMesh::Node* end_node = NULL;

        for (unsigned int c=0; c<parent_elem->n_children(); c++)
          {
            libMesh::Elem* rayfire_child = this->get_rayfire_elem(parent_elem->child_ptr(c)->id());

            if (rayfire_child)
              {
                if (!start_node)
                  {
                    start_node = rayfire_child->node_ptr(0);
                    end_node = rayfire_child->node_ptr(1);
                  }
                else
                  {
                    if ( (this->_origin - *(rayfire_child->node_ptr(0))).norm() < (this->_origin - *start_node).norm() )
                      start_node = rayfire_child->node_ptr(0);

                    if ( (this->_origin - *(rayfire_child->node_ptr(1))).norm() > (this->_origin - *end_node).norm() )
                      end_node = rayfire_child->node_ptr(1);
                  }

                rayfire_child->set_refinement_flag(libMesh::Elem::RefinementState::INACTIVE);
              }
          } // for c

        // make sure we found nodes
        libmesh_assert(start_node);
        libmesh_assert(end_node);

        // add a new rayfire elem
        libMesh::Elem* elem = _mesh->add_elem(new libMesh::Edge2);
        elem->set_node(0) = _mesh->node_ptr(start_node->id());
        elem->set_node(1) = _mesh->node_ptr(end_node->id());

        libmesh_assert_less( (*(elem->node_ptr(0))-_origin).norm(),  (*(elem->node_ptr(1))-_origin).norm());

        // add new rayfire elem to the map
        _elem_id_map[parent_elem->id()] = elem->id();
      }

  }

} //namespace GRINS
