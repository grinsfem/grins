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
#include "grins/rayfire_mesh.h"

// GRINS
#include "grins/math_constants.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/elem.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/fe.h"

namespace GRINS
{
  RayfireMesh::RayfireMesh(libMesh::Point& origin, libMesh::Real theta, libMesh::Real phi) :
    _dim(3),
    _origin(origin),
    _theta(theta),
    _phi(phi)
  {
    libmesh_not_implemented();
  }


  RayfireMesh::RayfireMesh(libMesh::Point& origin, libMesh::Real theta) :
    _dim(2),
    _origin(origin),
    _theta(theta),
    _phi(-7.0) // bounds on angles are +/- 2pi
  {
    if (std::abs(_theta) > 2.0*Constants::pi)
      libmesh_error_msg("Please supply a theta value between -2*pi and 2*pi");
  }

  void RayfireMesh::init(const libMesh::MeshBase& mesh_base)
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

    _mesh = new libMesh::Mesh(mesh_base.comm(),(unsigned char)1);

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

  }


  const libMesh::Elem* RayfireMesh::map_to_rayfire_elem(const libMesh::dof_id_type elem_id)
  {
    return this->get_rayfire_elem(elem_id);
  }


  void RayfireMesh::elem_ids_in_rayfire(std::vector<libMesh::dof_id_type>& id_vector)
  {
    std::map<libMesh::dof_id_type,libMesh::dof_id_type>::iterator it = _elem_id_map.begin();
    for(; it != _elem_id_map.end(); it++)
      {
        if (_mesh->elem(it->second)->active())
          id_vector.push_back(it->first);
      }
  }


  void RayfireMesh::reinit(const libMesh::MeshBase& mesh_base)
  {
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
        const libMesh::Elem* main_elem = mesh_base.elem(it->first);
        libmesh_assert(main_elem);

        if (main_elem->parent())
          {
            if (main_elem->parent()->refinement_flag() == libMesh::Elem::RefinementState::JUST_COARSENED)
              elems_to_coarsen.push_back(main_elem);
          }

        libMesh::Elem::RefinementState state = main_elem->refinement_flag();

        if (state == libMesh::Elem::RefinementState::INACTIVE)
          {
            if (main_elem->has_children())
              if (main_elem->child(0)->refinement_flag() == libMesh::Elem::RefinementState::JUST_REFINED)
                elems_to_refine.push_back(std::pair<const libMesh::Elem*, libMesh::Elem*>(main_elem,_mesh->elem(it->second)));
          }
      }

    // refine the elements that need it
    for (unsigned int i=0; i<elems_to_refine.size(); i++)
      this->refine(elems_to_refine[i].first, elems_to_refine[i].second);

    // coarsen the elements that need it
    for (unsigned int i=0; i<elems_to_coarsen.size(); i++)
      this->coarsen(elems_to_coarsen[i]);

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
        libMesh::UniquePtr<libMesh::Elem> edge_elem = start_elem->build_edge(s);
        valid |= edge_elem->contains_point(_origin);
      }

    if (!valid)
      libmesh_error_msg("The supplied origin point is not on the boundary of the starting element");
  }


  const libMesh::Elem* RayfireMesh::get_start_elem(const libMesh::MeshBase& mesh_base)
  {
    const libMesh::Elem* start_elem = NULL;

    libMesh::UniquePtr<libMesh::PointLocatorBase> locator = mesh_base.sub_point_locator();
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
      if (_mesh->elem(it->second)->refinement_flag() != libMesh::Elem::RefinementState::INACTIVE)
        retval = _mesh->elem(it->second);

    return retval;
  }


  const libMesh::Elem* RayfireMesh::get_next_elem(const libMesh::Elem* cur_elem, libMesh::Point& start_point, libMesh::Point& next_point, bool same_parent)
  {
    libmesh_assert(cur_elem);

    libMesh::Point intersection_point;

    // loop over all sides of the elem and check each one for intersection
    for (unsigned int s=0; s<cur_elem->n_sides(); s++)
      {
        const libMesh::UniquePtr<libMesh::Elem> edge_elem = cur_elem->build_edge(s);

        // Using the default tol can cause a false positive when start_point is near a node,
        // causing this loop to skip over an otherwise valid edge to check
        if (edge_elem->contains_point(start_point,libMesh::TOLERANCE*0.1))
          continue;

        bool converged = this->newton_solve_intersection(start_point,edge_elem.get(),intersection_point);

        if (converged)
          {
            if ( this->check_valid_point(intersection_point,start_point,*edge_elem,next_point) )
              return this->get_correct_neighbor(intersection_point,cur_elem,s,same_parent);
          }
        else
          continue;

      } // for s

    return NULL; // no intersection
  }


  bool RayfireMesh::check_valid_point(libMesh::Point& intersection_point, libMesh::Point& start_point, libMesh::Elem& edge_elem, libMesh::Point& next_point)
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


  const libMesh::Elem* RayfireMesh::get_correct_neighbor(libMesh::Point& end_point, const libMesh::Elem* cur_elem, unsigned int side, bool same_parent)
  {
    libmesh_assert(cur_elem);

    // check if the intersection point is a vertex
    bool is_vertex = false;
    libMesh::Node * vertex = NULL;
    for(unsigned int n=0; n<cur_elem->n_nodes(); n++)
      {
        if ((cur_elem->get_node(n))->absolute_fuzzy_equals(end_point))
          {
            is_vertex = true;
            vertex = cur_elem->get_node(n);
            break;
          }
      }

    if (is_vertex)
      {
        // rayfire goes through vertex

        // check elem neighbors first
        for (unsigned int s=0; s<cur_elem->n_sides(); s++)
          {
            libMesh::UniquePtr<libMesh::Elem> side_elem = cur_elem->build_side(s);
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
                            return neighbor;
                      }
                    else
                      continue;
                  }
                else
                  if (this->rayfire_in_elem(end_point,neighbor))
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

    return NULL;
  }


  bool RayfireMesh::newton_solve_intersection(libMesh::Point& initial_point, const libMesh::Elem* edge_elem, libMesh::Point& intersection_point)
  {
    libmesh_assert(edge_elem);

    unsigned int iter_max = 20; // max iterations

    // the number of shape functions needed for the edge_elem
    unsigned int n_sf = libMesh::FE<1,libMesh::LAGRANGE>::n_shape_functions(edge_elem->type(),edge_elem->default_order());

    // starting point on the elem
    libMesh::Real x0 = initial_point(0);
    libMesh::Real y0 = initial_point(1);

    // shape functions and derivatives w.r.t reference coordinate
    std::vector<libMesh::Real> phi(n_sf);
    std::vector<libMesh::Real> dphi(n_sf);

    // Newton iteration step
    libMesh::Real d_xi;

    // tan(theta) is the slope, so precompute since it is used repeatedly
    libMesh::Real tan_theta = std::tan(_theta);

    // Initial guess is center of the edge_elem
    libMesh::Real xi = 0.0;

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

        libMesh::Real X=0.0, Y=0.0, dX=0.0, dY=0.0;

        for(unsigned int i=0; i<phi.size(); i++)
          {
            X  += (*(edge_elem->get_node(i)))(0) * phi[i];
            dX += (*(edge_elem->get_node(i)))(0) * dphi[i];
            Y  += (*(edge_elem->get_node(i)))(1) * phi[i];
            dY += (*(edge_elem->get_node(i)))(1) * dphi[i];
          }

        libMesh::Real  f = tan_theta*(X-x0) - (Y-y0);
        libMesh::Real df = tan_theta*(dX) - dY;

        d_xi = f/df;

        if(std::abs(d_xi) < libMesh::TOLERANCE)
          {
            // convergence
            intersection_point(0) = X;
            intersection_point(1) = Y;
            return true;
          }
        else
          xi -= d_xi;

      } // for it

    // no convergence
    return false;
  }


  void RayfireMesh::refine(const libMesh::Elem* main_elem, libMesh::Elem* rayfire_elem)
  {
    libmesh_assert(main_elem);
    libmesh_assert(rayfire_elem);

    libmesh_assert_equal_to(main_elem->refinement_flag(),libMesh::Elem::RefinementState::INACTIVE);

    // these nodes cannot change
    libMesh::Node* start_node = rayfire_elem->get_node(0);
    libMesh::Node* end_node   = rayfire_elem->get_node(1);

    // set the rayfire_elem as INACTIVE
    rayfire_elem->set_refinement_flag(libMesh::Elem::RefinementState::INACTIVE);

    // find which child elem we start with
    libMesh::dof_id_type start_child = libMesh::invalid_uint;
    for(unsigned int i=0; i<main_elem->n_children(); i++)
      {
        if ( (main_elem->child(i))->contains_point(*start_node) )
          {
            // check if the start point is a vertex
            bool is_vertex = false;
            for(unsigned int n=0; n<main_elem->child(i)->n_nodes(); n++)
              {
                if ((main_elem->child(i)->get_node(n))->absolute_fuzzy_equals(*start_node))
                  {
                    is_vertex = true;
                    break;
                  }
              }

            if (is_vertex)
              {
                if ( this->rayfire_in_elem(*start_node, main_elem->child(i)) )
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
    const libMesh::Elem* prev_elem = main_elem->child(start_child);

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

        // add end point as node on the rayfire mesh
        new_node = _mesh->add_point(end_point);
        libMesh::Elem* elem = _mesh->add_elem(new libMesh::Edge2);

        elem->set_node(0) = prev_node;
        elem->set_node(1) = new_node;

        libmesh_assert_less( (*(elem->get_node(0))-_origin).norm(),  (*(elem->get_node(1))-_origin).norm());

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
            libMesh::Elem* rayfire_child = this->get_rayfire_elem(parent_elem->child(c)->id());

            if (rayfire_child)
              {
                if (!start_node)
                  {
                    start_node = rayfire_child->get_node(0);
                    end_node = rayfire_child->get_node(1);
                  }
                else
                  {
                    if ( (this->_origin - *(rayfire_child->get_node(0))).norm() < (this->_origin - *start_node).norm() )
                      start_node = rayfire_child->get_node(0);

                    if ( (this->_origin - *(rayfire_child->get_node(1))).norm() > (this->_origin - *end_node).norm() )
                      end_node = rayfire_child->get_node(1);
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

        libmesh_assert_less( (*(elem->get_node(0))-_origin).norm(),  (*(elem->get_node(1))-_origin).norm());

        // add new rayfire elem to the map
        _elem_id_map[parent_elem->id()] = elem->id();
      }

  }

} //namespace GRINS
