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
#include "grins/multiphysics_sys.h"
#include "grins/assembly_context.h"
#include "grins/materials_parsing.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"
#include "libmesh/point_locator_list.h"
#include "libmesh/elem.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/analytic_function.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"

namespace GRINS
{
  RayfireMesh::RayfireMesh(libMesh::Point& origin, libMesh::Real theta, libMesh::Real phi) : 
    _dim(3),_origin(origin),_theta(theta),_phi(phi)
  {
    libmesh_not_implemented();
  }
 
 
  RayfireMesh::RayfireMesh(libMesh::Point& origin, libMesh::Real theta) : 
    _dim(2),_origin(origin),_theta(theta),_phi(0)
  {}


  void RayfireMesh::init(const libMesh::MeshBase& mesh_base)
  {
    _mesh = new libMesh::Mesh(mesh_base.comm(),(unsigned char)1);
    
    unsigned int node_id = 0;

    libMesh::Point* start_point = new libMesh::Point(_origin);

    // get first element 
    libMesh::UniquePtr<libMesh::PointLocatorBase> locator = mesh_base.sub_point_locator();
    const libMesh::Elem* start_elem = (*locator)(_origin);
        
    if (!start_elem)
      libmesh_error_msg("Origin is not on mesh");
        
    // ensure the origin is on a boundary element
    // AND on the boundary of said element
    _check_origin_on_boundary(start_elem);
        
    // add the origin point to the point list
    _mesh->add_point(*start_point,node_id++);

    libMesh::Point* end_point = new libMesh::Point();

    const libMesh::Elem* next_elem;
    const libMesh::Elem* prev_elem = start_elem;
    
    // calculate the end point and
    // get the second elem in the rayfire
    next_elem = _get_next_elem(start_elem,start_point,end_point);
    
    // add end point as node on the rayfire mesh    
    _mesh->add_point(*end_point,node_id);
    libMesh::Elem* elem = _mesh->add_elem(new libMesh::Edge2);
    elem->set_node(0) = _mesh->node_ptr(node_id-1);
    elem->set_node(1) = _mesh->node_ptr(node_id++);
    
    // add new rayfire elem to the map  
    _elem_id_map[prev_elem->id()] = elem;            
    start_point = end_point;
    prev_elem = next_elem;

    // continue looping until we reach a boundary
    while(next_elem)
    {
      next_elem = _get_next_elem(prev_elem,start_point,end_point);
            
      _mesh->add_point(*end_point,node_id);
      libMesh::Elem* elem = _mesh->add_elem(new libMesh::Edge2);
      elem->set_node(0) = _mesh->node_ptr(node_id-1);
      elem->set_node(1) = _mesh->node_ptr(node_id++);
            
      _elem_id_map[prev_elem->id()] = elem;            
      start_point = end_point;
      prev_elem = next_elem;     
    } 
  }
  
  
  const libMesh::Elem* RayfireMesh::translate(const libMesh::dof_id_type elem_id)
  {
    std::map<libMesh::dof_id_type,libMesh::Elem*>::iterator it;
    it = _elem_id_map.find(elem_id);
    if (it != _elem_id_map.end())
        return it->second;
            
    return NULL;
  }
 

  // private functions
 
  void RayfireMesh::_check_origin_on_boundary(const libMesh::Elem* start_elem)
  {
    // first, make sure the elem is on a boundary
    if ( !(start_elem->on_boundary()) )
        libmesh_error_msg("The supplied origin point is not on a boundary element");
        
    // second, check all boundary sides of the elem
    // to see if one of them conatins the origin point
    bool valid = false;
        
    for (unsigned int s=0; s<start_elem->n_sides(); s++)
    {
      // neighbor() returns NULL on boundary elems
      if ( start_elem->neighbor(s) )
        continue;
            
      // we found a boundary elem, so make an edge and see if it contains the origin                
      libMesh::UniquePtr<libMesh::Elem> edge_elem = start_elem->build_edge(s);
      valid |= edge_elem->contains_point(_origin);
    }
        
    if (!valid)
      libmesh_error_msg("The supplied origin point is not on the boundary of the starting element");
  }
  
  
  const libMesh::Elem* RayfireMesh::_get_next_elem(const libMesh::Elem* cur_elem, libMesh::Point* start_point, libMesh::Point* next_point)
  {
    libMesh::Point* intersection_point = new libMesh::Point();
    
    // loop over all sides of the elem and check each one for intersection
    for (unsigned int s=0; s<cur_elem->n_sides(); s++)
    {
      const libMesh::UniquePtr<libMesh::Elem> edge_elem = cur_elem->build_edge(s);
      if (edge_elem->contains_point(*start_point))
        continue;

      bool converged = _newton_solver(*start_point,edge_elem.get(),intersection_point);
            
      if (converged)
      {
        if ( _check_valid_point(*intersection_point,*start_point,*edge_elem,next_point) )
          return _get_correct_neighbor(*intersection_point,cur_elem,s);
      }
      else
        continue;

    } // for s

    return NULL; // no intersection          
  }

    
  bool RayfireMesh::_check_valid_point(libMesh::Point& intersection_point, libMesh::Point& start_point, libMesh::Elem& edge_elem, libMesh::Point* next_point)
  {
    bool is_not_start = !(intersection_point.absolute_fuzzy_equals(start_point));
    bool is_on_edge = edge_elem.contains_point(intersection_point);
        
    if ( is_not_start && is_on_edge )
    {
      (*next_point)(0) = intersection_point(0);
      (*next_point)(1) = intersection_point(1);
      return true;           
    }
        
    return false;
  }


  const libMesh::Elem* RayfireMesh::_get_correct_neighbor(libMesh::Point& end_point, const libMesh::Elem* cur_elem, unsigned int side)
  {
    // check if side is a boundary
    if( !(cur_elem->neighbor(side)) )
      return NULL;
  
    // check if the intersection point is a vertex
    bool is_vertex = false;
    for(unsigned int n=0; n<4; n++)
        is_vertex |= (cur_elem->get_node(n))->absolute_fuzzy_equals(end_point);
        
    if (is_vertex)
    {
      // rayfire goes through vertex
      
      // get all elems that share this vertex
      std::set<const libMesh::Elem*> elem_set;
      cur_elem->find_point_neighbors(end_point,elem_set);
      std::set<const libMesh::Elem *>::const_iterator       it  = elem_set.begin();
      const std::set<const libMesh::Elem *>::const_iterator end = elem_set.end();

      // iterate over each elem
      for (; it != end; ++it)
      {
        const libMesh::Elem* elem = *it;
        
        if (elem == cur_elem) // skip the current elem
          continue;
        
        // move a little bit along the rayfire
        // and see if we are in the elem
        libMesh::Real L = elem->hmin();
        L *= 0.1;
        
        
        
        // parametric representation of rayfire line       
        libMesh::Real x = end_point(0) + L*std::cos(_theta);
        libMesh::Real y = end_point(1) + L*std::sin(_theta);
        
        if ( elem->contains_point(libMesh::Point(x,y)) )
          return elem;
      }
      
    }
    else
    {
      // not a vertex, so just get the elem on that side
      return cur_elem->neighbor(side);
    }
    
    libmesh_error_msg("We shouldn't be here...");    
    return NULL;
  }


  bool RayfireMesh::_newton_solver(libMesh::Point& initial_point, const libMesh::Elem* edge_elem, libMesh::Point* intersection_point)
  {
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
    libMesh::Real d_ksi;
    
    // tan(theta) is the slope, so precompute since it is used repeatedly  
    libMesh::Real tan_theta = std::tan(_theta);
    
    // Initial guess is center of the edge_elem
    libMesh::Real ksi = 0.0;
    
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
                                                         ksi);
                                                        
        dphi[i] = libMesh::FE<1,libMesh::LAGRANGE>::shape_deriv(edge_elem->type(),
                                                                edge_elem->default_order(),
                                                                i,
                                                                0, // const unsigned int libmesh_dbg_varj
                                                                ksi);
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
            
      d_ksi = f/df;
            
      if(std::abs(d_ksi) < libMesh::TOLERANCE)
      {
        // convergence
        (*intersection_point)(0) = X;
        (*intersection_point)(1) = Y;
        return true;         
      }
      else
        ksi -= d_ksi;
        
    } // for it  
        
    // no convergence
    return false;
  }

} //namespace GRINS
