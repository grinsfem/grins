// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014 Paul T. Bauman, Roy H. Stogner
// Copyright (C) 2010-2013 The PECOS Development Team                        

// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation; either version 2.1
// of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

// system
#include <limits>
#include <cmath>

// libmesh
#include "libmesh/libmesh_base.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_output.h" // for MeshSerializer
#include "libmesh/equation_systems.h"
#include "libmesh/explicit_system.h"
#include "libmesh/boundary_mesh.h"
#include "libmesh/node.h"
#include "libmesh/elem.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/dense_vector.h"
#include "libmesh/fe_base.h"
#include "libmesh/dof_map.h"
#include "libmesh/point.h"

// local
#include "grins/distance_function.h"


// anonymous namespace for implementation details -
// gives file scope and prevents symbol clashing at link time
// in case there is an e.g. "Line" declared somewhere else too.
namespace
{
  using namespace libMesh;

  //***************************************************
  // A helper data structure. Only used internally.
  //***************************************************
  struct Line
  {
    Point p0;
    Point p1;
  };



  //***************************************************
  // Some helper functions.  Only used internally.
  //***************************************************

  //---------------------------------------------------
  // Project pt onto line.
  //
  // Returns xi, the location of projected point in
  // the reference space of the line.
  //
  template<int dim>
  void ProjectPointToLine( Point pt, Line line, Real& xi )
  {
    libmesh_assert( (dim==2) || (dim==3) );

    Point& p0 = line.p0;
    Point& p1 = line.p1;

    Real dx[dim];
    Real xavg[dim];
    Real length2=0;
    for ( unsigned int ii=0; ii<dim; ii++ ) {
      dx[ii] = p1(ii) - p0(ii);
      length2 += dx[ii]*dx[ii];
      xavg[ii] = 0.5*( p1(ii) + p0(ii) );
    }

    // xi is the point in reference space that minimizes
    // the distance between the point pt and the line
    // defined by the points in the Line struct
    xi = 0.0;
    for ( unsigned int ii=0; ii<dim; ii++ ) {
      xi += 2.0*( pt(ii) - xavg[ii] )*dx[ii];
    }
    xi /= length2;
  }

  //---------------------------------------------------
  // Compute distance from pt to line segment
  //
  template<int dim>
  void
  DistanceToSegment( Point pt, Line line, Real& distance )
  {
    libmesh_assert( (dim==2) || (dim==3) );

    Real xi;
    ProjectPointToLine<dim> (pt, line, xi);

    // if xi is outside edge, bring to nearest node
    if( xi < -1.0 ) xi = -1.0;
    if( xi >  1.0 ) xi =  1.0;

    // interpolate between p0 and p1 to compute location in physical space
    Real xint[dim];
    for ( unsigned int ii=0; ii<dim; ii++ ) {
      xint[ii] = line.p0(ii)*(-0.5*(xi-1.0)) + line.p1(ii)*(0.5*(xi+1.0));
    }

    // compute distance
    distance = 0.0;
    for ( unsigned int ii=0; ii<dim; ii++ ) {
      distance += (pt(ii) - xint[ii])*(pt(ii) - xint[ii]);
    }
    distance = sqrt(distance);

    libmesh_assert( (std::isfinite(distance)) && (distance>=0.0) );
  }

} // end anonymous namespace



namespace GRINS {

//***************************************************
// DistanceFunction class functions
//***************************************************

//---------------------------------------------------
// Constructor
//
//DistanceFunction::DistanceFunction (EquationSystems& equation_systems, const UnstructuredMesh& boundary_mesh)
  DistanceFunction::DistanceFunction (libMesh::EquationSystems &es_in, const libMesh::UnstructuredMesh &bm_in):
  _equation_systems (es_in),
  _boundary_mesh    (bm_in),
  _dist_fe          (libMesh::FEBase::build(_equation_systems.get_mesh().mesh_dimension(), libMesh::FEType(libMesh::FIRST, libMesh::LAGRANGE)))
{
  // Ensure that libmesh is ready to roll
  libmesh_assert(libMesh::initialized());

  // Add distance function system
  _equation_systems.add_system<libMesh::System>("distance_function");

  // Get reference to distance function system we just added
  libMesh::System& sys = _equation_systems.get_system<libMesh::System>("distance_function");

  // Add distance function variable
  sys.add_variable("distance", libMesh::FIRST);

  // Attach initialization function
  //sys.attach_init_object(*this);

}

//---------------------------------------------------
// Compute distance function
//
void DistanceFunction::initialize ()
{
  // Call the compute function
  this->compute();
  return;
}


//---------------------------------------------------
// Compute distance from input node to boundary_mesh
//
libMesh::Real DistanceFunction::node_to_boundary (const libMesh::Node* node)
{
  // Ensure that node is not NULL
  libmesh_assert( node != NULL );

  // Initialize distance to infinity
  libMesh::Real distance = std::numeric_limits<libMesh::Real>::infinity();

  // Get dimension
  const unsigned int dim = _equation_systems.get_mesh().mesh_dimension();
  libmesh_assert( (dim==2) || (dim==3) );

  // Get coordinates of node
  std::vector<libMesh::Real> xnode(dim);
  for( unsigned int ii=0; ii<dim; ii++ ) xnode[ii] = (*node)(ii);

  // This function will work on a distributed interior mesh, but won't
  // give correct results on a distributed boundary mesh.
  libmesh_assert(_boundary_mesh.is_serial());

  // First loop over the boundary mesh and compute the distance to the
  // nodes.  This will give us an idea of what elements to check more
  // closely.

  libMesh::MeshBase::const_element_iterator       el     = _boundary_mesh.active_elements_begin();
  const libMesh::MeshBase::const_element_iterator end_el = _boundary_mesh.active_elements_end();

  if (el==end_el)
    {
      std::cout << "There are no elements to iterate through!!!  Returning..." << std::endl;
      return distance;
    }

  const libMesh::Elem* elem_containing_min = NULL;
  libMesh::Real dmin = std::numeric_limits<Real>::infinity();

  for ( ; el != end_el; ++el) {

    const libMesh::Elem* belem = *el;

    for (unsigned int bnode=0; bnode<belem->n_nodes(); ++bnode)
      {
	const libMesh::Point& pbndry = belem->point(bnode);

	libMesh::Real dnode = 0.0;
	for (unsigned int idim=0; idim<dim; ++idim)
	  dnode += (xnode[idim] - pbndry(idim))*(xnode[idim] - pbndry(idim));

	dnode = std::sqrt(dnode);

	if (dnode<dmin)
	  {
	    dmin = dnode;
	    elem_containing_min = belem;
	  }
      }

  } // end element loop

  // grab the elements around the minimum
  std::set<const libMesh::Elem*> near_min_elems;
  elem_containing_min->find_point_neighbors (near_min_elems);

  // insert the element itself into the set... we must check it also
  near_min_elems.insert(elem_containing_min);

  std::set<const libMesh::Elem*>::iterator nm_el;

  // loop over the near_min_elems
  for (nm_el=near_min_elems.begin(); nm_el!=near_min_elems.end(); ++nm_el)
    {

      const libMesh::Elem* belem = *nm_el; //*el;

      // Ensure that elem defined by edge/face is linear
      libmesh_assert( belem->default_order() == FIRST );

      // Initialize distance to this edge/face to infinity
      libMesh::Real dedge = std::numeric_limits<libMesh::Real>::infinity();

      if ( dim==2 )
	{ // 2d

	  // For 2d mesh, boundary elems must be 1d.  Due to the way
	  // find_point_neighbors works, there may be some 2d elements
	  // in the near_min_elems set.  These elements are clearly
	  // not part of the boundary mesh, and thus we skip them.
	  //
	  if (belem->dim()!=1) continue;

	  libmesh_assert( belem->n_nodes() == 2 );
	  libmesh_assert( belem->type() == EDGE2 );

	  // Points defining the edge
	  const libMesh::Point& p0 = belem->point(0);
	  const libMesh::Point& p1 = belem->point(1);

	  Line line;
	  line.p0 = p0;
	  line.p1 = p1;

	  DistanceToSegment<2>(*node, line, dedge);

	  if ( dedge < distance ) distance = dedge;

	}
      else
	{ // 3d

	  // For 3d mesh, boundary elems must be 2d.  Due to the way
	  // find_point_neighbors works, there may be some 3d elements
	  // in the near_min_elems set.  These elements are clearly
	  // not part of the boundary mesh, and thus we skip them.
	  //
	  if (belem->dim()!=2) continue;

	  libmesh_assert( (belem->type() == TRI3) || (belem->type() == QUAD4) );

	  if ( belem->type() == TRI3 )
	    { // triangular boundary faces

	      // Find the point in the plane defined by the boundary nodes that
	      // minimizes the distance between the plane and the input node.
	      //
	      // Done in terms of the reference coordinates of the boundary element.
	      //
	      const libMesh::Point& p0 = belem->point(0);
	      const libMesh::Point& p1 = belem->point(1);
	      const libMesh::Point& p2 = belem->point(2);

	      const libMesh::Real x_xi = (p1(0) - p0(0)), x_et = (p2(0) - p0(0));
	      const libMesh::Real y_xi = (p1(1) - p0(1)), y_et = (p2(1) - p0(1));
	      const libMesh::Real z_xi = (p1(2) - p0(2)), z_et = (p2(2) - p0(2));

	      libMesh::Real A[2][2], b[2];

	      A[0][0] = x_xi*x_xi + y_xi*y_xi + z_xi*z_xi;
	      A[0][1] = x_xi*x_et + y_xi*y_et + z_xi*z_et;
	      A[1][0] = x_xi*x_et + y_xi*y_et + z_xi*z_et;
	      A[1][1] = x_et*x_et + y_et*y_et + z_et*z_et;

	      b[0] = (xnode[0] - p0(0))*x_xi + (xnode[1] - p0(1))*y_xi + (xnode[2] - p0(2))*z_xi;
	      b[1] = (xnode[0] - p0(0))*x_et + (xnode[1] - p0(1))*y_et + (xnode[2] - p0(2))*z_et;

	      const libMesh::Real detA = A[0][0]*A[1][1] - A[1][0]*A[0][1];
	      libmesh_assert( fabs(detA) > 0.0 ); // assert that A is not singular

	      const libMesh::Real xi = ( A[1][1]*b[0] - A[0][1]*b[1])/detA;
	      const libMesh::Real et = (-A[1][0]*b[0] + A[0][0]*b[1])/detA;


	      // If projection of node onto plane defined by boundary face
	      // (i.e., what we just computed) is not inside boundary face
	      // then we need to figure out closest point that is inside the
	      // boundary face
	      //
	      if ( (xi<0.0) || (et<0.0) || (et>1.0-xi) ) {

		libMesh::Real dtmp = std::numeric_limits<libMesh::Real>::infinity();

		// for each edge of boundary face, find distance from
		// input node to the segment defined by that edge
		//
		// the minimum of these distances minimizes the distance
		// to the node over the points in the boundary face
		//
		for ( unsigned int iedge=0; iedge<3; iedge++ ){

		  Line line;

		  // Get the endpoints of this edge
		  if      (iedge==0) { line.p0 = p1; line.p1 = p2; }
		  else if (iedge==1) { line.p0 = p0; line.p1 = p2; }
		  else   /*iedge==2*/{ line.p0 = p0; line.p1 = p1; }

		  DistanceToSegment<3>(*node, line, dtmp);

		  if( dtmp < dedge ) dedge = dtmp;

		}

	      } else { // projection is inside face, so we're good to go

		// Map from reference to physical space
		libMesh::Real xint[3];
		for ( unsigned int ii=0; ii<3; ii++ ) {
		  xint[ii] = p0(ii)*(1.0 - xi - et) + p1(ii)*xi + p2(ii)*et;
		}

		// compute distance
		dedge = 0.0;
		for ( unsigned int ii=0; ii<3; ii++ ) {
		  dedge += (xnode[ii] - xint[ii])*(xnode[ii] - xint[ii]);
		}
		dedge = sqrt(dedge);

	      }

	    }
	  else if ( belem->type() == QUAD4 )
	    {
	      //std::cout << "WARNING: this functionality is not well-tested.  Sorry." << std::endl;

	      libMesh::Real RTOL = 1e-10;
	      libMesh::Real ATOL = 1e-20;
	      const unsigned int ITER_MAX=1000;
	      unsigned int iter=0;

	      ComputeDistanceResidual res(belem, node);

	      ComputeDistanceJacobian jac;

	      libMesh::DenseVector<libMesh::Real> X(2); X(0) = X(1) = 0.0;
	      libMesh::DenseVector<libMesh::Real> dX(2);
	      libMesh::Real det;

	      libMesh::DenseVector<libMesh::Real> R(2);
	      libMesh::DenseMatrix<libMesh::Real> dRdX(2,2);
	      libMesh::DenseMatrix<libMesh::Real> dRdXinv(2,2);

	      // evaluate residual... maybe we don't have to iterate
	      res(X,R);

	      libMesh::Real R0 = R.l2_norm();

	      // iterate
	      while ( (R.l2_norm() > ATOL) && (R.l2_norm()/R0 > RTOL) && (iter<ITER_MAX) )
		{
		  // compute jacobian
		  jac(X, R, res, dRdX);

		  // invert jacobian
		  det = dRdX(0,0)*dRdX(1,1) - dRdX(0,1)*dRdX(1,0);

		  // protect against divide by zero
		  if(std::abs(det)<=1e-20)
		    {
		      std::cout << "WARNING: about to divide by zero." << std::endl;
		    }

		  dRdXinv(0,0) =  dRdX(1,1)/det;
		  dRdXinv(0,1) = -dRdX(0,1)/det;
		  dRdXinv(1,0) = -dRdX(1,0)/det;
		  dRdXinv(1,1) =  dRdX(0,0)/det;

		  // compute delta
		  dX(0) = -( dRdXinv(0,0)*R(0) + dRdXinv(0,1)*R(1) );
		  dX(1) = -( dRdXinv(1,0)*R(0) + dRdXinv(1,1)*R(1) );

		  // update
		  X += dX;

		  // if ( std::abs(X(0)) > 1.0 || std::abs(X(1)) > 1.0 )
		  //   {
		  //     std::cout << "WARNING: on iter = " << iter << ", xi is leaving the element!" << std::endl;
		  //   }

		  // recompute Residual
		  res(X,R);

		  // increment counter
		  iter++;
		}

	      // check that we converged
	      if (iter==ITER_MAX)
		{
		  // std::cerr << "Failed to converge distance function!!!" << std::endl;
		  // std::cerr << "Started with ||R|| = " << R0 << std::endl;
		  // std::cerr << "Finished " << iter << " iterations with ||R|| = " << R.l2_norm() << std::endl;
		  // std::cout << "Final location was xi = (" << X(0) << ", " << X(1) << ")." << std::endl;
		  // libmesh_error();

		  std::cout << "WARNING: Failed to converge distance function iteration.  Assuming that min is outside this face.  Use closest distance to edges to proceed." << std::endl;
		}
	      else
		{
		  //std::cout << "Converged!" << std::endl;
		}


	      if ( std::abs(X(0)) > 1.0 || std::abs(X(1)) > 1.0 )
		{
		  // converged to point outside the face... so min must
		  // be on edge, and luckily the edges are linear

		  libMesh::Real dtmp = std::numeric_limits<libMesh::Real>::infinity();

		  // for each edge of boundary face, find distance from
		  // input node to the segment defined by that edge
		  //
		  // the minimum of these distances minimizes the distance
		  // to the node over the points in the boundary face
		  //
		  for ( unsigned int iedge=0; iedge<4; iedge++ ){

		    Line line;

		    // Get the endpoints of this edge
		    if      (iedge==0) { line.p0 = belem->point(0); line.p1 = belem->point(1); }
		    else if (iedge==1) { line.p0 = belem->point(1); line.p1 = belem->point(2); }
		    else if (iedge==2) { line.p0 = belem->point(2); line.p1 = belem->point(3); }
		    else   /*iedge==3*/{ line.p0 = belem->point(3); line.p1 = belem->point(0); }

		    DistanceToSegment<3>(*node, line, dtmp);

		    if( dtmp < dedge ) dedge = dtmp;
		  }
		}
	      else
		{
		  // converged to point inside, so compute point in
		  // physical space and use it to evaluate the distance

		  // assuming first order lagrange basis here
		  libMesh::AutoPtr<libMesh::FEBase> fe( libMesh::FEBase::build(2, libMesh::FEType(libMesh::FIRST, libMesh::LAGRANGE)) );

		  std::vector<libMesh::Point> xi(1);
		  xi[0](0) = X(0);
		  xi[0](1) = X(1);

		  // grab basis functions (evaluated at qpts)
		  const std::vector<std::vector<libMesh::Real> > &basis = fe->get_phi();

		  // reinitialize finite element data at xi
		  fe->reinit(belem, &xi);

		  // interpolate location
		  libMesh::DenseVector<libMesh::Real> xx(3);
		  xx.zero();

		  for (unsigned int inode=0; inode<belem->n_nodes(); ++inode)
		    for (unsigned int idim=0; idim<3; ++idim)
		      xx(idim) += belem->point(inode)(idim) * basis[inode][0];

		  // compute distance
		  dedge = 0.0;
		  for ( unsigned int ii=0; ii<3; ii++ ) {
		    dedge += (xnode[ii] - xx(ii))*(xnode[ii] - xx(ii));
		  }
		  dedge = sqrt(dedge);

		}

	    }
	  else
	    { // higher-order faces not supported yet
	      std::cout << "My type is " << belem->type() << std::endl;
	      libmesh_not_implemented();
	    }

	  if ( dedge < distance ) distance = dedge;

	} // end if(2d)

    } // end loop over elements

  // make sure we are returning valid distance---i.e., 0 <= distance < infinity
  // But this may fail if we try to get the distance on a mesh with no
  // boundary
  //  libmesh_assert( (std::isfinite(distance)) && (distance>=0.0) );

  return distance;
}


//---------------------------------------------------
// Fill distance function data at mesh nodes
//
void DistanceFunction::compute ()
{
  // Get mesh
  const libMesh::MeshBase& mesh = _equation_systems.get_mesh();

  // Get reference to system and system number
  libMesh::System& system = _equation_systems.get_system<libMesh::System>("distance_function");
  const unsigned int sys_num = system.number();

  // The boundary mesh needs to all be on this processor for us to
  // calculate a correct distance function.  Since we don't need it to
  // be serial afterwards, we use a temporary serializer.
  {
  libMesh::MeshSerializer serialize(const_cast<libMesh::UnstructuredMesh&>(_boundary_mesh));

  // Loop over nodes in mesh
  libMesh::MeshBase::const_node_iterator node_it  = mesh.local_nodes_begin();
  const libMesh::MeshBase::const_node_iterator node_end = mesh.local_nodes_end();

  for ( ; node_it != node_end; ++node_it) {

    // Grab node
    const libMesh::Node* node = *node_it;

    // Compute distance to nearest point in boundary_mesh
    const libMesh::Real distance = DistanceFunction::node_to_boundary (node);

    // Stuff data into appropriate place in the system solution
    const unsigned int dof = node->dof_number(sys_num,0,0);
    system.solution->set (dof, distance);

  } // end loop over nodes

  } // end boundary mesh serialization

  system.solution->close();
  system.update();

  std::cout << "Distance Function computed." << std::endl;
}


//---------------------------------------------------
// Interpolate nodal data
//
libMesh::AutoPtr< libMesh::DenseVector<libMesh::Real> >
DistanceFunction::interpolate (const libMesh::Elem* elem, const std::vector<libMesh::Point>& qpts) const
{
  libmesh_assert( elem != NULL );    // can't interpolate in NULL elem
  libmesh_assert( qpts.size() > 0 ); // can't interpolate if no points requested

  // grab basis functions (evaluated at qpts)
  const std::vector<std::vector<libMesh::Real> > &phi = _dist_fe->get_phi();

  // reinitialize finite element data at qpts
  _dist_fe->reinit(elem, &qpts);

  // number of basis functions
  const unsigned int n_dofs = phi.size();

  // number of points
  const unsigned int n_pts = qpts.size();

  // instantiate auto_ptr to dense vector to hold results
  libMesh::AutoPtr< DenseVector<libMesh::Real> > ap( new libMesh::DenseVector<libMesh::Real>(qpts.size()) );
  (*ap).zero();

  // pull off distance function at nodes on this element
  libMesh::System& sys = _equation_systems.get_system<libMesh::System>("distance_function");
  const libMesh::DofMap& dof_map = sys.get_dof_map();

  std::vector<unsigned int> dof_ind;
  dof_map.dof_indices(elem, dof_ind);

  libMesh::DenseVector<libMesh::Real> nodal_dist;
  nodal_dist.resize(n_dofs);
  //dof_map.extract_local_vector( *(sys.solution), dof_ind, nodal_dist);
  dof_map.extract_local_vector( *(sys.current_local_solution), dof_ind, nodal_dist);

  for ( unsigned int idof=0; idof<n_dofs; idof++ ) {
    for ( unsigned int iqpt=0; iqpt<n_pts; iqpt++ ) {
      (*ap)(iqpt) += nodal_dist(idof) * phi[idof][iqpt];
    }
  }

  return ap;
}


} // end namespace GRINS
