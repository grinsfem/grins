//-----------------------------------------------------------------------bl-
//
// Fully Implicit Navier-Stokes (FIN-S)
//
// Copyright (C) 2002-2012 Benjamin S. Kirk, Roy H. Stogner,
//                         Todd A. Oliver, Paul T. Bauman
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
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
//
// $Id: distance_function.h 38624 2013-04-11 16:11:59Z benkirk $
//
//--------------------------------------------------------------------------

#ifndef GRINS_DISTANCE_FUNCTION_H
#define GRINS_DISTANCE_FUNCTION_H

// system
//#include <limits>

// libmesh
#include "libmesh/libmesh_base.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/fe_base.h"
#include "libmesh/system.h"

// FIN-S
#include "fins_config.h"

// Forward Declarations
namespace libMesh {
  class EquationSystems;
  class BoundaryMesh;
  class Node;
  class Elem;
}


namespace GRINS {

// This class provides the functionality to compute the distance to
// the nearest no slip wall boundary.
class DistanceFunction : public System::Initialization
{
public:

  /**
   * Constructor
   */
  DistanceFunction (EquationSystems &es_in, const UnstructuredMesh &bm_in);

  /**
   * Destructor
   */
  ~DistanceFunction () {};

  /**
   * Initializes the distance
   */
  virtual void initialize ();

  /**
   * Compute distance from input node to boundary_mesh
   */
  Real node_to_boundary (const Node* node);

  /**
   * Initialize "distance_function" equation system by computing
   * distance from each mesh node to the nearest point in boundary_mesh
   */
  void compute ();

  /**
   * Interpolate distance function to points qpts (in reference space) for element *elem
   */
  AutoPtr< DenseVector<Real> > interpolate (const Elem* elem, const std::vector<Point>& qts) const;



private:

  /**
   * Pointer to EquationSystems object
   */
  EquationSystems &_equation_systems;

  /**
   * Pointer to boundary mesh object
   */
  const UnstructuredMesh &_boundary_mesh;

  /**
   * Finite element to use for interpolation of distance.
   * For internal use only, so don't provide any access
   * to the outside world.
   *
   * Currently type is hardcoded to first order, Lagrange
   * (see constructor), but this could be easily changed.
   */
  AutoPtr<FEBase> _dist_fe;

};

/**
 * This class provides the functionality to compute finite difference
 * Jacobians required when computing the distance function
 */
class ComputeDistanceJacobian
{
public:

  // ctor
  ComputeDistanceJacobian(){}

  // dtor
  ~ComputeDistanceJacobian(){}

  /**
   * Finite-differenced Jacobian approximation.
   */
  template <class f>
  void operator()(const DenseVector<Real> &U,
		  const DenseVector<Real> & F0,
		  f &F,
		  DenseMatrix<Real> &dFdU)
  {

    Up = U;
    Fp.resize(U.size());

    dFdU.resize(U.size(),U.size());

    // compute dF_i/dU_j.
    for (unsigned int j=0; j<U.size(); j++)
      {
	// tell the user's function which component we are perturbing
	// at this iteration.  this information may be used to avoid
	// repeated calculations of values which are known to be
	// unaffected for this perturbation.
	//F.set_perturbed_component(j);

	// define the pertubation for this component
	const Real
	  pert  = 4.e-8*(std::max (std::abs(U(j)), 1.e-3)); // U can be 0, pertubation cannot

	Up(j) += pert;

	const Real                     // note this difference may not strictly be pert
	  invpert = 1./(Up(j) - U(j)); // due to truncation error in the preceeding sum.

	// evaluate F at the perturbed state
	F(Up,Fp);

	// remove perturbation for next evaluation
	Up(j) = U(j);

	for (unsigned int i=0; i<U.size(); i++)
	  dFdU(i,j) = (Fp(i) - F0(i))*invpert;
      }
  }

private:

  // work vectors
  DenseVector<Real>    Up, Fp, Um, Fm;
};



/**
 * This class provides the functionality to compute the residual
 * required when computing the distance function.  The residual here
 * corresponds to the derivative of the distance between the field
 * point and the boundary point wrt the boundary parameterization.
 * When this residual is zero, the distance is minimized, and we have
 * found the point we're looking for.
 */
class ComputeDistanceResidual
{
public:

  // ctor
  ComputeDistanceResidual(const Elem* belem, const Point* point) :
    _belem(*belem),
    _dim(_belem.dim()+1),
    _p(*point),
    _fe (FEBase::build(_belem.dim(), FEType(FIRST, LAGRANGE))),
    fe  (_fe.get())
  {
    // only supporting QUAD4 this way for now
    libmesh_assert(_belem.type()==QUAD4);
  }

  // dtor
  ~ComputeDistanceResidual(){}

  // Calculate the residual
  void operator()(const DenseVector<Real> &U,
		        DenseVector<Real> &F)
  {

    libmesh_assert(U.size()==_dim-1);
    libmesh_assert(F.size()==_dim-1);

    DenseVector<Real> xx(_dim);
    xx.zero();

    DenseMatrix<Real> xx_U(_dim, _dim-1);
    xx_U.zero();

    std::vector<Point> xi(1);

    xi[0](0) = U(0);
    xi[0](1) = U(1);

    // grab basis functions (evaluated at qpts)
    const std::vector<std::vector<Real> > &basis = fe->get_phi();
    const std::vector<std::vector<RealGradient> > &dbasis = fe->get_dphi();

    // reinitialize finite element data at xi
    fe->reinit(&_belem, &xi);

    // interpolate location and derivatives
    for (unsigned int inode=0; inode<_belem.n_nodes(); ++inode)
      for (unsigned int idim=0; idim<_dim; ++idim)
	{
	  xx(idim) += _belem.point(inode)(idim) * basis[inode][0];

	  for (unsigned int jdim=0; jdim<_dim-1; ++jdim)
	    xx_U(idim, jdim) += _belem.point(inode)(idim) * dbasis[inode][0](jdim);
	}


    // form the residual
    F.zero();
    for (unsigned int ires=0; ires<_dim-1; ++ires)
      for (unsigned int idim=0; idim<_dim; ++idim)
	F(ires) += ( xx(idim) - _p(idim) ) * xx_U(idim, ires);


  }

private:

  // Reference to boundary element
  const Elem& _belem;
  const unsigned int _dim;
  const Point& _p;

  AutoPtr<FEBase> _fe;
  FEBase *fe;

  // work vectors
  DenseVector<Real>    Up, Fp, Um, Fm;
};

} // end namespace FINS

#endif // __distance_function_h__
