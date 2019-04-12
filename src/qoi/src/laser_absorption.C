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

// GRINS
#include "grins/laser_absorption.h"
#include "grins/fem_function_and_derivative_base.h"
#include "grins/spectroscopic_transmission.h"
#include "grins/laser_intensity_profile_base.h"
#include "grins/math_constants.h"

// libMesh
#include "libmesh/edge_edge2.h"
#include "libmesh/face_quad9.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/serial_mesh.h"

namespace GRINS
{
  // 2D constructor
  LaserAbsorption::LaserAbsorption( const std::shared_ptr<FEMFunctionAndDerivativeBase<libMesh::Real>> & absorb,
                                    const libMesh::Point & top_origin, const libMesh::Point & centerline_origin,
                                    const libMesh::Point & bottom_origin,
                                    libMesh::Real theta, unsigned int n_quadrature_points,
                                    std::shared_ptr<LaserIntensityProfileBase> intensity_profile,    
                                    const std::string & qoi_name)
    : MultiQoIBase(qoi_name),
      _intensity_profile(intensity_profile)
  {
    // create an EDGE2 elem to represent the start of the laser beam
    std::shared_ptr<libMesh::Elem> elem( new libMesh::Edge2() );
    elem->set_node(0) = new libMesh::Node(bottom_origin);
    elem->set_node(1) = new libMesh::Node(top_origin);

    // now use QGauss to identify the quadratures weights and points on this "laser" elem
    libMesh::Order order = (libMesh::Order)(2*n_quadrature_points - 1);
    libMesh::QGauss qbase(elem->dim(),order);
    qbase.init(elem->type(),order);

    std::unique_ptr< libMesh::FEGenericBase<libMesh::Real> > fe = libMesh::FEGenericBase<libMesh::Real>::build(elem->dim(),libMesh::FEType(libMesh::FIRST,libMesh::LAGRANGE));
    fe->attach_quadrature_rule( &qbase );

    const std::vector<libMesh::Point> & quadrature_xyz = fe->get_xyz();

    fe->reinit(elem.get());

    _quadrature_weights = qbase.get_weights();

    // init the intensity profile object
    _intensity_profile->init(quadrature_xyz,centerline_origin);

    // create an internal SpectroscopicTramsmission object at each quadrature node
    std::shared_ptr<RayfireMesh> rayfire;
    for(unsigned int p = 0; p < quadrature_xyz.size(); ++p)
      {
        libMesh::Point origin = quadrature_xyz[p];
        rayfire.reset( new RayfireMesh(origin,theta) );
        SpectroscopicTransmission spec(absorb,rayfire,qoi_name,false);

        this->add_qoi(spec);
      }

  }

  // 3D constructor
  LaserAbsorption::LaserAbsorption(const std::shared_ptr<FEMFunctionAndDerivativeBase<libMesh::Real>> & absorb,
                    const libMesh::Point & top_origin, const libMesh::Point & centerline_origin,
                    const libMesh::Point & bottom_origin,
                    libMesh::Real theta, libMesh::Real phi, unsigned int n_quadrature_points,
                    std::shared_ptr<LaserIntensityProfileBase> intensity_profile,   
                    const std::string & qoi_name)
    : MultiQoIBase(qoi_name),
      _intensity_profile(intensity_profile)
  {
    libMesh::Point P0(centerline_origin);
    libMesh::Point P1(top_origin);
    libMesh::Point P2(bottom_origin);

    libMesh::Point P01(P1-P0);
    libMesh::Point P02(P2-P0);
    libMesh::Point n = P01.cross(P02);

    libMesh::Point w = n.cross(P01);

    libMesh::Real a = std::acos( (P01*P02) / (P01.norm() * P02.norm()));
    if ( (a < libMesh::TOLERANCE) || ( std::abs(a - Constants::pi) < libMesh::TOLERANCE) )
      libmesh_error_msg("Error in LaserAbsorption: top_origin, centerline_origin, and bottom_origin points cannot be colinear");

    libMesh::Real radius = P01.norm();

    // make sure the radius is constant
    libmesh_assert_less(std::abs(radius-P02.norm()),libMesh::TOLERANCE); 

    libMesh::Point p = P01;

    std::shared_ptr<libMesh::Elem> elem( new libMesh::Quad9() );

    elem->set_node(0) = new libMesh::Node(p);
    elem->get_node(0)->set_id(0);

    for (unsigned int s = 1; s < 8; ++s)
      {
        libMesh::Real angle = Constants::pi/2.0;

        if (s == 4)
          angle += (Constants::pi/4.0);

        libMesh::Real x1 = std::cos(angle)/radius;
        libMesh::Real x2 = std::sin(angle)/(w.norm());

        libMesh::Point node = radius*(x1*p + x2*w);

        elem->set_node(s) = new libMesh::Node(node);
        elem->get_node(s)->set_id(s);

        p = node;
      }

    elem->set_node(8) = new libMesh::Node(centerline_origin);
    elem->get_node(8)->set_id(8);

    // need a dummy ID to get through several asserts
    elem->set_id(0);

    // now use QGauss to identify the quadratures weights and points on this "laser" elem
    libMesh::Order order = (libMesh::Order)(2*n_quadrature_points - 1);
    libMesh::QGauss qbase(elem->dim(),order);
    qbase.init(elem->type(),order);

    std::unique_ptr< libMesh::FEGenericBase<libMesh::Real> > fe = libMesh::FEGenericBase<libMesh::Real>::build(elem->dim(),libMesh::FEType(libMesh::FIRST,libMesh::LAGRANGE));
    fe->attach_quadrature_rule( &qbase );

    const std::vector<libMesh::Point> & quadrature_xyz = fe->get_xyz();

    fe->reinit(elem.get());

    _quadrature_weights = qbase.get_weights();

    // init the intensity profile object
    _intensity_profile->init(quadrature_xyz,centerline_origin);

    // create an internal SpectroscopicTramsmission object at each quadrature node
    std::shared_ptr<RayfireMesh> rayfire;
    for(unsigned int p = 0; p < quadrature_xyz.size(); ++p)
      {
        libMesh::Point origin = quadrature_xyz[p];
        rayfire.reset( new RayfireMesh(origin,theta,phi) );
        SpectroscopicTransmission spec(absorb,rayfire,qoi_name,false);

        this->add_qoi(spec);
      }
  }

  QoIBase * LaserAbsorption::clone() const
  {
    return new LaserAbsorption(*this);
  }

  void LaserAbsorption::element_qoi(AssemblyContext & context, const unsigned int qoi_index)
  {
    // perhaps a bit hacky, but since the Context doesn't know about the internal SpectroscopicTransmission
    // QoIs, we have to use the Context as a way to pass the QoI contributions to the class
    // and store them internally in _qoi_vals
    libMesh::Number & qoi = context.get_qois()[qoi_index];
    for(unsigned int q = 0; q < this->n_qois(); ++q)
      {
        qoi = 0.0;

        _qois[q]->element_qoi(context,qoi_index);

        _qoi_vals[q] += qoi;
      }

    qoi = 0.0;

  }

  void LaserAbsorption::element_qoi_derivative(AssemblyContext & /*context*/, const unsigned int /*qoi_index*/)
  {
    // TODO
    libmesh_not_implemented();
  }

  void LaserAbsorption::parallel_op(const libMesh::Parallel::Communicator & communicator,
                                    libMesh::Number & sys_qoi, libMesh::Number & /*local_qoi*/)
  {
    // here we don't care about the local_qoi because we keep track of all the
    // QoI values internally. So we just do our absorption calculation
    // and set the sys_qoi and _qoi_value at the end

    // we first need to sum the _qoi_vals from all processors
    // and do the exponential to get the actual transmission value
    for (unsigned int i = 0; i < this->n_qois(); ++i)
      _qois[i]->parallel_op(communicator,_qoi_vals[i],_qoi_vals[i]);

    // calculate the initial laser intensity, Io*
    libMesh::Real initial = 0.0;
    for (unsigned int i = 0; i < this->n_qois(); ++i)
      {
        libMesh::Real intensity = _intensity_profile->intensity(i);
        libMesh::Real gq_weight = _quadrature_weights[i];

        initial += gq_weight*intensity;
      }

    // calculate the final laser intensity, If*
    libMesh::Real final = 0.0;
    for (unsigned int i = 0; i < this->n_qois(); ++i)
      {
        libMesh::Real gq_weight = _quadrature_weights[i];
        libMesh::Real intensity = _intensity_profile->intensity(i);
        libMesh::Real qoi_val = _qoi_vals[i];

        final += gq_weight*intensity*qoi_val;
      }

    // this is the actual QoI value we are looking for: the total absorbed laser intensity
    libMesh::Real absorbed_laser_intensity = 1.0 - final/initial;

    sys_qoi = absorbed_laser_intensity;
    QoIBase::_qoi_value = sys_qoi;
  }

}
