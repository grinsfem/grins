//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2014 Paul T. Bauman, Roy H. Stogner
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

#ifndef GRINS_HYPERELASTICITY_H
#define GRINS_HYPERELASTICITY_H

namespace GRINS
{
  template <typename StrainEnergy>
  class Hyperelasticity : public StressStrainLaw<Hyperelasticity<StrainEnergy> >
  {
  public:
    Hyperelasticity( const GetPot& input );
    virtual ~Hyperelasticity();

    // So we can make implementation private
    friend class StressStrainLaw<Hyperelasticity<StrainEnergy> >;

  private:

    void compute_stress_imp( const libMesh::TensorValue<libMesh::Real>& g_contra,
                             const libMesh::TensorValue<libMesh::Real>& G_contra,
                             const libMesh::TensorValue<libMesh::Real>& G_cov,
                             const libMesh::TensorValue<libMesh::Real>& strain,
                             unsigned int dim,
                             libMesh::TensorValue<libMesh::Real>& stress );

    StrainEnergy _W;

  };

} // end namespace GRINS

#endif // GRINS_HYPERELASTICITY_H
