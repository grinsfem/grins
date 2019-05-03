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

#ifndef GRINS_NONLINEAR_ELASTICITY_WEAK_FORM_H
#define GRINS_NONLINEAR_ELASTICITY_WEAK_FORM_H

namespace GRINS
{
  //! Class to encapsulate residual evaluations for hyperelasticity
  /*!
   * Several Physics classes will reuse these terms so let's encapsulate them
   * in an object. These will sit inside dof loops so we want them to be inlined.
   * This object supports both 2D (plane strain) and 3D evaluations. The functions
   * are differentiated by the number of residual vectors passed to the function (2 or 3)
   * and similarly for the Jacobian.
   */
  template<typename StrainEnergy>
  class HyperelasticityWeakForm
  {
  public:

    HyperelasticityWeakForm() = default;
    ~HyperelasticityWeakForm() = default;

    //! 2D version, consistent with plane strain assumption
    void evaluate_internal_stress_residual( const libMesh::Tensor & P,
                                            const libMesh::RealGradient & dphi_i_times_JxW,
                                            libMesh::Number & Fu,
                                            libMesh::Number & Fv ) const;

    //! Full 3D version
    void evaluate_internal_stress_residual( const libMesh::Tensor & P,
                                            const libMesh::RealGradient & dphi_i_times_JxW,
                                            libMesh::Number & Fu,
                                            libMesh::Number & Fv,
                                            libMesh::Number & Fw ) const;


    //! 2D version, consistent with plane strain assumption
    void evaluate_internal_stress_jacobian( const libMesh::Tensor & S,
                                            const libMesh::Tensor & F,
                                            const libMesh::RealGradient & dphi_i_times_JxW,
                                            const libMesh::RealGradient & dphi_j,
                                            const CartesianHyperlasticity<StrainEnergy> & stress_law,
                                            libMesh::Number & Kuu,
                                            libMesh::Number & Kuv,
                                            libMesh::Number & Kvu,
                                            libMesh::Number & Kvv  ) const;

    //! Full 3D version
    void evaluate_internal_stress_jacobian( const libMesh::Tensor & S,
                                            const libMesh::Tensor & F,
                                            const libMesh::RealGradient & dphi_i_times_JxW,
                                            const libMesh::RealGradient & dphi_j,
                                            const CartesianHyperlasticity<StrainEnergy> & stress_law,
                                            libMesh::Number & Kuu,
                                            libMesh::Number & Kuv,
                                            libMesh::Number & Kuw,
                                            libMesh::Number & Kvu,
                                            libMesh::Number & Kvv,
                                            libMesh::Number & Kvw,
                                            libMesh::Number & Kwu,
                                            libMesh::Number & Kwv,
                                            libMesh::Number & Kww ) const;

  };

  template<typename StrainEnergy>
  inline
  void HyperelasticityWeakForm<StrainEnergy>::evaluate_internal_stress_residual
  ( const libMesh::Tensor & P,
    const libMesh::RealGradient & dphi_i_times_JxW,
    libMesh::Number & Fu,
    libMesh::Number & Fv ) const
  {
    for( int alpha = 0; alpha < 2 /*dimension*/; alpha++)
      {
        libMesh::Real c = dphi_i_times_JxW(alpha);
        Fu += P(0,alpha)*c;
        Fv += P(1,alpha)*c;
      }
  }

  template<typename StrainEnergy>
  inline
  void HyperelasticityWeakForm<StrainEnergy>::evaluate_internal_stress_residual
  ( const libMesh::Tensor & P,
    const libMesh::RealGradient & dphi_i_times_JxW,
    libMesh::Number & Fu,
    libMesh::Number & Fv,
    libMesh::Number & Fw ) const
  {
    for( int alpha = 0; alpha < 3 /*dimension*/; alpha++)
      {
        libMesh::Real c = dphi_i_times_JxW(alpha);
        Fu += P(0,alpha)*c;
        Fv += P(1,alpha)*c;
        Fw += P(2,alpha)*c;
      }
  }

  template<typename StrainEnergy>
  inline
  void HyperelasticityWeakForm<StrainEnergy>::evaluate_internal_stress_jacobian
  ( const libMesh::Tensor & S,
    const libMesh::Tensor & F,
    const libMesh::RealGradient & dphi_i_times_JxW,
    const libMesh::RealGradient & dphi_j,
    const CartesianHyperlasticity<StrainEnergy> & stress_law,
    libMesh::Number & Kuu,
    libMesh::Number & Kuv,
    libMesh::Number & Kvu,
    libMesh::Number & Kvv ) const
  {
    // Compute the  derivative term without the elasticity tensor
    libMesh::Number term1 = dphi_j*(S*dphi_i_times_JxW);

    Kuu += term1;
    Kvv += term1;

    // Now compute the derivative term with the elasticity tensor
    const int dim = 2;
    for( int I = 0; I < dim; I++)
      for( int J = 0; J < dim; J++)
        for( int K = 0; K < dim; K++)
          for( int L = 0; L < dim; L++)
            {
              libMesh::Number Cijkl = stress_law.elasticity_tensor(I,J,K,L);

              libMesh::Real c0 = dphi_i_times_JxW(J)*Cijkl*dphi_j(L);

              Kuu += F(0,I)*F(0,K)*c0;
              Kuv += F(0,I)*F(1,K)*c0;
              Kvu += F(1,I)*F(0,K)*c0;
              Kvv += F(1,I)*F(1,K)*c0;
            }
  }


  template<typename StrainEnergy>
  inline
  void HyperelasticityWeakForm<StrainEnergy>::evaluate_internal_stress_jacobian
  ( const libMesh::Tensor & S,
    const libMesh::Tensor & F,
    const libMesh::RealGradient & dphi_i_times_JxW,
    const libMesh::RealGradient & dphi_j,
    const CartesianHyperlasticity<StrainEnergy> & stress_law,
    libMesh::Number & Kuu,
    libMesh::Number & Kuv,
    libMesh::Number & Kuw,
    libMesh::Number & Kvu,
    libMesh::Number & Kvv,
    libMesh::Number & Kvw,
    libMesh::Number & Kwu,
    libMesh::Number & Kwv,
    libMesh::Number & Kww ) const
  {
    // Compute the  derivative term without the elasticity tensor
    libMesh::Number term1 = dphi_j*(S*dphi_i_times_JxW);

    Kuu += term1;
    Kvv += term1;
    Kww += term1;

    // Now compute the derivative term with the elasticity tensor
    const int dim = 3;
    for( int I = 0; I < dim; I++)
      for( int J = 0; J < dim; J++)
        for( int K = 0; K < dim; K++)
          for( int L = 0; L < dim; L++)
            {
              libMesh::Number Cijkl = stress_law.elasticity_tensor(I,J,K,L);

              libMesh::Real c0 = dphi_i_times_JxW(J)*Cijkl*dphi_j(L);

              Kuu += F(0,I)*F(0,K)*c0;
              Kuv += F(0,I)*F(1,K)*c0;
              Kuw += F(0,I)*F(2,K)*c0;
              Kvu += F(1,I)*F(0,K)*c0;
              Kvv += F(1,I)*F(1,K)*c0;
              Kvw += F(1,I)*F(2,K)*c0;
              Kwu += F(2,I)*F(0,K)*c0;
              Kwv += F(2,I)*F(1,K)*c0;
              Kww += F(2,I)*F(2,K)*c0;
            }
  }

} // end namespace GRINS

#endif // GRINS_NONLINEAR_ELASTICITY_WEAK_FORM_H
