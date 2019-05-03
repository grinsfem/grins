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

#ifndef GRINS_INCOMPRESSIBLE_HYPERELASTICITY_WEAK_FORM_H
#define GRINS_INCOMPRESSIBLE_HYPERELASTICITY_WEAK_FORM_H

#include "grins/hyperelasticity_weak_form.h"

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
  class IncompressibleHyperelasticityWeakForm : public HyperelasticityWeakForm<StrainEnergy>
  {
  public:

    IncompressibleHyperelasticityWeakForm() = default;
    virtual ~IncompressibleHyperelasticityWeakForm() = default;

    //! 2D version, consistent with plane strain assumption
    void evaluate_pressure_stress_residual( libMesh::Number J,
                                            libMesh::Number press,
                                            const libMesh::Tensor & F_times_Cinv,
                                            const libMesh::RealGradient & dphi_i_times_JxW,
                                            libMesh::Number & Fu,
                                            libMesh::Number & Fv ) const;

    //! Full 3D version
    void evaluate_pressure_stress_residual( libMesh::Number J,
                                            libMesh::Number press,
                                            const libMesh::Tensor & F_times_Cinv,
                                            const libMesh::RealGradient & dphi_i_times_JxW,
                                            libMesh::Number & Fu,
                                            libMesh::Number & Fv,
                                            libMesh::Number & Fw ) const;

    void evaluate_pressure_stress_displacement_jacobian( libMesh::Number J,
                                                         libMesh::Number press,
                                                         const libMesh::Tensor & F,
                                                         const libMesh::Tensor & Cinv,
                                                         const libMesh::Tensor & F_times_Cinv,
                                                         const libMesh::RealGradient & dphi_i_times_JxW,
                                                         const libMesh::RealGradient & dphi_j,
                                                         libMesh::Number & Kuu,
                                                         libMesh::Number & Kuv,
                                                         libMesh::Number & Kvu,
                                                         libMesh::Number & Kvv ) const;

    void evaluate_pressure_stress_displacement_jacobian( libMesh::Number J,
                                                         libMesh::Number press,
                                                         const libMesh::Tensor & F,
                                                         const libMesh::Tensor & Cinv,
                                                         const libMesh::Tensor & F_times_Cinv,
                                                         const libMesh::RealGradient & dphi_i_times_JxW,
                                                         const libMesh::RealGradient & dphi_j,
                                                         libMesh::Number & Kuu,
                                                         libMesh::Number & Kuv,
                                                         libMesh::Number & Kuw,
                                                         libMesh::Number & Kvu,
                                                         libMesh::Number & Kvv,
                                                         libMesh::Number & Kvw,
                                                         libMesh::Number & Kwu,
                                                         libMesh::Number & Kwv,
                                                         libMesh::Number & Kww ) const;

    void evaluate_pressure_stress_pressure_jacobian( libMesh::Number J,
                                                     const libMesh::Tensor & F_times_Cinv,
                                                     const libMesh::Real & p_phi_j,
                                                     const libMesh::RealGradient & dphi_i_times_JxW,
                                                     libMesh::Number & Kup,
                                                     libMesh::Number & Kvp ) const;

    void evaluate_pressure_stress_pressure_jacobian( libMesh::Number J,
                                                     const libMesh::Tensor & F_times_Cinv,
                                                     const libMesh::Real & p_phi_j,
                                                     const libMesh::RealGradient & dphi_i_times_JxW,
                                                     libMesh::Number & Kup,
                                                     libMesh::Number & Kvp,
                                                     libMesh::Number & Kwp ) const;

    void evaluate_pressure_constraint_residual(libMesh::Number J,
                                               const libMesh::Real & phi_times_JxW,
                                               libMesh::Number & Fp) const;



  protected:

    //! Volumetric stored energy function
    /*!
     * We use the form U(J) = 1/2(ln(J)^2)
     */
    libMesh::Number U( libMesh::Number J ) const;

    //! First derivative of volumetric stored energy function
    libMesh::Number dU( libMesh::Number J ) const;

    //! Second derivative of volumetric stored energy function
    libMesh::Number d2U( libMesh::Number J ) const;

  };


  template<typename StrainEnergy>
  inline
  void IncompressibleHyperelasticityWeakForm<StrainEnergy>::evaluate_pressure_stress_residual
  ( libMesh::Number J,
    libMesh::Number press,
    const libMesh::Tensor & F_times_Cinv,
    const libMesh::RealGradient & dphi_i_times_JxW,
    libMesh::Number & Fu,
    libMesh::Number & Fv ) const
  {

    for( int alpha = 0; alpha < 2 /*dim*/; alpha++)
      {
        libMesh::Number c = press*J*dphi_i_times_JxW(alpha);

        Fu += F_times_Cinv(0,alpha)*c;
        Fv += F_times_Cinv(1,alpha)*c;
      }
  }

  template<typename StrainEnergy>
  inline
  void IncompressibleHyperelasticityWeakForm<StrainEnergy>::evaluate_pressure_stress_residual
  ( libMesh::Number J,
    libMesh::Number press,
    const libMesh::Tensor & F_times_Cinv,
    const libMesh::RealGradient & dphi_i_times_JxW,
    libMesh::Number & Fu,
    libMesh::Number & Fv,
    libMesh::Number & Fw ) const
  {
    for( int alpha = 0; alpha < 3 /*dim*/; alpha++)
      {
        libMesh::Number c = press*J*dphi_i_times_JxW(alpha);

        Fu += F_times_Cinv(0,alpha)*c;
        Fv += F_times_Cinv(1,alpha)*c;
        Fw += F_times_Cinv(2,alpha)*c;
      }
  }

  template<typename StrainEnergy>
  inline
  void IncompressibleHyperelasticityWeakForm<StrainEnergy>::evaluate_pressure_stress_pressure_jacobian
  ( libMesh::Number J,
    const libMesh::Tensor & F_times_Cinv,
    const libMesh::Real & p_phi_j,
    const libMesh::RealGradient & dphi_i_times_JxW,
    libMesh::Number & Kup,
    libMesh::Number & Kvp ) const
  {
    for( int alpha = 0; alpha < 2 /*dim*/; alpha++)
      {
        libMesh::Number c = p_phi_j*J*dphi_i_times_JxW(alpha);

        Kup += F_times_Cinv(0,alpha)*c;
        Kvp += F_times_Cinv(1,alpha)*c;
      }
  }


  template<typename StrainEnergy>
  inline
  void IncompressibleHyperelasticityWeakForm<StrainEnergy>::evaluate_pressure_stress_displacement_jacobian
  ( libMesh::Number J,
    libMesh::Number press,
    const libMesh::Tensor & F,
    const libMesh::Tensor & Cinv,
    const libMesh::Tensor & F_times_Cinv,
    const libMesh::RealGradient & dphi_i_times_JxW,
    const libMesh::RealGradient & dphi_j,
    libMesh::Number & Kuu,
    libMesh::Number & Kuv,
    libMesh::Number & Kvu,
    libMesh::Number & Kvv ) const
  {
    // residual is
    // p*J*F_times_Cinv*dphi_i*JxW
    // So we have three term: J, then F, then Cinv

    libMesh::Number pJ = press*J;

    // We'll need these in a couple of places
    libMesh::Gradient FCinv_i = F_times_Cinv*dphi_i_times_JxW;
    libMesh::Gradient FCinv_j = pJ*(F_times_Cinv*dphi_j);
    libMesh::Number dphi_j_Cinv_phi_i = pJ*dphi_j*(Cinv*dphi_i_times_JxW);

    // J term: p*dJ/du*F*Cinv*dphi_i*JxW
    // --> p*J*( (F*Cinv)*dphi_j * (F*Cinv)*dphi_i )*JxW
    {
      Kuu += FCinv_j(0)*FCinv_i(0);
      Kuv += FCinv_j(1)*FCinv_i(0);
      Kvu += FCinv_j(0)*FCinv_i(1);
      Kvv += FCinv_j(1)*FCinv_i(1);
    }

    // F term: p*J*dF/du*Cinv*dphi_i*JxW
    // --> p*J*dphi_j*Cinv*dphi_i*JxW
    {
      libMesh::Number term2 = dphi_j_Cinv_phi_i;
      Kuu += term2;
      Kvv += term2;
    }

    // Cinv term: p*J*F*dCinv/du*dphi_i*JxW
    {
      Kuu -= FCinv_j(0)*FCinv_i(0);
      Kuv -= FCinv_j(0)*FCinv_i(1);
      Kvu -= FCinv_j(1)*FCinv_i(0);
      Kvv -= FCinv_j(1)*FCinv_i(1);

      const int dim = 2;
      libMesh::Number term3 = dphi_j_Cinv_phi_i;
      for (int M = 0; M < dim; M++)
        for (int K=0; K < dim; K++)
          {
            Kuu -= term3*F(0,K)*Cinv(K,M)*F(0,M);
            Kuv -= term3*F(0,K)*Cinv(K,M)*F(1,M);
            Kvu -= term3*F(1,K)*Cinv(K,M)*F(0,M);
            Kvv -= term3*F(1,K)*Cinv(K,M)*F(1,M);
          }
    }
  }

  template<typename StrainEnergy>
  inline
  void IncompressibleHyperelasticityWeakForm<StrainEnergy>::evaluate_pressure_stress_displacement_jacobian
  ( libMesh::Number J,
    libMesh::Number press,
    const libMesh::Tensor & F,
    const libMesh::Tensor & Cinv,
    const libMesh::Tensor & F_times_Cinv,
    const libMesh::RealGradient & dphi_i_times_JxW,
    const libMesh::RealGradient & dphi_j,
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
    // residual is
    // p*J*F_times_Cinv*dphi_i*JxW
    // So we have three term: J, then F, then Cinv

    libMesh::Number pJ = press*J;

    // We'll need these in a couple of places
    libMesh::Gradient FCinv_i = F_times_Cinv*dphi_i_times_JxW;
    libMesh::Gradient FCinv_j = pJ*(F_times_Cinv*dphi_j);
    libMesh::Number dphi_j_Cinv_phi_i = pJ*dphi_j*(Cinv*dphi_i_times_JxW);

    // J term: p*dJ/du*F*Cinv*dphi_i*JxW
    // --> p*J*( (F*Cinv)*dphi_j * (F*Cinv)*dphi_i )*JxW
    {
      Kuu += FCinv_j(0)*FCinv_i(0);
      Kuv += FCinv_j(1)*FCinv_i(0);
      Kuw += FCinv_j(2)*FCinv_i(0);
      Kvu += FCinv_j(0)*FCinv_i(1);
      Kvv += FCinv_j(1)*FCinv_i(1);
      Kvw += FCinv_j(2)*FCinv_i(1);
      Kwu += FCinv_j(0)*FCinv_i(2);
      Kwv += FCinv_j(1)*FCinv_i(2);
      Kww += FCinv_j(2)*FCinv_i(2);
    }

    // F term: p*J*dF/du*Cinv*dphi_i*JxW
    // --> p*J*dphi_j*Cinv*dphi_i*JxW
    {
      libMesh::Number term2 = dphi_j_Cinv_phi_i;
      Kuu += term2;
      Kvv += term2;
      Kww += term2;
    }

    // Cinv term: p*J*F*dCinv/du*dphi_i*JxW
    {
      Kuu -= FCinv_j(0)*FCinv_i(0);
      Kuv -= FCinv_j(0)*FCinv_i(1);
      Kuw -= FCinv_j(0)*FCinv_i(2);
      Kvu -= FCinv_j(1)*FCinv_i(0);
      Kvv -= FCinv_j(1)*FCinv_i(1);
      Kvw -= FCinv_j(1)*FCinv_i(2);
      Kwu -= FCinv_j(2)*FCinv_i(0);
      Kwv -= FCinv_j(2)*FCinv_i(1);
      Kww -= FCinv_j(2)*FCinv_i(2);

      const int dim = 3;
      libMesh::Number term3 = dphi_j_Cinv_phi_i;
      for (int M = 0; M < dim; M++)
        for (int K=0; K < dim; K++)
          {
            Kuu -= term3*F(0,K)*Cinv(K,M)*F(0,M);
            Kuv -= term3*F(0,K)*Cinv(K,M)*F(1,M);
            Kuw -= term3*F(0,K)*Cinv(K,M)*F(2,M);
            Kvu -= term3*F(1,K)*Cinv(K,M)*F(0,M);
            Kvv -= term3*F(1,K)*Cinv(K,M)*F(1,M);
            Kvw -= term3*F(1,K)*Cinv(K,M)*F(2,M);
            Kwu -= term3*F(2,K)*Cinv(K,M)*F(0,M);
            Kwv -= term3*F(2,K)*Cinv(K,M)*F(1,M);
            Kww -= term3*F(2,K)*Cinv(K,M)*F(2,M);
          }
    }
  }

  template<typename StrainEnergy>
  inline
  void IncompressibleHyperelasticityWeakForm<StrainEnergy>::evaluate_pressure_stress_pressure_jacobian
  ( libMesh::Number J,
    const libMesh::Tensor & F_times_Cinv,
    const libMesh::Real & p_phi_j,
    const libMesh::RealGradient & dphi_i_times_JxW,
    libMesh::Number & Kup,
    libMesh::Number & Kvp,
    libMesh::Number & Kwp ) const
  {
    for( int alpha = 0; alpha < 3 /*dim*/; alpha++)
      {
        libMesh::Number c = p_phi_j*J*dphi_i_times_JxW(alpha);

        Kup += F_times_Cinv(0,alpha)*c;
        Kvp += F_times_Cinv(1,alpha)*c;
        Kwp += F_times_Cinv(2,alpha)*c;
      }
  }

  template<typename StrainEnergy>
  inline
  void IncompressibleHyperelasticityWeakForm<StrainEnergy>::evaluate_pressure_constraint_residual
  ( libMesh::Number J,
    const libMesh::Real & phi_times_JxW,
    libMesh::Number & Fp ) const
  {
    libMesh::Number Up = this->dU(J);
    Fp += Up*phi_times_JxW;
  }

  template<typename StrainEnergy>
  inline
  libMesh::Number IncompressibleHyperelasticityWeakForm<StrainEnergy>::U( libMesh::Number J ) const
  {
    libMesh::Number lnJ = std::log(J);
    return 0.5*lnJ*lnJ;
  }

  template<typename StrainEnergy>
  inline
  libMesh::Number IncompressibleHyperelasticityWeakForm<StrainEnergy>::dU( libMesh::Number J ) const
  {
    return std::log(J)/J;
  }

  template<typename StrainEnergy>
  inline
  libMesh::Number IncompressibleHyperelasticityWeakForm<StrainEnergy>::d2U( libMesh::Number J ) const
  {
    return (1.0-std::log(J))/(J*J);
  }

} // end namespace GRINS

#endif // GRINS_INCOMPRESSIBLE_HYPERELASTICITY_WEAK_FORM_H
