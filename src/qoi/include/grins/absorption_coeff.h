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


#ifndef GRINS_ABSORPTION_COEFF_H
#define GRINS_ABSORPTION_COEFF_H

// libMesh
#include "libmesh/fem_function_base.h"
#include "libmesh/auto_ptr.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/hitran.h"
#include "grins/single_variable.h"
#include "grins/multicomponent_variable.h"
#include "grins/variable_warehouse.h"
#include "grins/fem_function_and_derivative_base.h"

namespace GRINS
{
  /*!
    Evaluates the Beer-Lambert Law at a given point in space. It is intended to be used with the IntegratedFunction class for QoI evaluation.

    We use the differential form for calculating the absorbance \f$ \alpha_{\nu} \f$

    \f$ \alpha_{\nu} = -\frac{dI_{\nu}}{I_{\nu}} = k_{\nu} dx \f$

    This class calculates the <i>spectral absorption coefficient</i>, denoted above as \f$ k_{\nu} \f$, which is passed back to IntegratedFunction::element_qoi()
    to evaluate the integral Beer-Lambert law

    \f$ \frac{I_{\nu}}{I_{\nu}^0} = \exp\left\{- \int_0^L k_{\nu} dx\right\} \f$

    This class operates internally in [cm], [K], and [atm] since those are the base units used in the HITRAN data.
    In addition, the coefficients used in calculating the Voigt profile require these units to be used.

    However, values given to the physics class(es) must be in standard SI units [m] and [Pa].

    A chemistry library (Antioch or Cantera) is also required.
  */
  template<typename Chemistry>
  class AbsorptionCoeff : public FEMFunctionAndDerivativeBase<libMesh::Real>
  {
  public:

    /*!
      @param chem AntiochChemistry or CanteraMixture object
      @param hitran A HITRAN object
      @param nu_min The minimum wavenumber to use in calculating \f$ k_{\nu} \f$, inclusive
      @param nu_max The maximum wavenumber to use in calculating \f$ k_{\nu} \f$, inclusive
      @param desired_nu Wavenumber at which to calculate the absorption, [\f$ cm^{-1} \f$]
      @param species The string representing the species of interest (much match species given in input file)
      @param termo_pressure The thermodynamic pressure (in [Pa]), or -1.0 if non-constant
    */
    AbsorptionCoeff(SharedPtr<Chemistry> & chem, SharedPtr<HITRAN> & hitran,
                    libMesh::Real nu_min, libMesh::Real nu_max,
                    libMesh::Real desired_nu, const std::string & species,
                    libMesh::Real thermo_pressure);

    //! Calculate the absorption coefficient at a quadratue point
    virtual libMesh::Real operator()(const libMesh::FEMContext & context, const libMesh::Point & qp_xyz, const libMesh::Real t);

    //! Not used
    virtual void operator()( const libMesh::FEMContext & context,
                             const libMesh::Point & p,
                             const libMesh::Real time,
                             libMesh::DenseVector<libMesh::Real> & output);

    //! Calculate the derivative of the absorption coefficient at a QP with respect to all state variables
    virtual void derivatives( libMesh::FEMContext & context,
                              const libMesh::Point & qp_xyz,
                              const libMesh::Real & JxW,
                              const unsigned int qoi_index,
                              const libMesh::Real time);

    //! Not used
    virtual libMesh::UniquePtr<libMesh::FEMFunctionBase<libMesh::Real> > clone() const;

  protected:
    //! Antioch/Cantera object
    SharedPtr<Chemistry> _chemistry;

    //! HITRAN
    SharedPtr<HITRAN> _hitran;

    //! Desired wavenumber [cm^-1]
    libMesh::Real _nu;

    PrimitiveTempFEVariables & _T_var;
    PressureFEVariable & _P_var;
    SpeciesMassFractionsVariable & _Y_var;

    //! Reference temperature [K]
    libMesh::Real _T0;

    //! Reference pressure [atm]
    libMesh::Real _Pref;

    //! Second radiation coefficient [cm K]
    libMesh::Real _rad_coeff;

    //! Index of minimum wavenumber
    unsigned int _min_index;

    //! Index of maximum wavenumber
    unsigned int _max_index;

    //! Thermodynamic Pressure [atm]
    libMesh::Real _thermo_pressure;

    //! Flag for whether Thermodynamic Pressure is calculated or constant
    bool _calc_thermo_pressure;

    //! Index for the species of interest
    unsigned int _species_idx;

    //! 2D coefficient matrix for approximating the Voigt profile
    std::vector<std::vector<libMesh::Real> > _voigt_coeffs;

    //! Absorption coefficient [cm^-1]
    libMesh::Real kv(libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> Y, unsigned int i);

    //! Absorption coefficient temperature derivative
    libMesh::Real d_kv_dT(libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> Y, unsigned int i);

    //! Absorption coefficient pressure derivative
    libMesh::Real d_kv_dP(libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> Y, unsigned int i);

    //! Absorption coefficient mass fraction derivative
    libMesh::Real d_kv_dY(libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> Y, unsigned int i);

    //! Linestrength [cm^-2 atm^-1]
    libMesh::Real Sw(libMesh::Real T, libMesh::Real P, unsigned int i);

    //! Linestrength temperature derivative
    libMesh::Real dS_dT(libMesh::Real T, libMesh::Real P, unsigned int i);

    //! Linestrength pressure derivative
    libMesh::Real dS_dP(libMesh::Real T, libMesh::Real P, unsigned int i);

    //! Doppler broadening [cm^-1]
    libMesh::Real nu_D(libMesh::Real T, libMesh::Real P, unsigned int i);

    //! Doppler broadening temperature derivative
    libMesh::Real d_nuD_dT(libMesh::Real T, libMesh::Real P, unsigned int i);

    //! Doppler broadening pressure derivative
    libMesh::Real d_nuD_dP(libMesh::Real T, unsigned int i);

    //! Collisional broadening [cm^-1]
    libMesh::Real nu_C(libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> Y, unsigned int i);

    //! Collisional broadening temperature derivative
    libMesh::Real d_nuC_dT(libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> Y, unsigned int i);

    //! Collisional broadening pressure derivative
    libMesh::Real d_nuC_dP(libMesh::Real T, std::vector<libMesh::Real> Y, unsigned int i);

    //! Collisional broadening mass fraction derivative
    libMesh::Real d_nuC_dY(libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> Y, unsigned int i);

    //! Calculate the Voigt profile [cm^-1]
    /*!
      See reference:

      Implementation of an efficient analytical approximation to the Voigt function for photoemission lineshape analysis\n
      McLean A, Mitchell C, Swanston D\n
      Journal of Electron Spectroscopy and Related Phenomena 1994 vol: 69 (2) pp: 125-132
    */
    libMesh::Real voigt(libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> Y, unsigned int i);

    //! Voigt profile temperature derivative
    libMesh::Real d_voigt_dT(libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> Y, unsigned int i);

    //! Voigt profile pressure derivative
    libMesh::Real d_voigt_dP(libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> Y, unsigned int i);

    //! Voigt profile mass fraction derivative
    libMesh::Real d_voigt_dY(libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> Y, unsigned int i);

    //! Initialize the coeff matrix for calculating the Voigt profile
    void init_voigt();

    //! Voigt a parameter
    libMesh::Real voigt_a(libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> Y, unsigned int i);

    //! Voigt a parameter temperature derivative
    libMesh::Real d_voigt_a_dT(libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> Y, unsigned int i);

    //! Voigt a parameter pressure derivative
    libMesh::Real d_voigt_a_dP(libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> Y, unsigned int i);

    //! Voigt a parameter mass fraction derivative
    libMesh::Real d_voigt_a_dY(libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> Y, unsigned int i);

    //! Voigt w parameter
    libMesh::Real voigt_w(libMesh::Real T, libMesh::Real P, unsigned int i);

    //! Voigt w parameter temperature derivative
    libMesh::Real d_voigt_w_dT(libMesh::Real T, libMesh::Real P, unsigned int i);

    //! Voigt w parameter pressure derivative
    libMesh::Real d_voigt_w_dP(libMesh::Real T, libMesh::Real P, unsigned int i);

    //! Pressure shift of linecenter wavenumber
    libMesh::Real get_nu(libMesh::Real P, unsigned int i);

    //! Derivative of pressure-shifted linecenter wavenumber
    libMesh::Real d_nu_dP(unsigned int i);

    //! Mole fraction derivative with respect to mass fraction
    libMesh::Real dX_dY(std::vector<libMesh::Real> Y);

    //! Partition Function derivative (finite difference)
    libMesh::Real dQ_dT(libMesh::Real T, unsigned int iso);

    //! User should not call empty constructor
    AbsorptionCoeff();
  };

}
#endif //GRINS_ABSORPTION_COEFF_H
