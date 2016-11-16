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

// GRINS
#include "grins/assembly_context.h"
#include "grins/hitran.h"
#include "grins/single_variable.h"
#include "grins/multicomponent_variable.h"
#include "grins/variable_warehouse.h"

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
  class AbsorptionCoeff : public libMesh::FEMFunctionBase<libMesh::Real>
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

    //! Not used
    virtual libMesh::UniquePtr<libMesh::FEMFunctionBase<libMesh::Real> > clone() const;

  private:
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
    libMesh::Real kv(libMesh::Real P,libMesh::Real T, libMesh::Real X, libMesh::Real M);

    //! Doppler broadening [cm^-1]
    libMesh::Real nu_D(libMesh::Real nu, libMesh::Real T,libMesh::Real M);

    //! Collisional broadening [cm^-1]
    libMesh::Real nu_C(libMesh::Real T,libMesh::Real X,libMesh::Real P,unsigned int index);

    //! Calculate the Voigt profile [cm^-1]
    /*!
      See reference:

      Implementation of an efficient analytical approximation to the Voigt function for photoemission lineshape analysis\n
      McLean A, Mitchell C, Swanston D\n
      Journal of Electron Spectroscopy and Related Phenomena 1994 vol: 69 (2) pp: 125-132
    */
    libMesh::Real voigt(libMesh::Real nu_D, libMesh::Real nu_c, libMesh::Real nu);

    //! Initialize the coeff matrix for calculating the Voigt profile
    void init_voigt();
    
    //! User should not call empty constructor
    AbsorptionCoeff();
  };

}
#endif //GRINS_ABSORPTION_COEFF_H
