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
#include "grins/absorption_coeff.h"
#include "grins/variable_warehouse.h"
#include "grins/physical_constants.h"
#include "grins/math_constants.h"

#if GRINS_HAVE_ANTIOCH
#include "grins/antioch_chemistry.h"
#endif

#if GRINS_HAVE_CANTERA
#include "grins/cantera_mixture.h"
#endif

namespace GRINS
{
  template<typename Chemistry>
  AbsorptionCoeff<Chemistry>::AbsorptionCoeff(SharedPtr<Chemistry> & chem, SharedPtr<HITRAN> & hitran,
                                              libMesh::Real nu_min, libMesh::Real nu_max,
                                              libMesh::Real desired_nu, const std::string & species,
                                              libMesh::Real thermo_pressure)
    : _chemistry(chem),
      _hitran(hitran),
      _nu(desired_nu),
      _T_var(GRINSPrivate::VariableWarehouse::get_variable_subclass<PrimitiveTempFEVariables>("Temperature")),
      _P_var(GRINSPrivate::VariableWarehouse::get_variable_subclass<PressureFEVariable>("Pressure")),
      _Y_var(GRINSPrivate::VariableWarehouse::get_variable_subclass<SpeciesMassFractionsVariable>("SpeciesMassFractions")),
      _T0(296), // [K]
      _Pref(1), // [atm]
      _rad_coeff(1.43887752) // [cm K]
  {
    // sanity checks
    if ( (nu_min>nu_max) || (desired_nu>nu_max) || (desired_nu<nu_min) )
      {
        std::stringstream ss;
        ss <<"Invalid specification of wavenumber range:" <<std::endl
           <<"nu_min: " <<nu_min <<std::endl
           <<"nu_max: " <<nu_max <<std::endl
           <<"desired_nu: " <<desired_nu <<std::endl;   
        libmesh_error_msg(ss.str());
      }
    
    _species_idx = _chemistry->species_index(species);
    unsigned int data_size = _hitran->get_data_size();
    
    bool min_flag = false;
    
    for (unsigned int i=0; i<data_size; i++)
      {
        if (_hitran->nu0(i) > nu_min)
          {
            _min_index = i;
            min_flag = true;
            break;
          }
      }
      
    if (!min_flag)
      {
        std::stringstream ss;
        ss <<"Minimum wavenumber " <<nu_min <<" is greater than the maximum wavenumber in provided HITRAN data";
        libmesh_error_msg(ss.str());
      }
    
    bool max_flag = false;
    
    for (unsigned int i=data_size-1; i>=0; i--)
      {
        if (_hitran->nu0(i) < nu_max)
          {
            _max_index = i;
            max_flag = true;
            break;
          }
      }
      
    if (!max_flag)
      _max_index = data_size-1;

    if (thermo_pressure == -1.0) {
      _calc_thermo_pressure = true;
      _thermo_pressure = -1.0; // not needed in this case
    } else {
      _calc_thermo_pressure = false;
      _thermo_pressure = thermo_pressure;
    }

    this->init_voigt();
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::operator()(const libMesh::FEMContext& context, const libMesh::Point& qp_xyz, const libMesh::Real /*t*/)
  {
    libMesh::Real T,p,thermo_p; // temperature, hydrostatic pressure, thermodynamic pressure
    std::vector<libMesh::Real> Y(_chemistry->n_species()); // mass fractions

    context.point_value(_T_var.T(), qp_xyz, T); // [K]

    if (_calc_thermo_pressure) {
      libmesh_not_implemented();
    } else {
      thermo_p = _thermo_pressure;
    }

    context.point_value(_P_var.p(), qp_xyz, p); // [Pa]

    libMesh::Real P = p + thermo_p; // total pressure [Pa]
    libmesh_assert_greater(P,0.0);

    // all mass fractions needed to get M_mix
    for (unsigned int s=0; s<_chemistry->n_species(); s++)
      context.point_value(_Y_var.species(s), qp_xyz, Y[s]);

    libMesh::Real M = _chemistry->M(_species_idx); // [kg/mol]
    libMesh::Real M_mix = _chemistry->M_mix(Y); // [kg/mol]
    libMesh::Real X = _chemistry->X(_species_idx,M_mix,Y[_species_idx]);

    P /= 101325.0; // convert to [atm]

    return this->kv(P,T,X,M);

  }

  template<typename Chemistry>
  void AbsorptionCoeff<Chemistry>::operator()( const libMesh::FEMContext& /*context*/,
                                          const libMesh::Point& /*p*/,
                                          const libMesh::Real /*time*/,
                                          libMesh::DenseVector<libMesh::Real>& /*output*/)
  {
    libmesh_not_implemented();
  }

  template<typename Chemistry>
  libMesh::UniquePtr<libMesh::FEMFunctionBase<libMesh::Real> > AbsorptionCoeff<Chemistry>::clone() const
  {
    libmesh_not_implemented();
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::kv(libMesh::Real P,libMesh::Real T,libMesh::Real X,libMesh::Real M)
  {
    libMesh::Real kv = 0.0;
    
    for (unsigned int i=_min_index; i<=_max_index; i++)
      {
        // isotopologue 
        unsigned int iso = _hitran->isotopologue(i);
        
        // linecenter wavenumber
        libMesh::Real nu0 = _hitran->nu0(i);

        // linestrength
        libMesh::Real sw0 = _hitran->sw(i);

        // partition function at reference temp
        libMesh::Real QT0 = _hitran->partition_function(_T0,iso);

        // partition function at current temp
        libMesh::Real QT = _hitran->partition_function(T,iso);

        // lower state energy of transition
        libMesh::Real E = _hitran->elower(i);

        // air pressure-induced line shift
        libMesh::Real d_air = _hitran->delta_air(i);

        // pressure shift of the linecenter wavenumber
        libMesh::Real nu = nu0 + d_air*(P/_Pref);

        // linestrength 
        libMesh::Real S = sw0 * (QT0/QT) * std::exp(-E*_rad_coeff*( (1.0/T) - (1.0/_T0) )) * ( 1.0-std::exp(-_rad_coeff*nu/T) ) * pow(1.0-std::exp(-_rad_coeff*nu/_T0),-1.0);

        // convert linestrength units to [cm^-2 atm^-1]
        libMesh::Real loschmidt = (101325.0*1.0e-6)/(T*Constants::Boltzmann);
        S = S*loschmidt;

        // collisional FWHM [cm^-1]
        libMesh::Real nu_c = this->nu_C(T,X,P,i);

        // Doppler FWHM [cm^-1]
        libMesh::Real nu_D = this->nu_D(nu,T,M);

        // Voigt profile [cm^-1]
        libMesh::Real phi_V = this->voigt(nu_D,nu_c,nu);

        // absorption coefficient [cm^-1]
        kv += S*P*X*phi_V;
      }
      
    return kv;
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::nu_D(libMesh::Real nu,libMesh::Real T,libMesh::Real M)
  {
    libMesh::Real k = Constants::Boltzmann;
    libMesh::Real c = Constants::c_vacuum;
    libMesh::Real NA = Constants::Avogadro;

    return (nu/c)*std::sqrt( ( 8.0*k*T*std::log(2.0) )/( M/NA ) );
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::nu_C(libMesh::Real T,libMesh::Real X,libMesh::Real P, unsigned int index)
  {
    libMesh::Real g_self = _hitran->gamma_self(index);
    libMesh::Real g_air = _hitran->gamma_air(index);
    libMesh::Real n = _hitran->n_air(index);

    return 2.0*P*pow(_T0/T,n) * ( X*g_self + (1.0-X)*g_air );
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::voigt(libMesh::Real nu_D, libMesh::Real nu_c, libMesh::Real nu)
  {
    //repeated term
    libMesh::Real root_ln2 = std::sqrt(std::log(2.0));

    libMesh::Real a = root_ln2*nu_c/nu_D;
    libMesh::Real w = 2*root_ln2*(_nu-nu)/nu_D;

    // Voigt coefficient
    libMesh::Real V = 0.0;

    for(int i=0; i<4; i++)
     {
        libMesh::Real Ai = _voigt_coeffs[0][i];
        libMesh::Real Bi = _voigt_coeffs[1][i];
        libMesh::Real Ci = _voigt_coeffs[2][i];
        libMesh::Real Di = _voigt_coeffs[3][i];

        V += ( Ci*(a-Ai) + Di*(w-Bi) )/( (a-Ai)*(a-Ai) + (w-Bi)*(w-Bi) );
     }

    libMesh::Real phi_V = (2.0*root_ln2)/(std::sqrt(Constants::pi)*nu_D)*V;

    return phi_V;
  }

  template<typename Chemistry>
  void AbsorptionCoeff<Chemistry>::init_voigt() {
    _voigt_coeffs.resize(4);
    for (int i=0; i<4; i++)
        _voigt_coeffs[i].resize(4);

    _voigt_coeffs[0][0] = -1.2150;
    _voigt_coeffs[0][1] = -1.3509;
    _voigt_coeffs[0][2] = -1.2150;
    _voigt_coeffs[0][3] = -1.3509;
    _voigt_coeffs[1][0] =  1.2359;
    _voigt_coeffs[1][1] =  0.3786;
    _voigt_coeffs[1][2] = -1.2359;
    _voigt_coeffs[1][3] = -0.3786;
    _voigt_coeffs[2][0] = -0.3085;
    _voigt_coeffs[2][1] =  0.5906;
    _voigt_coeffs[2][2] = -0.3085;
    _voigt_coeffs[2][3] =  0.5906;
    _voigt_coeffs[3][0] =  0.0210;
    _voigt_coeffs[3][1] = -1.1858;
    _voigt_coeffs[3][2] = -0.0210;
    _voigt_coeffs[3][3] =  1.1858;
  }

#if GRINS_HAVE_ANTIOCH
template class AbsorptionCoeff<AntiochChemistry>;
#endif

#if GRINS_HAVE_CANTERA
template class AbsorptionCoeff<CanteraMixture>;
#endif
} //namespace GRINS
