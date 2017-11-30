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

// GRINS
#include "grins/variable_warehouse.h"
#include "grins/physical_constants.h"
#include "grins/math_constants.h"

#if GRINS_HAVE_ANTIOCH
#include "grins/antioch_chemistry.h"
#endif

#if GRINS_HAVE_CANTERA
#include "grins/cantera_mixture.h"
#endif

// libMesh
#include "libmesh/fe.h"
#include "libmesh/elem.h"

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
      _rad_coeff(Constants::second_rad_const * 100) // [cm K]
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
    START_LOG("operator()","AbsorptionCoeff");
    libMesh::Real T,p,thermo_p; // temperature, hydrostatic pressure, thermodynamic pressure
    std::vector<libMesh::Real> Y(_chemistry->n_species()); // mass fractions
    
    for (unsigned int s=0; s<_chemistry->n_species(); s++)
      context.point_value(_Y_var.species(s), qp_xyz, Y[s]);

    context.point_value(_T_var.T(), qp_xyz, T); // [K]

    if (_calc_thermo_pressure) {
      libmesh_not_implemented();
    } else {
      thermo_p = _thermo_pressure;
    }

    context.point_value(_P_var.p(), qp_xyz, p); // [Pa]

    libMesh::Real P = p + thermo_p; // total pressure [Pa]
    libmesh_assert_greater(P,0.0);

    libMesh::Real kv = 0.0;

    for (unsigned int i=_min_index; i<=_max_index; i++)
      kv += this->kv(T,P,Y,i);

    STOP_LOG("operator()","AbsorptionCoeff");
    return kv;
  }

  template<typename Chemistry>
  void AbsorptionCoeff<Chemistry>::operator()( const libMesh::FEMContext & /*context*/,
                                               const libMesh::Point & /*p*/,
                                               const libMesh::Real /*time*/,
                                               libMesh::DenseVector<libMesh::Real> & /*output*/)
  {
    libmesh_not_implemented();
  }

  template<typename Chemistry>
  void AbsorptionCoeff<Chemistry>::derivatives( libMesh::FEMContext & /*context*/,
                                                const libMesh::Point & /*qp_xyz*/,
                                                const libMesh::Real & /*JxW*/,
                                                const unsigned int /*qoi_index*/,
                                                const libMesh::Real /*time*/)
  {
    // TODO
    libmesh_not_implemented();
  }

  template<typename Chemistry>
  libMesh::UniquePtr<libMesh::FEMFunctionBase<libMesh::Real> > AbsorptionCoeff<Chemistry>::clone() const
  {
    libmesh_not_implemented();
  }


  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::kv(libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> Y, unsigned int i)
  {
    // linestrength
    libMesh::Real S = this->Sw(T,P,i);

    // Voigt profile [cm^-1]
    libMesh::Real phi_V = this->voigt(T,P,Y,i);

    libMesh::Real M_mix = _chemistry->M_mix(Y);
    libMesh::Real X = _chemistry->X(_species_idx,M_mix,Y[_species_idx]); 

    // absorption coefficient [cm^-1]
    return S*(P/Constants::atmosphere_Pa)*X*phi_V;
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::d_kv_dT(libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> Y, unsigned int i)
  {
    libMesh::Real dS = dS_dT(T,P,i);
    libMesh::Real dV = d_voigt_dT(T,P,Y,i);

    libMesh::Real S = this->Sw(T,P,i);
    libMesh::Real V = this->voigt(T,P,Y,i);
    libMesh::Real M_mix = _chemistry->M_mix(Y);
    libMesh::Real X = _chemistry->X(_species_idx,M_mix,Y[_species_idx]);

    return X*(P/Constants::atmosphere_Pa) * ( S*dV + dS*V );
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::d_kv_dP(libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> Y, unsigned int i)
  {
    libMesh::Real dS = dS_dP(T,P,i);
    libMesh::Real dV = d_voigt_dP(T,P,Y,i);

    libMesh::Real S = this->Sw(T,P,i);
    libMesh::Real V = this->voigt(T,P,Y,i);
    libMesh::Real M_mix = _chemistry->M_mix(Y);
    libMesh::Real X = _chemistry->X(_species_idx,M_mix,Y[_species_idx]);

    return X * ( S*(P/Constants::atmosphere_Pa)*dV + S*(1.0/Constants::atmosphere_Pa)*V + dS*P*V );
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::d_kv_dY(libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> Y, unsigned int i)
  {
    libMesh::Real dX = dX_dY(Y);
    libMesh::Real dV = d_voigt_dY(T,P,Y,i);

    libMesh::Real S = this->Sw(T,P,i);
    libMesh::Real V = this->voigt(T,P,Y,i);
    libMesh::Real M_mix = _chemistry->M_mix(Y);
    libMesh::Real X = _chemistry->X(_species_idx,M_mix,Y[_species_idx]);
    
    return S*(P/Constants::atmosphere_Pa) * ( X*dV + dX*V );
  }


  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::Sw(libMesh::Real T, libMesh::Real P, unsigned int i)
  {
    // isotopologue
    unsigned int iso = _hitran->isotopologue(i);

    // linecenter wavenumber
    libMesh::Real nu = this->get_nu(P,i);

    // linestrength
    libMesh::Real sw0 = _hitran->sw(i);

    // lower state energy of transition
    libMesh::Real E = _hitran->elower(i);

    // partition function at reference temp
    libMesh::Real QT0 = _hitran->partition_function(_T0,iso);

    // partition function at current temp
    libMesh::Real QT = _hitran->partition_function(T,iso);

    libMesh::Real a = sw0;
    libMesh::Real b = (QT0/QT);
    libMesh::Real c = std::exp(-E*_rad_coeff*( (1.0/T) - (1.0/_T0) ));
    libMesh::Real d = ( 1.0-std::exp(-_rad_coeff*nu/T) );
    libMesh::Real e = std::pow(1.0-std::exp(-_rad_coeff*nu/_T0),-1.0);

    libMesh::Real S = a*b*c*d*e;

    // convert linestrength units to [cm^-2 atm^-1]
    libMesh::Real loschmidt = (Constants::atmosphere_Pa)/(T*Constants::Boltzmann*1.0e6);
    S = S*loschmidt;

    return S;
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::dS_dT(libMesh::Real T, libMesh::Real P, unsigned int i)
  {
    libMesh::Real iso = _hitran->isotopologue(i);
    libMesh::Real sw = _hitran->sw(i);
    libMesh::Real E = _hitran->elower(i);
    libMesh::Real QT0 = _hitran->partition_function(_T0,iso);
    libMesh::Real QT = _hitran->partition_function(T,iso);
    libMesh::Real dQT = this->dQ_dT(T,iso);
    libMesh::Real nu = this->get_nu(P,i);
    
    libMesh::Real constants = sw * QT0 * std::pow( (1.0-std::exp(-_rad_coeff*nu/_T0)), -1.0) * Constants::atmosphere_Pa/(Constants::Boltzmann*1.0e6);
    libMesh::Real A  = 1.0/T;
    libMesh::Real dA = -1.0/(T*T);

    libMesh::Real B  = 1.0/QT;
    libMesh::Real dB = -1.0/(QT*QT) * dQT;

    libMesh::Real C  = std::exp(-_rad_coeff*E* (1.0/T - 1.0/_T0) );
    libMesh::Real dC = std::exp(-_rad_coeff*E* (1.0/T - 1.0/_T0) ) * (_rad_coeff*E/(T*T));

    libMesh::Real D  = 1.0 - std::exp(-_rad_coeff*nu/T);
    libMesh::Real dD = -std::exp(-_rad_coeff*nu/T) * (_rad_coeff*nu/(T*T));

    return constants * ( A*B*C*dD + A*B*dC*D + A*dB*C*D + dA*B*C*D );
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::dS_dP(libMesh::Real T, libMesh::Real P, unsigned int i)
  {
    libMesh::Real iso = _hitran->isotopologue(i);
    libMesh::Real sw = _hitran->sw(i);
    libMesh::Real E = _hitran->elower(i);
    libMesh::Real QT0 = _hitran->partition_function(_T0,iso);
    libMesh::Real QT = _hitran->partition_function(T,iso);
    libMesh::Real nu = this->get_nu(P,i);
    libMesh::Real dnu = this->d_nu_dP(i);

    libMesh::Real constants = sw * QT0/QT * std::exp(-_rad_coeff*E*(1.0/T - 1.0/_T0)) * Constants::atmosphere_Pa/(Constants::Boltzmann*1.0e6 * T);

    libMesh::Real A  = 1.0 - std::exp(-_rad_coeff*nu/T);
    libMesh::Real dA = -std::exp(-_rad_coeff*nu/T) * (-_rad_coeff/T) * dnu;

    libMesh::Real B  = std::pow( (1.0-std::exp(-_rad_coeff*nu/T)), -1.0);
    libMesh::Real dB = -std::pow( (1.0-std::exp(-_rad_coeff*nu/T)), -2.0) * -std::exp(-_rad_coeff*nu/_T0) * (-_rad_coeff/_T0) * dnu;

    return constants * ( A*dB + dA*B );
  }


  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::nu_D(libMesh::Real T, libMesh::Real P, unsigned int i)
  {
    libMesh::Real k = Constants::Boltzmann;
    libMesh::Real c = Constants::c_vacuum;
    libMesh::Real NA = Constants::Avogadro;
    libMesh::Real M = _chemistry->M(_species_idx);
    libMesh::Real nu = this->get_nu(P,i);

    return (nu/c)*std::sqrt( ( 8.0*k*T*std::log(2.0) )/( M/NA ) );
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::d_nuD_dT(libMesh::Real T, libMesh::Real P, unsigned int i)
  {
    libMesh::Real k = Constants::Boltzmann;
    libMesh::Real c = Constants::c_vacuum;
    libMesh::Real NA = Constants::Avogadro;
    libMesh::Real M = _chemistry->M(_species_idx);
    libMesh::Real nu = this->get_nu(P,i);

    return (nu/c) * std::sqrt( (8.0*k*std::log(2.0))/(M/NA) ) * 0.5 * 1.0/(std::sqrt(T));
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::d_nuD_dP(libMesh::Real T, unsigned int i)
  {
    libMesh::Real k = Constants::Boltzmann;
    libMesh::Real c = Constants::c_vacuum;
    libMesh::Real NA = Constants::Avogadro;
    libMesh::Real M = _chemistry->M(_species_idx);
    libMesh::Real dnu = this->d_nu_dP(i);

    return (1.0/c)*std::sqrt( ( 8.0*k*T*std::log(2.0) )/( M/NA ) ) * dnu;
  }


  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::nu_C(libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> Y, unsigned int i)
  {
    libMesh::Real g_self = _hitran->gamma_self(i);
    libMesh::Real g_air = _hitran->gamma_air(i);
    libMesh::Real n = _hitran->n_air(i);

    libMesh::Real M_mix = _chemistry->M_mix(Y);
    libMesh::Real X = _chemistry->X(_species_idx,M_mix,Y[_species_idx]);

    return 2.0*(P/Constants::atmosphere_Pa)*std::pow(_T0/T,n)*( X*g_self + (1.0-X)*g_air );
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::d_nuC_dT(libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> Y, unsigned int i)
  {
    libMesh::Real g_self = _hitran->gamma_self(i);
    libMesh::Real g_air = _hitran->gamma_air(i);
    libMesh::Real n = _hitran->n_air(i);

    libMesh::Real M_mix = _chemistry->M_mix(Y);
    libMesh::Real X = _chemistry->X(_species_idx,M_mix,Y[_species_idx]);

    return 2.0*(P/Constants::atmosphere_Pa) * n*std::pow(_T0/T,n-1.0)*(-_T0/(T*T)) * ( X*g_self + (1.0-X)*g_air );
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::d_nuC_dP(libMesh::Real T, std::vector<libMesh::Real> Y, unsigned int i)
  {
    libMesh::Real g_self = _hitran->gamma_self(i);
    libMesh::Real g_air = _hitran->gamma_air(i);
    libMesh::Real n = _hitran->n_air(i);

    libMesh::Real M_mix = _chemistry->M_mix(Y);
    libMesh::Real X = _chemistry->X(_species_idx,M_mix,Y[_species_idx]);

    return 2.0*(1.0/Constants::atmosphere_Pa)*std::pow(_T0/T,n)*( X*g_self + (1.0-X)*g_air );
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::d_nuC_dY(libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> Y, unsigned int i)
  {
    libMesh::Real g_self = _hitran->gamma_self(i);
    libMesh::Real g_air = _hitran->gamma_air(i);
    libMesh::Real n = _hitran->n_air(i);
    libMesh::Real dX = dX_dY(Y);

    return 2.0*(P/Constants::atmosphere_Pa)*std::pow(_T0/T,n)*( dX*g_self - dX*g_air );
  }


  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::voigt(libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> Y, unsigned int i)
  {
    libMesh::Real nu_D = this->nu_D(T,P,i);

    libMesh::Real a = this->voigt_a(T,P,Y,i);
    libMesh::Real w = this->voigt_w(T,P,i);

    // Voigt coefficient
    libMesh::Real V = 0.0;

    for(int i=0; i<4; i++)
      {
        libMesh::Real Ai = _voigt_coeffs[0][i];
        libMesh::Real Bi = _voigt_coeffs[1][i];
        libMesh::Real Ci = _voigt_coeffs[2][i];
        libMesh::Real Di = _voigt_coeffs[3][i];

        libMesh::Real aAi = a-Ai;
        libMesh::Real wBi = w-Bi;
        libMesh::Real aAi2 = std::pow(a-Ai,2.0);
        libMesh::Real wBi2 = std::pow(w-Bi,2.0);

        V += ( Ci*aAi + Di*wBi )/( aAi2 + wBi2 );
      }

    libMesh::Real phi_V = (2.0*std::sqrt(std::log(2.0)))/(std::sqrt(Constants::pi)*nu_D)*V;

    return phi_V;
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::d_voigt_dT(libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> Y, unsigned int i)
  {
    libMesh::Real nu_D = this->nu_D(T,P,i);
    libMesh::Real dnu_D = d_nuD_dT(T,P,i);

    libMesh::Real a  = this->voigt_a(T,P,Y,i);
    libMesh::Real da = this->d_voigt_a_dT(T,P,Y,i);

    libMesh::Real w  = this->voigt_w(T,P,i);
    libMesh::Real dw = this->d_voigt_w_dT(T,P,i); 

    // Voigt coefficient
    libMesh::Real V  = 0.0;
    libMesh::Real dV = 0.0;

    for(int i=0; i<4; i++)
      {
        libMesh::Real Ai = _voigt_coeffs[0][i];
        libMesh::Real Bi = _voigt_coeffs[1][i];
        libMesh::Real Ci = _voigt_coeffs[2][i];
        libMesh::Real Di = _voigt_coeffs[3][i];

        libMesh::Real aAi = a-Ai;
        libMesh::Real wBi = w-Bi;
        libMesh::Real aAi2 = std::pow(a-Ai,2.0);
        libMesh::Real wBi2 = std::pow(w-Bi,2.0);

        V  += ( Ci*aAi + Di*wBi )/( aAi2 + wBi2 );

        dV += ( (aAi2+wBi2)*(Ci*da + Di*dw) - (Ci*aAi + Di*wBi)*(2.0*aAi*da + 2.0*wBi*dw) )/( std::pow(aAi2+wBi2,2.0) );
      }

    libMesh::Real constants = 2.0*std::sqrt(std::log(2.0))/std::sqrt(Constants::pi);
    libMesh::Real C  = 1.0/nu_D;
    libMesh::Real dC = -1.0/(nu_D*nu_D) * dnu_D;

    return constants * ( C*dV + dC*V );
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::d_voigt_dP(libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> Y, unsigned int i)
  {
    libMesh::Real nu_D = this->nu_D(T,P,i);
    libMesh::Real dnu_D = d_nuD_dP(T,i);

    libMesh::Real a  = this->voigt_a(T,P,Y,i);
    libMesh::Real da = this->d_voigt_a_dP(T,P,Y,i);

    libMesh::Real w  = this->voigt_w(T,P,i);
    libMesh::Real dw = this->d_voigt_w_dP(T,P,i);

    // Voigt coefficient
    libMesh::Real V  = 0.0;
    libMesh::Real dV = 0.0;

    for(int i=0; i<4; i++)
      {
        libMesh::Real Ai = _voigt_coeffs[0][i];
        libMesh::Real Bi = _voigt_coeffs[1][i];
        libMesh::Real Ci = _voigt_coeffs[2][i];
        libMesh::Real Di = _voigt_coeffs[3][i];

        libMesh::Real aAi = a-Ai;
        libMesh::Real wBi = w-Bi;
        libMesh::Real aAi2 = std::pow(a-Ai,2.0);
        libMesh::Real wBi2 = std::pow(w-Bi,2.0);

        V  += ( Ci*aAi + Di*wBi )/( aAi2 + wBi2 );

        dV += ( (aAi2+wBi2)*(Ci*da + Di*dw) - (Ci*aAi + Di*wBi)*(2.0*aAi*da + 2.0*wBi*dw) )/( std::pow(aAi2+wBi2,2.0) );
      }

    libMesh::Real constants = 2.0*std::sqrt(std::log(2.0))/std::sqrt(Constants::pi);
    libMesh::Real C  = 1.0/nu_D;
    libMesh::Real dC = -1.0/(nu_D*nu_D) * dnu_D;

    return constants * ( C*dV + dC*V );
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::d_voigt_dY(libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> Y, unsigned int i)
  {
    libMesh::Real nu_D = this->nu_D(T,P,i);

    libMesh::Real a  = this->voigt_a(T,P,Y,i);
    libMesh::Real da = this->d_voigt_a_dY(T,P,Y,i);

    libMesh::Real w  = this->voigt_w(T,P,i);

    // Voigt coefficient
    libMesh::Real dV = 0.0;

    for(int i=0; i<4; i++)
      {
        libMesh::Real Ai = _voigt_coeffs[0][i];
        libMesh::Real Bi = _voigt_coeffs[1][i];
        libMesh::Real Ci = _voigt_coeffs[2][i];
        libMesh::Real Di = _voigt_coeffs[3][i];

        libMesh::Real aAi = a-Ai;
        libMesh::Real wBi = w-Bi;
        libMesh::Real aAi2 = std::pow(a-Ai,2.0);
        libMesh::Real wBi2 = std::pow(w-Bi,2.0);

        dV += ( (aAi2+wBi2)*(Ci*da) - (Ci*aAi + Di*wBi)*(2.0*aAi*da) )/( std::pow(aAi2+wBi2,2.0) );
      }

    libMesh::Real constants = 2.0*std::sqrt(std::log(2.0))/(std::sqrt(Constants::pi)*nu_D);

    return constants * dV;
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

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::voigt_a(libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> Y, unsigned int i)
  {
    libMesh::Real nu_c = this->nu_C(T,P,Y,i);
    libMesh::Real nu_D = this->nu_D(T,P,i);

    return std::sqrt(std::log(2.0))*nu_c/nu_D;
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::d_voigt_a_dT(libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> Y, unsigned int i)
  {
    libMesh::Real nu_c = this->nu_C(T,P,Y,i);
    libMesh::Real nu_D = this->nu_D(T,P,i);
    libMesh::Real dnu_c = d_nuC_dT(T,P,Y,i);
    libMesh::Real dnu_D = d_nuD_dT(T,P,i);

    return std::sqrt(std::log(2.0)) * (nu_D*dnu_c - nu_c*dnu_D)/(nu_D*nu_D);
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::d_voigt_a_dP(libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> Y, unsigned int i)
  {
    libMesh::Real nu_c = this->nu_C(T,P,Y,i);
    libMesh::Real nu_D = this->nu_D(T,P,i);
    libMesh::Real dnu_c = d_nuC_dP(T,Y,i);
    libMesh::Real dnu_D = d_nuD_dP(T,i);

    return std::sqrt(std::log(2.0)) * (nu_D*dnu_c - nu_c*dnu_D)/(nu_D*nu_D);
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::d_voigt_a_dY(libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> Y, unsigned int i)
  {
    libMesh::Real nu_D = this->nu_D(T,P,i);
    libMesh::Real dnu_c = d_nuC_dY(T,P,Y,i);

    return std::sqrt(std::log(2.0))/nu_D * dnu_c;
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::voigt_w(libMesh::Real T, libMesh::Real P, unsigned int i)
  {
    libMesh::Real nu = this->get_nu(P,i);
    libMesh::Real nu_D = this->nu_D(T,P,i);

    return 2.0*std::sqrt(std::log(2.0))*(_nu-nu)/nu_D;
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::d_voigt_w_dT(libMesh::Real T, libMesh::Real P, unsigned int i)
  {
    libMesh::Real nu = this->get_nu(P,i);
    libMesh::Real nu_D = this->nu_D(T,P,i);
    libMesh::Real dnu_D = d_nuD_dT(T,P,i);

    return 2.0*std::sqrt(std::log(2.0))*(_nu-nu)/(nu_D*nu_D) * -dnu_D;
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::d_voigt_w_dP(libMesh::Real T, libMesh::Real P, unsigned int i)
  {
    libMesh::Real nu = this->get_nu(P,i);
    libMesh::Real dnu = d_nu_dP(i);
    libMesh::Real nu_D = this->nu_D(T,P,i);
    libMesh::Real dnu_D = d_nuD_dP(T,i);

    return 2.0*std::sqrt(std::log(2.0)) * ( (nu_D)*(-dnu) - (_nu-nu)*(dnu_D) )/(nu_D*nu_D);
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::get_nu(libMesh::Real P, unsigned int i)
  {
    libMesh::Real nu = _hitran->nu0(i);
    libMesh::Real d_air = _hitran->delta_air(i);

    return nu + d_air*((P/Constants::atmosphere_Pa)/_Pref);
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::d_nu_dP(unsigned int i)
  {
    libMesh::Real d_air = _hitran->delta_air(i);
    
    return d_air/_Pref * 1.0/Constants::atmosphere_Pa;
  }

  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::dX_dY(std::vector<libMesh::Real> Y)
  {
    libMesh::Real MW = _chemistry->M(_species_idx);
    libMesh::Real MW_mix = _chemistry->M_mix(Y);

    libMesh::Real Ys = Y[_species_idx];

    return 1.0/MW * ( MW_mix - Ys*(MW_mix*MW_mix)/MW );
  }


  template<typename Chemistry>
  libMesh::Real AbsorptionCoeff<Chemistry>::dQ_dT(libMesh::Real T, unsigned int iso)
  {
    libMesh::Real deriv = _hitran->partition_function_derivative(T,iso);
    return deriv;
  }

#if GRINS_HAVE_ANTIOCH
  template class AbsorptionCoeff<AntiochChemistry>;
#endif

#if GRINS_HAVE_CANTERA
  template class AbsorptionCoeff<CanteraMixture>;
#endif
} //namespace GRINS
