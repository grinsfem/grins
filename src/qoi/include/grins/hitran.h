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

#ifndef GRINS_HITRAN_H
#define GRINS_HITRAN_H


// C++
#include <vector>
#include <string>

// libMesh
#include "libmesh/libmesh.h"

namespace GRINS
{

  //! HITRAN interface object
  /*!
    This class provides an interface to utilize data from the HITRAN Spectroscopic Database.
    Data can be pulled directly from the HITRAN website, or through the HITRAN python API, HAPI.

    Two files must be provided to this class, one with various data members, and one with values
    of the partition function at various temperatures.

    The data file must contain comma-separated values for individual spectrscopic lines in the following format:\n
    Isotopologue,Linecenter,Linestrength,Air-broadened Half Width,Self-Broadened Half Width,Lower State Energy,Temperature Coefficient,Air-Pressure Shift
    
    The partition function file must contain comma-separated values for each isotopologue at the same temperature values.
    The temperature limits themselves are specified as inputs to the constructor.

    Note the partition function data <b>must</b> explicitly include a value at the given reference temperature,
    which is currently 296K for HITRAN.

    The basic units for values are those given by HITRAN, which currently are [cm] and [atm]
  */
  class HITRAN
  {
  public:
    //! Constructor
    /*!
      @param data_file A file containing values for the required quantities, formatted as described above
      @param partition_function_file A file containing partition function values, formatted as described above
      @param T_min The minimum temperature for the partition function values, inclusive
      @param T_max The maximum temperature for the partition function values, inclusive
      @param T_step The step between successive temperature values
    */
    HITRAN(const std::string & data_file, const std::string & partition_function_file,
           libMesh::Real T_min, libMesh::Real T_max, libMesh::Real T_step);

    //! Isotopologue ID
    unsigned int isotopologue(unsigned int index);    

    //! Linecenter wavenumber [\f$ cm^{-1} \f$]
    libMesh::Real nu0(unsigned int index);

    //! Linestrength [\f$ \frac{cm^{-1}}{molecule \cdot cm^{-2}} \f$]
    libMesh::Real sw(unsigned int index);

    //! Air-broadening half-width [\f$ cm^{-1} atm^{-1} \f$]
    libMesh::Real gamma_air(unsigned int index);

    //! Self-broadening half-wdith [\f$ cm^{-1} atm^{-1} \f$]
    libMesh::Real gamma_self(unsigned int index);

    //! Lower state energy of transition [\f$ cm^{-1} \f$]
    libMesh::Real elower(unsigned int index);

    //! Temperature coefficient []
    libMesh::Real n_air(unsigned int index);

    //! Air pressure-induced line shift [\f$ cm^{-1} atm^{-1} \f$]
    libMesh::Real delta_air(unsigned int index);

    //! Returns the value of the partition function
    //! of the given isotopologue at the given temperature
    libMesh::Real partition_function(libMesh::Real T, unsigned int iso);

    //! Return the data size
    unsigned int get_data_size();

  protected:
    libMesh::Real _Tmin, _Tmax, _Tstep;
    
    //! Size of spectroscopic data
    int _data_size;

    //! Size of partition function data
    int _q_size;

    //! Reference temperature (296K)
    libMesh::Real _T0;

    //! Value of partition function at reference temperature for all isotopologues.
    //! Cached since it is used frequently
    std::vector<libMesh::Real> _qT0;

    // Vectors of data values
    std::vector<unsigned int> _isotop;
    std::vector<libMesh::Real> _nu;
    std::vector<libMesh::Real> _sw;
    std::vector<libMesh::Real> _gamma_air;
    std::vector<libMesh::Real> _gamma_self;
    std::vector<libMesh::Real> _elower;
    std::vector<libMesh::Real> _n;
    std::vector<libMesh::Real> _delta_air;

    //! Vector for partition function values for all isotopologues
    std::vector<std::vector<libMesh::Real>> _qT;

    //! Find the index into _T corresponding to the given temperature
    int _T_index(libMesh::Real T);

    //! Search through the partition function data to get value for given temperature
    libMesh::Real search_partition_function(libMesh::Real T, unsigned int iso);

    //! Linear interpolation helper function
    libMesh::Real interpolate_values( int index_r, libMesh::Real T_star, const std::vector<libMesh::Real> & y) const;

    //! User should not call empty constructor
    HITRAN();
  };

}

#endif // GRINS_HITRAN_H
