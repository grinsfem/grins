
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

// this class
#include "grins/hitran.h"

// libMesh

// GRINS
#include "grins/string_utils.h"

// C++
#include <fstream>
#include <algorithm>

namespace GRINS
{
  HITRAN::HITRAN(const std::string & data_file, const std::string & partition_function_file,
                 libMesh::Real T_min, libMesh::Real T_max, libMesh::Real T_step)
  : _Tmin(T_min),
    _Tmax(T_max),
    _Tstep(T_step),
    _T0(296.0)
  {
    // sanity checks on temperature range specification
    if ( (T_min<0.0) || (T_min>=T_max) || (T_step<=0.0) || (T_min>_T0) || (T_max<_T0) )
      {
        std::stringstream ss;
        ss <<"Invalid specification of temperature range:" <<std::endl;
        ss <<"T_min: " <<T_min <<std::endl;
        ss <<"T_max: " <<T_max <<std::endl;
        ss <<"T_step: " <<T_step <<std::endl;       
        libmesh_error_msg(ss.str());
      }
  
    // open data file
    std::ifstream hitran_file;
    hitran_file.open(data_file);
    if (!hitran_file.is_open())
      {
        std::stringstream ss;
        ss <<"Unable to open provided hitran_file: " <<data_file <<std::endl;    
        libmesh_error_msg(ss.str());
      }
    
    while(!hitran_file.eof())
      {
        std::string line;
        getline(hitran_file,line);
        
        if (line == "")
          continue;
        
        std::vector<libMesh::Real> vals;
        
        GRINS::StringUtilities::split_string_real(line,",",vals);
        
        libmesh_assert_equal_to(vals.size(),8);
        
        // HITRAN gives isotopologue numbers starting at 1
        // shift to zero-based to make indexing easier
        _isotop.push_back(static_cast<unsigned int>(vals[0]-1));
        
        _nu.push_back(vals[1]);
        _sw.push_back(vals[2]);
        _gamma_air.push_back(vals[3]);
        _gamma_self.push_back(vals[4]);
        _elower.push_back(vals[5]);
        _n.push_back(vals[6]);
        _delta_air.push_back(vals[7]);      
      }

    // sanity checks
    libmesh_assert_equal_to( _isotop.size(), _nu.size() );
    libmesh_assert_equal_to( _nu.size(), _sw.size() );
    libmesh_assert_equal_to( _sw.size(), _gamma_air.size() );
    libmesh_assert_equal_to( _gamma_air.size(), _gamma_self.size() );
    libmesh_assert_equal_to( _gamma_self.size(), _elower.size() );
    libmesh_assert_equal_to( _elower.size(), _n.size() );
    libmesh_assert_equal_to( _n.size(), _delta_air.size() );

    // save data length and close HITRAN data file
    _data_size = _isotop.size(); // all data vectors are the same length
    hitran_file.close();

    // file with partition function values
    std::ifstream qT_file;
    qT_file.open(partition_function_file);
    if (!qT_file.is_open())
      {
        std::stringstream ss;
        ss <<"Unable to open provided partition_function_file: " <<partition_function_file <<std::endl;     
        libmesh_error_msg(ss.str());
      }
    
    // number of temperature values
    unsigned int num_T = (T_max-T_min)/T_step + 1;
    
    // read the partition function values
    unsigned int counter = 0;
    
    while(!qT_file.eof())
      {
        std::string line;
        getline(qT_file,line);
        
        if (line == "")
          continue;
        
        _qT.push_back(std::vector<libMesh::Real>());
        
        GRINS::StringUtilities::split_string_real(line,",",_qT[counter]);
        
        // we should have a partition function value for each temperature
        libmesh_assert_equal_to(num_T,_qT[counter].size());
        
        counter++;   
      }

    // save length and close partition sum file
    _q_size = num_T;
    qT_file.close();
    
    // cache the partition function values at the referece temperature
    for(unsigned int i=0; i<_qT.size(); i++)
      _qT0.push_back(this->search_partition_function(_T0,i));

    return;
  }

  unsigned int HITRAN::get_data_size()
  {
    return _data_size;
  }

  unsigned int HITRAN::isotopologue(unsigned int index)
  {
    libmesh_assert_less(index,_isotop.size());
    return _isotop[index];
  }

  libMesh::Real HITRAN::nu0(unsigned int index)
  {
    libmesh_assert_less(index,_nu.size());
    return _nu[index];
  }

  libMesh::Real HITRAN::sw(unsigned int index)
  {
    libmesh_assert_less(index,_sw.size());
    return _sw[index];
  }

  libMesh::Real HITRAN::gamma_air(unsigned int index)
  {
    libmesh_assert_less(index,_gamma_air.size());
    return _gamma_air[index];
  }

  libMesh::Real HITRAN::gamma_self(unsigned int index)
  {
    libmesh_assert_less(index,_gamma_self.size());
    return _gamma_self[index];
  }

  libMesh::Real HITRAN::elower(unsigned int index)
  {
    libmesh_assert_less(index,_elower.size());
    return _elower[index];
  }

  libMesh::Real HITRAN::n_air(unsigned int index)
  {
    libmesh_assert_less(index,_n.size());
    return _n[index];
  }

  libMesh::Real HITRAN::delta_air(unsigned int index)
  {
    libmesh_assert_less(index,_delta_air.size());
    return _delta_air[index];
  }

  libMesh::Real HITRAN::partition_function(libMesh::Real T, unsigned int iso)
  {
    libMesh::Real qt;

    if (T==_T0)
      qt = _qT0[iso];
    else
      qt = this->search_partition_function(T,iso);

    return qt;
  }

  libMesh::Real HITRAN::search_partition_function(libMesh::Real T, unsigned int iso)
  {
    libMesh::Real retval = -1.0;

    int i = _T_index(T);

    if (i >= 0)
        retval = this->interpolate_values(i,T,_qT[iso]);
    else
      {
        std::stringstream ss;
        ss <<"Error: Temperature " <<T <<"K does not exist in the given partition sum data" <<std::endl;
        libmesh_error_msg(ss.str());
      }

    return retval;
  }

  int HITRAN::_T_index(libMesh::Real T)
  {
    unsigned int index = std::ceil((T-_Tmin)/_Tstep) + 1;
    return index;
  }
  
  libMesh::Real HITRAN::interpolate_values( int index_r, libMesh::Real T_star, const std::vector<libMesh::Real> & y) const
  {
    if ( (T_star>_Tmax) || (T_star<_Tmin) )
      {
        std::stringstream ss;
        ss <<"Error, temperature T=" <<T_star <<" is outside the specified range of provided partition function values" <<std::endl;
        libmesh_error_msg(ss.str());
      }
      
    libMesh::Real T = _Tmin + (_Tstep*index_r);
    return y[index_r-1] + ( (T_star-(T-_Tstep))*(y[index_r]-y[index_r-1]) )/(T-(T-_Tstep));
  }

}
