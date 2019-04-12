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


#ifndef GRINS_GAS_MIXTURE_H
#define GRINS_GAS_MIXTURE_H

// libMesh forward declarations
class GetPot;

namespace GRINS
{
  template<typename Chemistry, typename Thermo, typename Transport, typename Kinetics>
  class GasMixture
  {
  public:

    GasMixture( const GetPot& input );
    ~GasMixture();

    const Chemistry& get_chemistry() const;

    const Thermo& get_thermo() const;

    const Transport& get_transport() const;

    const Kinetics& get_kinetics() const;

  protected:

    Chemistry _chemistry;

    Thermo _thermo;

    Transport _transport;

    Kinetics _kinetics;

  private:

    GasMixture();

  };

  /* ------------------------- Inline Functions -------------------------*/
  template<typename C, typename Th, typename Tr, typename K>
  inline
  const Chemistry& GasMixture<C,Th,Tr,K>::get_chemistry() const
  {
    return _chemistry;
  }

  template<typename C, typename Th, typename Tr, typename K>
  inline
  const Thermo& GasMixture<C,Th,Tr,K>::get_thermo() const
  {
    return _thermo;
  }

  template<typename C, typename Th, typename Tr, typename K>
  inline
  const Transport& GasMixture<C,Th,Tr,K>::get_transport() const
  {
    return _transport;
  }

  template<typename C, typename Th, typename Tr, typename K>
  inline
  const Kinetics& GasMixture<C,Th,Tr,K>::get_kinetics() const
  {
    return _transport;
  }

} // end namespace GRINS

#endif // GRINS_GAS_MIXTURE_H
