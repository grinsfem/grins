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

#ifndef GRINS_CATALYCITY_FACTORY_ABSTRACT_H
#define GRINS_CATALYCITY_FACTORY_ABSTRACT_H

// GRINS
#include "grins/factory_with_getpot.h"
#include "grins/catalycity_base.h"

namespace GRINS
{
  // According to the standard, we need a declaration of the
  // specialization which precedes any automatic instantiation.
  template<> const GetPot* FactoryWithGetPot<CatalycityBase>::_input;

  class CatalycityFactoryAbstract : public FactoryWithGetPot<CatalycityBase>
  {
  public:
    CatalycityFactoryAbstract( const std::string& physics_name )
      : FactoryWithGetPot<CatalycityBase>(physics_name)
    {}

    ~CatalycityFactoryAbstract(){};

    static void set_section( const std::string& section )
    { _section = section; }

  protected:

    virtual std::unique_ptr<CatalycityBase> build_catalycity( const GetPot& input,
                                                              const std::string& section ) =0;

    //! Helper function to reduce code duplication
    virtual void check_state() const;

    //! Helper function to reduce code duplication
    virtual void reset_state();

    static std::string _section;

  private:

    virtual std::unique_ptr<CatalycityBase> create();
  };

  inline
  std::unique_ptr<CatalycityBase> CatalycityFactoryAbstract::create()
  {
    this->check_state();

    std::unique_ptr<CatalycityBase> new_catalycity = this->build_catalycity( *_input, _section );

    this->reset_state();

    return new_catalycity;
  }

  inline
  void CatalycityFactoryAbstract::check_state() const
  {
    if( !_input )
      libmesh_error_msg("ERROR: must call set_getpot() before building Catalycity!");

    if( _section == std::string("DIE!") )
      libmesh_error_msg("ERROR: must call set_section() before building Catalycity!");
  }

  inline
  void CatalycityFactoryAbstract::reset_state()
  {
    // Reset for error checking
    _section = std::string("DIE!");
  }

} // end namespace GRINS

#endif // GRINS_CATALYCITY_FACTORY_ABSTRACT_H
