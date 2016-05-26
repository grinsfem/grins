//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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

#ifndef GRINS_FACTORY_ABSTRACT_H
#define GRINS_FACTORY_ABSTRACT_H

// C++
#include <string>
#include <map>

// libMesh
#include "libmesh/auto_ptr.h" // libMesh::UniquePtr

namespace GRINS
{
  //! Abstract factory.
  /*! Copied from libMesh::Factory. The main difference is the helper method
      to fetch the factory (this was moved out of the build() function in libMesh::Factory).
      That is useful so subclasses can additional static methods
      that dispatch to virtual methods of the subclass. */
  template<typename Base>
  class FactoryAbstract
  {
  public:

    virtual ~FactoryAbstract(){}

    //! Use this method to build objects of type Base.
    static libMesh::UniquePtr<Base> build(const std::string & name);

    //! Subclasses implement the actual construction of the Base object in create().
    virtual libMesh::UniquePtr<Base> create () = 0;

  protected:

    //! Constructor is protected. Use the build() method to construct Base objects.
    FactoryAbstract(const std::string & name);

    //! Helper method that looks up the factory and returns it if present, or error if it's not.
    static FactoryAbstract<Base>& get_factory(const std::string & name);

    //! Like get_factory, but will downcast to DerivedType
    template<typename DerivedType>
    static DerivedType& get_factory_subclass(const std::string & name);

    static std::map<std::string, FactoryAbstract<Base>*>& factory_map();

  };

  template <class Base>
  inline
  FactoryAbstract<Base>::FactoryAbstract(const std::string & name)
  {
    // Make sure we haven't already added this name
    // to the map
    libmesh_assert (!factory_map().count(name));

    factory_map()[name] = this;
  }

  template <class Base>
  inline
  libMesh::UniquePtr<Base> FactoryAbstract<Base>::build(const std::string & name)
  {
    FactoryAbstract<Base>& factory = get_factory(name);
    return factory.create();
  }

  template <class Base>
  inline
  FactoryAbstract<Base>& FactoryAbstract<Base>::get_factory(const std::string & name)
  {
    if (!factory_map().count(name))
      {
        std::stringstream ss;
        ss << "ERROR: Tried to build an unknown FactoryAbstract type: "
           << name << std::endl
           << "valid options are:" << std::endl;

        for( typename std::map<std::string,FactoryAbstract<Base>*>::const_iterator
               it = factory_map().begin(); it != factory_map().end(); ++it )
          ss << "  " << it->first << std::endl;

        libmesh_error_msg(ss.str());
      }

    FactoryAbstract<Base>* factory = factory_map()[name];
    return *factory;
  }

  template <class Base>
  template<typename DerivedType>
  inline
  DerivedType& FactoryAbstract<Base>::get_factory_subclass(const std::string & name)
  {
    FactoryAbstract<Base>& raw_factory = get_factory(name);
    DerivedType& factory = libMesh::cast_ref<DerivedType&>(raw_factory);
    return factory;
  }

} // end namespace GRINS

#endif // GRINS_FACTORY_ABSTRACT_H
