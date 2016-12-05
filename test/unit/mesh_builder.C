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

#include "grins_config.h"

#ifdef GRINS_HAVE_CPPUNIT

#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include "test_comm.h"
#include "grins_test_paths.h"

// GRINS
#include "grins/grins_enums.h"
#include "grins/mesh_builder.h"

// libMesh
#include "libmesh/elem.h"

// Ignore warnings from auto_ptr in CPPUNIT_TEST_SUITE_END()
#include <libmesh/ignore_warnings.h>

namespace GRINSTesting
{
  class MeshBuilderTest : public CppUnit::TestCase
  {
  public:
    CPPUNIT_TEST_SUITE( MeshBuilderTest );

    CPPUNIT_TEST( test_build_1d_mesh );
    CPPUNIT_TEST( test_build_2d_mesh );
    CPPUNIT_TEST( test_build_3d_mesh );

    CPPUNIT_TEST_SUITE_END();

  public:

    void test_build_1d_mesh()
    {
      std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/mesh_build_1d.in";
      GetPot input(filename);
      GRINS::SharedPtr<libMesh::UnstructuredMesh> mesh = this->build_mesh(input);
      CPPUNIT_ASSERT_EQUAL((libMesh::dof_id_type)22,mesh->n_elem());
      this->test_elem_type(*mesh,GRINSEnums::EDGE2);
    }

    void test_build_2d_mesh()
    {
      std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/mesh_build_2d.in";
      GetPot input(filename);
      GRINS::SharedPtr<libMesh::UnstructuredMesh> mesh = this->build_mesh(input);
      CPPUNIT_ASSERT_EQUAL((libMesh::dof_id_type)100,mesh->n_elem());
      this->test_elem_type(*mesh,GRINSEnums::QUAD9);

    }

    void test_build_3d_mesh()
    {
      std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/mesh_build_3d.in";
      GetPot input(filename);
      GRINS::SharedPtr<libMesh::UnstructuredMesh> mesh = this->build_mesh(input);
      CPPUNIT_ASSERT_EQUAL((libMesh::dof_id_type)125,mesh->n_elem());
      this->test_elem_type(*mesh,GRINSEnums::HEX8);
    }

  private:

    GRINS::SharedPtr<libMesh::UnstructuredMesh> build_mesh( const GetPot& input )
    {
      GRINS::MeshBuilder mesh_builder;
      return mesh_builder.build( input, *TestCommWorld );
    }

    void test_elem_type( const libMesh::MeshBase& mesh, GRINSEnums::ElemType elem_type_expected )
    {
      for( libMesh::MeshBase::const_element_iterator e = mesh.active_elements_begin();
           e != mesh.active_elements_end(); ++ e )
        {
          const libMesh::Elem* elem = *e;
          GRINSEnums::ElemType elem_type_computed = elem->type();
          CPPUNIT_ASSERT_EQUAL( elem_type_expected, elem_type_computed);
        }
    }
  };

  CPPUNIT_TEST_SUITE_REGISTRATION( MeshBuilderTest );

} // end namespace GRINSTesting

#endif // GRINS_HAVE_CPPUNIT
