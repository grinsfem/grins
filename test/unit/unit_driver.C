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
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
#include <libmesh/restore_warnings.h>
#endif // GRINS_HAVE_CPPUNIT

#include <libmesh/libmesh.h>

#include "test_comm.h"

int main(int argc, char **argv)
{
#ifdef GRINS_HAVE_CPPUNIT

  // Initialize the library.  This is necessary because the library
  // may depend on a number of other libraries (i.e. MPI  and Petsc)
  // that require initialization before use.
  libMesh::LibMeshInit init(argc, argv);
  TestCommWorld = &init.comm();

  CppUnit::TextUi::TestRunner runner;
  CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();
  runner.addTest( registry.makeTest() );

  // If the tests all succeed, report success
  if (runner.run())
    return 0;

  // If any test fails report failure
  return 1;

#else
  // If we don't have CPPUnit, report we skipped
  // 77 return code tells Automake we skipped this.
  return 77;
#endif // GRINS_HAVE_CPPUNIT
}

#ifdef GRINS_HAVE_CPPUNIT
libMesh::Parallel::Communicator *TestCommWorld;
#endif // GRINS_HAVE_CPPUNIT
