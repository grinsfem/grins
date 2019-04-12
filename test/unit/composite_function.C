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

// This test is in libMesh now, but let's leave it in GRINS until
// we're sure that's working.

// libMesh
#include "libmesh/composite_function.h"
#include "libmesh/const_function.h"
#include "libmesh/dense_vector.h"

using namespace libMesh;

int main( /*int argc, char* argv[]*/ )
{
  std::vector<std::vector<unsigned int> > index_sets(4);
  index_sets[0].resize(2);
  index_sets[0][0] = 3;
  index_sets[0][1] = 4;
  index_sets[1].resize(3);
  index_sets[1][0] = 0;
  index_sets[1][1] = 1;
  index_sets[1][2] = 2;
  index_sets[2].resize(3);
  index_sets[2][0] = 0;
  index_sets[2][1] = 2;
  index_sets[2][2] = 4;
  index_sets[3].resize(5);
  index_sets[3][0] = 5;
  index_sets[3][1] = 1;
  index_sets[3][2] = 3;
  index_sets[3][3] = 6;
  index_sets[3][4] = 7;

  CompositeFunction<libMesh::Real> composite_outer;

  {
    CompositeFunction<libMesh::Real> composite_inner;
    composite_inner.attach_subfunction
      (ConstFunction<libMesh::Real>(1), index_sets[0]);
    composite_inner.attach_subfunction
      (ConstFunction<libMesh::Real>(2), index_sets[1]);
    composite_outer.attach_subfunction
      (composite_inner, index_sets[3]);

    libMesh::DenseVector<libMesh::Real> test_one(5);

    composite_inner(libMesh::Point(0), 0, test_one);

    if (test_one(0) != 2)
      return 1;
    if (test_one(1) != 2)
      return 1;
    if (test_one(2) != 2)
      return 1;
    if (test_one(3) != 1)
      return 1;
    if (test_one(4) != 1)
      return 1;
  }
  composite_outer.attach_subfunction
    (ConstFunction<libMesh::Real>(3), index_sets[2]);

  libMesh::DenseVector<libMesh::Real> test_two(8);
  composite_outer(libMesh::Point(0), 0, test_two);

  if (test_two(0) != 3)
    return 1;
  if (test_two(2) != 3)
    return 1;
  if (test_two(4) != 3)
    return 1;
  if (test_two(5) != 2)
    return 1;
  if (test_two(1) != 2)
    return 1;
  if (test_two(3) != 2)
    return 1;
  if (test_two(6) != 1)
    return 1;
  if (test_two(7) != 1)
    return 1;

  return 0;
}
