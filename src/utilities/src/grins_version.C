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


#include "grins/grins_version.h"

namespace GRINS
{

  void grins_version_stdout()
  {
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "GRINS Package: Version = " << GRINS_LIB_VERSION;
    std::cout << " (" << get_grins_version() << ")" << std::endl << std::endl;

    std::cout << GRINS_LIB_RELEASE << std::endl << std::endl;

    std::cout << "Build Date   = " << GRINS_BUILD_DATE     << std::endl;
    std::cout << "Build Host   = " << GRINS_BUILD_HOST     << std::endl;
    std::cout << "Build User   = " << GRINS_BUILD_USER     << std::endl;
    std::cout << "Build Arch   = " << GRINS_BUILD_ARCH     << std::endl;
    std::cout << "Build Rev    = " << GRINS_BUILD_VERSION  << std::endl << std::endl;

    std::cout << "C++ Config   = " << GRINS_CXX << " " << GRINS_CXXFLAGS << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;

    return;
  }

  int get_grins_version()
  {
    /* Note: return format follows the versioning convention xx.yy.zz where

       xx = major version number
       yy = minor version number
       zz = micro version number

       For example:
       v.   0.23  -> 002300 = 2300
       v   0.23.1 -> 002301 = 2301
       v. 10.23.2 -> 102302         */

    int major_version = 0;
    int minor_version = 0;
    int micro_version = 0;

#ifdef GRINS_MAJOR_VERSION
    major_version = GRINS_MAJOR_VERSION;
#endif

#ifdef GRINS_MINOR_VERSION
    minor_version = GRINS_MINOR_VERSION;
#endif

#ifdef GRINS_MICRO_VERSION
    micro_version = GRINS_MICRO_VERSION;
#endif

    return major_version*10000 + minor_version*100 + micro_version;
  }

} // namespace GRINS
