//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef GRINS_POSTPROCESSING_FACTORY_H
#define GRINS_POSTPROCESSING_FACTORY_H

#include "boost/tr1/memory.hpp"

// libMesh
#include "libmesh/getpot.h"

// GRINS
#include "grins/postprocessed_quantities.h"

namespace GRINS
{

  //! This object handles constructing the postprocessing object to be used.
  /*! To allow the user to easily extend the postprocesing capabilities,
      the postprocessing construction is handled in this object. */
  class PostprocessingFactory
  {
  public:
    
    PostprocessingFactory();
    virtual ~PostprocessingFactory();

    virtual std::tr1::shared_ptr<PostProcessedQuantities<Real> > build(const GetPot& input);

  };

} // namespace GRINS

#endif // GRINS_POSTPROCESSING_FACTORY_H
