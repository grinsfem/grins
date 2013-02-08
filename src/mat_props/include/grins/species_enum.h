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

#ifndef GRINS_SPECIES_ENUM_H
#define GRINS_SPECIES_ENUM_H

namespace GRINS
{
  /*!
   * Simple species enumeration list.  By convention, 
   * standard periodic-table names are used, and for 
   * charged particles 'm' replaces '-', and 'p' 
   * replaces '+'.  For example, 'NO+' is represented
   * as 'NOp'. Originally taken from FIN-S.
   */
  enum Species 
    {
      Air = 0,  // Thermally perfect air (Cv=Cv(T))
      CPAir,    // Calorically perfect air (Cv=constant)
      Ar,  Arp, // Argon and its positively charged ion
      C,   Cp,  // Atomic Carbon and its positively charged ion  	
      C2,   	
      C2H,  	
      C2H2, 	
      C3,    	
      CF,   	
      CF2,  	
      CF3,  	
      CF4,  	
      CH,   	
      CH2,  	
      CH3,  	
      CH4,  	
      Cl,    	
      Cl2,   	
      CN,  CNp,  	
      CO,  COp,  	
      CO2,  	
      F,    	 	
      F2,   	
      H,  Hp,   	
      H2, H2p,  	
      H2O,
      H2O2,
      HCl,   	
      HCN,  	
      He, Hep, // Helium and its positively charged ion
      HO2,
      N,  Np,  // Atomic Nitrogen and its positively charged ion 	 	
      N2, N2p, // Molecular Nitrogen and its positively charged ion 	       	
      CPN2,    // Calorically perfect N2 (Cv=constant)
      Ne,   	 // Neon
      NCO,  	
      NH, NHp,  	
      NH2,  	
      NH3,  	
      NO, NOp, // Nitric Oxide and its positively charged ion
      NO2,  	
      O,  Op,  // Atomic Oxygen and its positively charged ion
      O2, O2p, // Molecular Oxygen and its positively charged ion	
      OH,   	
      Si,   	
      SiO,  	
      e        // Free electron
    };

} // namespace GRINS

#endif // GRINS_SPECIES_ENUM_H
