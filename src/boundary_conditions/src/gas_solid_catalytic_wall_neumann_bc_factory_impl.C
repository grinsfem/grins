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

// This class
#include "grins/gas_solid_catalytic_wall_neumann_bc_factory_impl.h"

// GRINS
#include "grins/string_utils.h"
#include "grins/chemistry_builder.h"

#ifdef GRINS_HAVE_CANTERA
#include "grins/cantera_mixture.h"
#endif

#ifdef GRINS_HAVE_ANTIOCH
#include "grins/antioch_chemistry.h"
#endif

namespace GRINS
{
  // To avoid compiler warnings without GRINS or Cantera
#if defined(GRINS_HAVE_ANTIOCH) || defined(GRINS_HAVE_CANTERA)
  std::shared_ptr<NeumannBCAbstract>
  GasSolidCatalyticWallNeumannBCFactoryImpl::build_catalytic_wall
  ( const GetPot& input, const std::string& reaction,std::shared_ptr<CatalycityBase>& gamma_ptr,
    const std::vector<VariableIndex>& species_vars,const std::string& material,
    VariableIndex T_var,libMesh::Real p0,const std::string& thermochem_lib )
#else
    std::shared_ptr<NeumannBCAbstract>
    GasSolidCatalyticWallNeumannBCFactoryImpl::build_catalytic_wall
                                              ( const GetPot& /*input*/, const std::string& reaction,std::shared_ptr<CatalycityBase>& /*gamma_ptr*/,
                                                const std::vector<VariableIndex>& /*species_vars*/,const std::string& /*material*/,
                                                VariableIndex /*T_var*/,libMesh::Real /*p0*/,const std::string& thermochem_lib )
#endif
  {
    std::string gas_reactant;
    std::string solid_reactant;
    std::string product;
    this->parse_reactants_and_product(reaction,gas_reactant,solid_reactant,product);

    // Now construct the Neumann BC
    std::shared_ptr<NeumannBCAbstract> catalytic_wall;

    ChemistryBuilder chem_builder;

    if( thermochem_lib == "cantera" )
      {
#ifdef GRINS_HAVE_CANTERA
        std::unique_ptr<CanteraMixture> chem_uptr;
        chem_builder.build_chemistry(input,material,chem_uptr);

        /*! \todo Update the API for the catalytic walls to take a unique_ptr to avoid this garbage.*/
        std::shared_ptr<CanteraMixture> chem_ptr(chem_uptr.release());

        this->build_wall_ptr<CanteraMixture>(chem_ptr,gamma_ptr,gas_reactant,solid_reactant,
                                             product,species_vars,T_var,p0,catalytic_wall);
#else
        libmesh_error_msg("Error: Cantera not enabled in this configuration. Reconfigure using --with-cantera option.");
#endif
      }
    else if( thermochem_lib == "antioch" )
      {
#ifdef GRINS_HAVE_ANTIOCH
        std::unique_ptr<AntiochChemistry> chem_uptr;
        chem_builder.build_chemistry(input,material,chem_uptr);

        /*! \todo Update the API for the catalytic walls to take a unique_ptr to avoid this garbage.*/
        std::shared_ptr<AntiochChemistry> chem_ptr(chem_uptr.release());

        this->build_wall_ptr<AntiochChemistry>(chem_ptr,gamma_ptr,gas_reactant,solid_reactant,
                                               product,species_vars,T_var,p0,catalytic_wall);
#else
        libmesh_error_msg("Error: Antioch not enabled in this configuration. Reconfigure using --with-antioch option.");
#endif
      }
    else
      libmesh_error_msg("ERROR: Invalid thermochemistry library "+thermochem_lib+"!");

    return catalytic_wall;
  }

  void GasSolidCatalyticWallNeumannBCFactoryImpl::parse_reactants_and_product( const std::string& reaction,
                                                                               std::string& gas_reactant,
                                                                               std::string& solid_reactant,
                                                                               std::string& product ) const
  {
    /* We are expecting reactions of the form
       X+Y(s)->Z  or
       Y(s)+X->X
       So, first we'll split on the "->", then split the reactants up and
       figure out which is the gas species and which is the solid species. */
    std::vector<std::string> partners;
    StringUtilities::split_string(reaction, "->", partners);

    const std::string pre_split_reactants = partners[0];
    product = partners[1];

    std::vector<std::string> split_reactants;
    StringUtilities::split_string(pre_split_reactants, "+", split_reactants);

    // We can only handle two reactants currently
    if( split_reactants.size() != 2 )
      {
        std::string error_msg = "ERROR: Currently, GasSolidCatalyticWall boundary condition only supports\n";
        error_msg += "       reactions of the form X+Y(s)->Z or Y(s)+X->X. Found ";
        error_msg += StringUtilities::T_to_string<unsigned int>(split_reactants.size())+" reactants.\n";
        libmesh_error_msg(error_msg);
      }


    // Check if the first reactant is the solid one
    if( split_reactants[0].find("(s)") == split_reactants[0].npos )
      {
        // If not found, check the second reactant and error out if not found.
        if( split_reactants[1].find("(s)") == split_reactants[1].npos )
          {
            std::string error_msg = "ERROR: could not find solid reactant for GasSolidCatalyticWall!\n";
            error_msg += "       Found reactants "+split_reactants[0]+", "+split_reactants[1]+"\n";
            libmesh_error_msg(error_msg);
          }
        else
          {
            gas_reactant = split_reactants[0];
            solid_reactant = split_reactants[1].substr(0,split_reactants[1].find("(s)"));
          }
      }
    // Found (s) in the first reactant
    else
      {
        // Check that there's not 2 solid reactants
        if( split_reactants[1].find("(s)") != split_reactants[1].npos )
          {
            std::string error_msg = "ERROR: can have only one solid reactant for GasSolidCatalyticWall!\n";
            error_msg += "       Found reactants "+split_reactants[0]+", "+split_reactants[1]+"\n";
            libmesh_error_msg(error_msg);
          }

        gas_reactant = split_reactants[1];
        solid_reactant = split_reactants[0].substr(0,split_reactants[0].find("(s)"));
      }
  }

} // end namespace GRINS
