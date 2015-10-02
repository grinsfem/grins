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
#include "grins/parsed_conductivity.h"

//GRINS
#include "grins/common.h"
#include "grins/grins_physics_names.h"
#include "grins/materials_parsing.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{

  ParsedConductivity::ParsedConductivity( const GetPot& input )
    : ParsedPropertyBase(),
      ParameterUser("ParsedConductivity")
  {
    // Warning about this constructor being deprecated
    {
      std::string warning = "WARNING: Use of this constructor is DEPRECATED.\n";
      warning += "         Please update to use constructor with input material name.\n";
      grins_warning(warning);
    }

    this->set_parameter(this->_func, input,
                        "Materials/Conductivity/k",
                        "DIE!");

    std::string conductivity_function = input("Materials/Conductivity/k",std::string("0"));

    if( !this->check_func_nonzero(conductivity_function) )
      {
        libmesh_error_msg("ERROR: Detected '0' function for ParsedConductivity!");
      }
  }

  ParsedConductivity::ParsedConductivity( const GetPot& input, const std::string& material )
    : ParsedPropertyBase(),
      ParameterUser("ParsedConductivity")
  {
    std::string conductivity_function;

    // If we have the new version, we parse that
    if( input.have_variable("Materials/"+material+"/ThermalConductivity/value") )
      {
        this->set_parameter(this->_func, input,
                            "Materials/"+material+"/ThermalConductivity/value",
                            "DIE!");

        conductivity_function = input("Materials/"+material+"/ThermalConductivity/value",std::string("0"));
      }
    // If we have the old DEPRECATED version, use that
    else if( input.have_variable("Materials/Conductivity/k") )
      {
        MaterialsParsing::dep_input_warning( "Materials/Conductivity/k",
                                             "ThermalConductivity" );

        this->set_parameter(this->_func, input,
                            "Materials/Conductivity/k",
                            "DIE!");

        conductivity_function = input("Materials/Conductivity/k",std::string("0"));
      }
    // If we don't have either, that's an error
    else
      {
        libmesh_error_msg("Error: Could not find either Materials/"+material+"/ThermalConductivity/value or Materials/Conductivity/k");
      }

    if( !this->check_func_nonzero(conductivity_function) )
      {
        libmesh_error_msg("ERROR: Detected '0' function for ParsedConductivity!");
      }
  }

  ParsedConductivity::~ParsedConductivity()
  {
    return;
  }

} // namespace GRINS
