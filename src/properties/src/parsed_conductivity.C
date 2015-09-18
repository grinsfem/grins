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

  ParsedConductivity::~ParsedConductivity()
  {
    return;
  }

} // namespace GRINS
