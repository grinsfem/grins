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

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/parsed_function.h"

namespace GRINS
{

  ParsedConductivity::ParsedConductivity( const GetPot& input )
    : ParsedPropertyBase(),
      ParameterUser("ParsedConductivity")
  {
    if( !input.have_variable("Materials/Conductivity/k") )
      {
        libmesh_error_msg("ERROR: No conductivity has been specified!");
      }
    else
      {
        std::string conductivity_function = input("Materials/Conductivity/k",std::string("0"));

        if( !this->check_func_nonzero(conductivity_function) )
          {
            libmesh_error_msg("ERROR: Detected '0' function for ParsedConductivity!");
          }

        this->_func.reset(new libMesh::ParsedFunction<libMesh::Number>(conductivity_function));
      }
  }

  ParsedConductivity::~ParsedConductivity()
  {
    return;
  }

} // namespace GRINS
