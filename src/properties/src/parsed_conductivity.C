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

namespace GRINS
{

  ParsedConductivity::ParsedConductivity( const GetPot& input )
    : ParsedPropertyBase(),
      ParameterUser("ParsedConductivity")
  {
    this->set_parameter(this->_func, input,
                        "Materials/Conductivity/k",
                        "DIE!");
  }

  ParsedConductivity::~ParsedConductivity()
  {
    return;
  }

} // namespace GRINS
