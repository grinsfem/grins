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
#include "grins/grins_physics_names.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/parsed_function.h"

namespace GRINS
{

   ParsedConductivity::ParsedConductivity( const GetPot& input )    
    {
      if( !input.have_variable("Materials/Conductivity/k") )
       {
         std::cerr<<"No conductivity has been specified."<<std::endl;
      
         libmesh_error();
       }
      else
       {
         std::string conductivity_function = input("Materials/Conductivity/k",std::string("0"));

         k.reset(new libMesh::ParsedFunction<libMesh::Number>(conductivity_function));

         if (conductivity_function == "0")
            {
              std::cerr << "Warning! Zero Conductivity specified!" << std::endl;

              libmesh_error();
            }
       }

      return;
      }

  ParsedConductivity::~ParsedConductivity()
  {
    return;
  }

} // namespace GRINS
