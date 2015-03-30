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
#include "grins/parsed_viscosity.h"

//GRINS
#include "grins/grins_physics_names.h"
#include "grins/common.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/parsed_function.h"

namespace GRINS
{

   ParsedViscosity::ParsedViscosity( const GetPot& input ) :
     ParameterUser("ParsedViscosity")
    {

      // Warning about this constructor being deprecated
      {
        std::string warning = "WARNING: Use of this constructor is DEPRECATED.\n";
        warning += "         Please update to use constructor with input material name.\n";
        grins_warning(warning);
      }

      if( !input.have_variable("Materials/Viscosity/mu") )
       {
         std::cerr<<"No viscosity has been specified."<<std::endl;
      
         libmesh_error();
       }
      else
       {
         std::string viscosity_function = input("Materials/Viscosity/mu",std::string("0"));

         mu.reset(new libMesh::ParsedFunction<libMesh::Number>(viscosity_function));

         if (viscosity_function == "0")
            {
              std::cerr << "Warning! Zero Viscosity specified!" << std::endl;

              libmesh_error();
            }
       }

      return;
      }

  ParsedViscosity::~ParsedViscosity()
  {
    return;
  }

} // namespace GRINS
