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

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{

  ParsedViscosity::ParsedViscosity( const GetPot& input )
    : _mu( input( "Materials/Viscosity/mu", 1.0 ) )
  {
    if( !input.have_variable("Materials/Viscosity/mu") )
      {
        libmesh_warning("No Materials/Viscosity/mu specified!\n");

	// Try and get the viscosity from other specifications
	_mu = input("Physics/"+incompressible_navier_stokes+"/mu", 1.0);
	
      }

    return;
  }

  ParsedViscosity::~ParsedViscosity()
  {
    return;
  }

} // namespace GRINS
