/**
 * @brief this file provides definitions for the off-topic internal function
 *        definitions.
 */

#include "data.hpp"
#include "internal.hpp"

namespace cuda
{
   namespace internal
   {
      /*
       * Initlialization functions
       */

      bool __initialize_atoms ();
      bool __initialize_fields ();
      bool __initialize_cells ();
      bool __initialize_materials ();
      bool __initialize_topology ();

   } /* internal */
} /* cuda */
