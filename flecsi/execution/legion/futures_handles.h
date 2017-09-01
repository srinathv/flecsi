/*~--------------------------------------------------------------------------~*
*  @@@@@@@@  @@           @@@@@@   @@@@@@@@ @@
* /@@/////  /@@          @@////@@ @@////// /@@
* /@@       /@@  @@@@@  @@    // /@@       /@@
* /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@
* /@@////   /@@/@@@@@@@/@@       ////////@@/@@
* /@@       /@@/@@//// //@@    @@       /@@/@@
* /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@
* //       ///  //////   //////  ////////  //
*
* Copyright (c) 2016 Los Alamos National Laboratory, LLC
* All rights reserved
*~--------------------------------------------------------------------------~*/

#ifndef flecsi_execution_legion_futures_handles_h
#define flecsi_execution_legion_futures_handles_h

//----------------------------------------------------------------------------//
//! @file
//! @date Initial file creation: Sep 1, 2017
//----------------------------------------------------------------------------//

#include <vector>
#include <type_traits>

#include "legion.h"
#include "arrays.h"

#include "flecsi/data/common/privilege.h"
#include "flecsi/utils/tuple_walker.h"
#include "flecsi/data/data_client_handle.h"
#include "flecsi/topology/mesh_types.h"

namespace flecsi {
namespace execution {

//----------------------------------------------------------------------------//
//! The futures_handles_t type can be called to walk task args after task
//! launch. This allows us to map physical regions to internal handle
//! buffers/accessors.
//!
//! @ingroup execution
//----------------------------------------------------------------------------//

struct futures_handles_t : public utils::tuple_walker__<futures_handles_t>
{

  //--------------------------------------------------------------------------//
  //! Construct an futures_handles_t instance.
  //!
  //! @param runtime The Legion task runtime.
  //! @param context The Legion task runtime context.
  //--------------------------------------------------------------------------//

  futures_handles_t(
    Legion::Runtime* runtime,
    Legion::Context& context,
    const std::vector<Legion::Future>& futures
  )
  :
    runtime(runtime),
    context(context),
    futures(futures),
    future(0)
  {
  } // futures_handles


  //-----------------------------------------------------------------------//
  // If this is a data handle, then simply skip it.
  //-----------------------------------------------------------------------//

  template<
    typename T,
    size_t EXCLUSIVE_PERMISSIONS,
    size_t SHARED_PERMISSIONS,
    size_t GHOST_PERMISSIONS
  >
  void
  handle(
    data_handle__<
      T,
      EXCLUSIVE_PERMISSIONS,
      SHARED_PERMISSIONS,
      GHOST_PERMISSIONS
    > & h
  )
  {
  } // handle

  template<
    typename T,
    size_t PERMISSIONS
  >
  void
  handle(
    data_client_handle__<T, PERMISSIONS> & h
  )
  {
  } // handle

  template<
    typename T
  >
  typename std::enable_if_t<!std::is_base_of<data_handle_base_t, T>::value>
  handle(
    T &
  )
  {
    future++;
  } // handle

  Legion::Runtime * runtime;
  Legion::Context & context;
  const std::vector<Legion::Future> & futures;
  size_t future;
}; // struct futures_handles_t

} // namespace execution 
} // namespace flecsi

#endif // flecsi_execution_legion_futures_handles_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
