/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#ifndef flecsi_partition_h
#define flecsi_partition_h

//----------------------------------------------------------------------------//
//! @file
//! @date Initial file creation: May 19, 2017
//----------------------------------------------------------------------------//

namespace flecsi {

  enum partition_t : size_t{
    exclusive = 0b001,
    shared    = 0b010,
    ghost     = 0b100,
    owned     = 0b011
  };

} // namespace flecsi

#endif // flecsi_partition_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
