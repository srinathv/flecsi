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

#ifndef flecsi_execution_legion_task_prolog_h
#define flecsi_execution_legion_task_prolog_h

//----------------------------------------------------------------------------//
//! @file
//! @date Initial file creation: May 19, 2017
//----------------------------------------------------------------------------//

#include <legion.h>
#include <vector>

#include "flecsi/data/data.h"
#include "flecsi/execution/context.h"
#include "flecsi/execution/legion/internal_field.h"

clog_register_tag(prolog);

namespace flecsi {
namespace execution {

  //--------------------------------------------------------------------------//
  //! The task_prolog_t type can be called to walk the task args after the
  //! task launcher is created, but before the task has run. This allows
  //! synchronization dependencies to be added to the execution flow.
  //!
  //! @ingroup execution
  //--------------------------------------------------------------------------//

  struct task_prolog_t : public utils::tuple_walker__<task_prolog_t>
  {

    //------------------------------------------------------------------------//
    //! Construct a task_prolog_t instance.
    //!
    //! @param runtime The Legion task runtime.
    //! @param context The Legion task runtime context.
    //------------------------------------------------------------------------//

    task_prolog_t(
      Legion::Runtime * runtime,
      Legion::Context & context,
      Legion::TaskLauncher & launcher
    )
    :
      runtime(runtime),
      context(context),
      launcher(launcher)
    {
    } // task_prolog_t

    //------------------------------------------------------------------------//
    //! FIXME: Need a description.
    //!
    //! @tparam T                     The data type referenced by the handle.
    //! @tparam EXCLUSIVE_PERMISSIONS The permissions required on the exclusive
    //!                               indices of the index partition.
    //! @tparam SHARED_PERMISSIONS    The permissions required on the shared
    //!                               indices of the index partition.
    //! @tparam GHOST_PERMISSIONS     The permissions required on the ghost
    //!                               indices of the index partition.
    //!
    //! @param runtime The Legion task runtime.
    //! @param context The Legion task runtime context.
    //------------------------------------------------------------------------//

    template<
      typename T,
      size_t EXCLUSIVE_PERMISSIONS,
      size_t SHARED_PERMISSIONS,
      size_t GHOST_PERMISSIONS
    >
    void
    handle_old(
      data_handle__<
        T,
        EXCLUSIVE_PERMISSIONS,
        SHARED_PERMISSIONS,
        GHOST_PERMISSIONS
      > & h
    )
    {
 if (!h.global && !h.color){
      auto& flecsi_context = context_t::instance();

      bool read_phase = false;
      bool write_phase = false;
      const int my_color = runtime->find_local_MPI_rank();

      read_phase = GHOST_PERMISSIONS != reserved;
      write_phase = (SHARED_PERMISSIONS == wo) || (SHARED_PERMISSIONS == rw);

      if(read_phase) {
        {
        clog_tag_guard(prolog);
        clog(trace) << "rank " << my_color << " READ PHASE PROLOGUE" <<
          std::endl;

        // As owner
        clog(trace) << "rank " << my_color << " arrives & advances " <<
          *(h.pbarrier_as_owner_ptr) << std::endl;
        } // scope

        // Phase WRITE
        h.pbarrier_as_owner_ptr->arrive(1);

        // Phase WRITE
        *(h.pbarrier_as_owner_ptr) = runtime->advance_phase_barrier(context,
          *(h.pbarrier_as_owner_ptr));

        const size_t _pbp_size = h.ghost_owners_pbarriers_ptrs.size();

        // As user
        for(size_t owner{0}; owner<_pbp_size; owner++) {

          {
          clog_tag_guard(prolog);
          clog(trace) << "rank " << my_color << " WAITS " <<
            *(h.ghost_owners_pbarriers_ptrs[owner]) << std::endl;

          clog(trace) << "rank " << my_color << " arrives & advances " <<
            *(h.ghost_owners_pbarriers_ptrs[owner]) << std::endl;
          } // scope

          Legion::RegionRequirement rr_shared(h.ghost_owners_lregions[owner],
            READ_ONLY, EXCLUSIVE, h.ghost_owners_lregions[owner]);

          Legion::RegionRequirement rr_ghost(h.ghost_lr,
            WRITE_DISCARD, EXCLUSIVE, h.color_region);

          // auto iitr = flecsi_context.field_info_map().find(
          //   { h.data_client_hash, h.index_space });

          // clog_assert(iitr != flecsi_context.field_info_map().end(),
          //   "invalid index space");

          auto ghost_owner_pos_fid = LegionRuntime::HighLevel::FieldID(
            internal_field::ghost_owner_pos);

          rr_ghost.add_field(ghost_owner_pos_fid);

          // for(auto& fitr: iitr->second) {
          //   const context_t::field_info_t & fi = fitr.second;
          //   if(!utils::hash::is_internal(fi.key)){
          //     rr_shared.add_field(fi.fid);
          //     rr_ghost.add_field(fi.fid);
          //   }
          // } // for

          rr_shared.add_field(h.fid);
          rr_ghost.add_field(h.fid);

          // TODO - circular dependency including internal_task.h
          auto constexpr key = flecsi::utils::const_string_t{
            EXPAND_AND_STRINGIFY(ghost_copy_task)}.hash();

          const auto ghost_copy_tid = flecsi_context.task_id<key>();

          // Anonymous struct for arguments in task launcer below.
          struct {
            size_t data_client_hash;
            size_t index_space;
            size_t owner;
          } args;

          args.data_client_hash = h.data_client_hash;
          args.index_space = h.index_space;
          args.owner = owner;

          Legion::TaskLauncher launcher(ghost_copy_tid,
            Legion::TaskArgument(&args, sizeof(args)));

          {
          clog_tag_guard(prolog);
          clog(trace) << "gid to lid map size = " <<
            h.global_to_local_color_map_ptr->size() << std::endl;
          } // scope

          launcher.add_future(Legion::Future::from_value(runtime,
            *(h.global_to_local_color_map_ptr)));

          // Phase READ
          launcher.add_region_requirement(rr_shared);
          launcher.add_region_requirement(rr_ghost);
          launcher.add_wait_barrier(*(h.ghost_owners_pbarriers_ptrs[owner]));

          // Phase WRITE
          launcher.add_arrival_barrier(*(h.ghost_owners_pbarriers_ptrs[owner]));

          // Execute the ghost copy task
          runtime->execute_task(context, launcher);

          // Phase WRITE
          *(h.ghost_owners_pbarriers_ptrs[owner]) =
            runtime->advance_phase_barrier(context,
            *(h.ghost_owners_pbarriers_ptrs[owner]));
        } // for
      } // read_phase

      if(write_phase) {
        {
        clog_tag_guard(prolog);
        clog(trace) << "rank " << runtime->find_local_MPI_rank() <<
          " WRITE PHASE PROLOGUE" << std::endl;
        clog(trace) << "rank " << my_color << " wait & arrival barrier " <<
          *(h.pbarrier_as_owner_ptr) << std::endl;
        } // scope

        // Phase WRITE
        launcher.add_wait_barrier(*(h.pbarrier_as_owner_ptr));

        // Phase READ
        launcher.add_arrival_barrier(*(h.pbarrier_as_owner_ptr));
      } // if
      }//end if
    } // handle

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
      if (!h.global && !h.color){
         auto& flecsi_context = context_t::instance();

         auto iitr = flecsi_context.field_info_map().find(
           { h.data_client_hash, h.index_space });

         clog_assert(iitr != flecsi_context.field_info_map().end(),
           "invalid index space");


         for(auto& fitr: iitr->second) {
           const context_t::field_info_t & fi = fitr.second;
           if(fi.fid == h.fid && utils::hash::is_internal(fi.key)){
             return;
           }
         } // for

         bool read_phase = false;
         bool write_phase = false;
         const int my_color = runtime->find_local_MPI_rank();

         read_phase = GHOST_PERMISSIONS != reserved;
         write_phase = (SHARED_PERMISSIONS == wo) || (SHARED_PERMISSIONS == rw);

         if(read_phase) {
          /*
           {
           clog_tag_guard(prolog);
           clog(trace) << "rank " << my_color << " READ PHASE PROLOGUE" <<
             std::endl;

           // As owner
           clog(trace) << "rank " << my_color << " arrives & advances " <<
             *(h.pbarrier_as_owner_ptr) << std::endl;
           } // scope
           */

           {
             // Phase WRITE
             h.pbarrier_as_owner_ptr->arrive(1);

             // Phase WRITE
             *(h.pbarrier_as_owner_ptr) =
                runtime->advance_phase_barrier(context,
                      *(h.pbarrier_as_owner_ptr));
           }

           const size_t _pbp_size = h.ghost_owners_pbarriers_ptrs.size();

           region_info_t* ri;

           auto ritr = region_info_map.find(h.index_space);
       
           if(ritr == region_info_map.end()){
              region_info_map[h.index_space] = region_info_t();
              ri = &region_info_map[h.index_space];
              ri->data_client_hash = h.data_client_hash;
              ri->shared_lrs = h.ghost_owners_lregions;
              ri->ghost_lr = h.ghost_lr;
              ri->color_lr = h.color_region;
              ri->global_to_local_color_map_ptr = 
               h.global_to_local_color_map_ptr;
           }
           else{
              ri = &ritr->second;
           }

           //   if(!utils::hash::is_internal(fi.key)){
           //     rr_shared.add_field(fi.fid);
           //     rr_ghost.add_field(fi.fid);

           ri->fids.push_back(h.fid);

           // As user
           for(size_t owner{0}; owner<_pbp_size; owner++) {
             ri->owners.push_back(owner);
             ri->barriers.push_back((h.ghost_owners_pbarriers_ptrs[owner]));

            /*
             {
             clog_tag_guard(prolog);
             clog(trace) << "rank " << my_color << " WAITS " <<
               *(h.ghost_owners_pbarriers_ptrs[owner]) << std::endl;

             clog(trace) << "rank " << my_color << " arrives & advances " <<
               *(h.ghost_owners_pbarriers_ptrs[owner]) << std::endl;
             } // scope
             */

{
             }

           } // for
         } // read_phase

         if(write_phase) {
           {
           clog_tag_guard(prolog);
           clog(trace) << "rank " << runtime->find_local_MPI_rank() <<
             " WRITE PHASE PROLOGUE" << std::endl;
           clog(trace) << "rank " << my_color << " wait & arrival barrier " <<
             *(h.pbarrier_as_owner_ptr) << std::endl;
           } // scope

           // ndm - is this correct?

           // Phase WRITE
           launcher.add_wait_barrier(*(h.pbarrier_as_owner_ptr));

           // Phase READ
           launcher.add_arrival_barrier(*(h.pbarrier_as_owner_ptr));
         } // if
      }//end if
    }

    //------------------------------------------------------------------------//
    //! FIXME: Need to document.
    //------------------------------------------------------------------------//

    template<
      typename T
    >
    static
    typename std::enable_if_t<!std::is_base_of<data_handle_base_t, T>::value>
    handle(
      T&
    )
    {
    } // handle

    void
    launch_ghost_copy()
    {
      struct index_space_t{
        size_t data_client_hash;
        size_t index_space;
        size_t num_owners;
        size_t owners[8];
        size_t owner_regions[8];
      };

      struct args_t {
        size_t num_index_spaces;
        index_space_t index_spaces[32];
      };

      size_t i = 0;

      args_t args;

      auto& flecsi_context = context_t::instance();

      struct shared_region_info{
        std::set<Legion::FieldID> fids;
        size_t region;
        Legion::LogicalRegion lr;
        Legion::LogicalRegion color_lr;
        Legion::PhaseBarrier* barrier; 
      };

      std::map<Legion::LogicalRegion, shared_region_info> rm;

      size_t num_regions = region_info_map.size();

      for(auto& itr : region_info_map){
        size_t index_space = itr.first;

        region_info_t& ri = itr.second;

        args.index_spaces[i].data_client_hash = ri.data_client_hash;
        args.index_spaces[i].index_space = itr.first;

        args.index_spaces[i].num_owners = ri.owners.size();

        std::memcpy(args.index_spaces[i].owners, &ri.owners[0],
          sizeof(size_t) * args.index_spaces[i].num_owners);

        for (size_t owner = 0; owner < ri.owners.size(); owner++){
          Legion::LogicalRegion ro = ri.shared_lrs[owner];

          auto ritr = rm.find(ro);
          if(ritr == rm.end()){
            shared_region_info si;
            si.region = num_regions + rm.size();
            si.fids.insert(ri.fids.begin(), ri.fids.end());
            si.lr = ro;
            si.color_lr = ri.color_lr;
            si.barrier = ri.barriers[owner];
            rm[ro] = si;
            args.index_spaces[i].owner_regions[owner] = si.region;
          }
          else{
            args.index_spaces[i].owner_regions[owner] = ritr->second.region;
            for(auto fid : ri.fids){
              ritr->second.fids.insert(fid);
            }
          }
        }//end owner 
/*
        Legion::RegionRequirement rr_ghost(ri.ghost_lr,
            WRITE_DISCARD, EXCLUSIVE, ri.color_lr);

        auto ghost_owner_pos_fid = LegionRuntime::HighLevel::FieldID(
            internal_field::ghost_owner_pos);

        rr_ghost.add_field(ghost_owner_pos_fid);

        for(auto fid : ri.fids){
          rr_ghost.add_field(fid);          
        }
        */
        launcher.add_future(Legion::Future::from_value(runtime,
          *(ri.global_to_local_color_map_ptr)));

        ++i;
       }

       std::map<size_t, shared_region_info> nm;
       for(auto& itr : rm){
         nm[itr.second.region] = std::move(itr.second);
       }

       for(auto& itr : nm){
          // Phase READ
        // launcher.add_region_requirement(rr_shared);
        //launcher.add_region_requirement(rr_ghost);
        launcher.add_wait_barrier(*itr.second.barrier);

        // Phase WRITE

        launcher.add_arrival_barrier(*itr.second.barrier);

        /*
         Legion::RegionRequirement rr_shared(itr.second.lr,
           READ_ONLY, EXCLUSIVE, itr.second.color_lr);

         for(auto fid : itr.second.fids){
           rr_shared.add_field(fid);
         }

         launcher.add_region_requirement(rr_shared);
         */
               // Phase WRITE
               *(itr.second.barrier) =
                 runtime->advance_phase_barrier(context,
                 *(itr.second.barrier));
                
       }

      // TODO - circular dependency including internal_task.h
      auto constexpr key = flecsi::utils::const_string_t{
        EXPAND_AND_STRINGIFY(ghost_copy_task)}.hash();

      const auto ghost_copy_tid = flecsi_context.task_id<key>();

      Legion::TaskLauncher launcher2(ghost_copy_tid,
        Legion::TaskArgument(&args, sizeof(args)));

      // Execute the ghost copy task
      runtime->execute_task(context, launcher2);

      /* need?
      for(auto& itr : region_info_map){
        region_info_t& ri = itr.second;
        ri.barrier = runtime->advance_phase_barrier(context,
                 ri.barrier);
      }//for
      */
    }

    struct region_info_t{
      size_t data_client_hash;
      std::vector<Legion::FieldID> fids;
      std::vector<Legion::LogicalRegion> shared_lrs; 
      Legion::LogicalRegion ghost_lr;
      Legion::LogicalRegion color_lr;
      std::vector<Legion::PhaseBarrier*> barriers; 
      const Legion::STL::map<LegionRuntime::Arrays::coord_t,
        LegionRuntime::Arrays::coord_t>* global_to_local_color_map_ptr;
      std::vector<size_t> owners;
    };

    Legion::Runtime* runtime;
    Legion::Context & context;
    Legion::TaskLauncher& launcher;
    std::map<size_t, region_info_t> region_info_map;
    bool first = true; 
  }; // struct task_prolog_t

} // namespace execution 
} // namespace flecsi

#endif // flecsi_execution_legion_task_prolog_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
