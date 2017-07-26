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

#ifndef flecsi_execution_legion_execution_policy_h
#define flecsi_execution_legion_execution_policy_h

//----------------------------------------------------------------------------//
//! @file
//! @date Initial file creation: Nov 15, 2015
//----------------------------------------------------------------------------//

#include <functional>
#include <memory>
#include <type_traits>

#include <cinchlog.h>
#include <legion.h>

#include "flecsi/execution/common/processor.h"
#include "flecsi/execution/context.h"
#include "flecsi/execution/legion/context_policy.h"
#include "flecsi/execution/legion/future.h"
#include "flecsi/execution/legion/runtime_state.h"
#include "flecsi/execution/legion/task_wrapper.h"
#include "flecsi/execution/legion/init_args.h"
#include "flecsi/execution/legion/task_prolog.h"
#include "flecsi/execution/legion/task_epilog.h"
#include "flecsi/utils/const_string.h"
#include "flecsi/data/data_handle.h"

namespace flecsi {
namespace execution {

//----------------------------------------------------------------------------//
// Execution policy.
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
//! The legion_execution_policy_t is the backend runtime execution policy
//! for Legion.
//!
//! @ingroup legion-execution
//----------------------------------------------------------------------------//

struct legion_execution_policy_t
{
  //--------------------------------------------------------------------------//
  //! The future__ type may be used for explicit synchronization of tasks.
  //!
  //! @tparam RETURN The return type of the task.
  //--------------------------------------------------------------------------//

  template<typename RETURN>
  using future__ = legion_future__<RETURN>;

#if 0
  //--------------------------------------------------------------------------//
  //! The task_wrapper__ type FIXME
  //!
  //! @tparam RETURN The return type of the task. FIXME
  //--------------------------------------------------------------------------//

  template<
    typename FUNCTOR_TYPE
  >
  using functor_task_wrapper__ =
    typename flecsi::execution::functor_task_wrapper__<FUNCTOR_TYPE>;
#endif

  //--------------------------------------------------------------------------//
  //! The runtime_state_t type identifies a public type for the high-level
  //! runtime interface to pass state required by the backend.
  //--------------------------------------------------------------------------//

  using runtime_state_t = legion_runtime_state_t;

  //--------------------------------------------------------------------------//
  //! Return the runtime state of the calling FleCSI task.
  //!
  //! @param task The calling task.
  //--------------------------------------------------------------------------//

  static
  runtime_state_t &
  runtime_state(
    void * task
  );

  //--------------------------------------------------------------------------//
  // Task interface.
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  //! This method allows the user to register a pure Legion task with
  //! the runtime. A task id will automatically be generated, and can be
  //! accessed via legion_context_policy_t::task_id using a valid
  //! task hash.
  //!
  //! @tparam KEY    A hash key identifying the task.
  //! @tparam RETURN The return type of the pure Legion task.
  //! @tparam TASK   The function pointer template type of the task.
  //!
  //! @param name The string name for the task. This can be set to any
  //!             valid std::string value.
  //--------------------------------------------------------------------------//

  template<
    size_t KEY,
    typename RETURN,
    RETURN (*TASK)(
      const Legion::Task *,
      const std::vector<Legion::PhysicalRegion> &,
      Legion::Context,
      Legion::Runtime *
    )
  >
  static
  bool
  register_legion_task(
    processor_type_t processor,
    launch_t launch,
    std::string name
  )
  {
    clog(info) << "Registering legion task " << KEY <<
      " " << name << std::endl;

    if(!context_t::instance().register_task(KEY, processor, launch, name,
      pure_task_wrapper__<RETURN, TASK>::registration_callback)) {
      clog(fatal) << "callback registration failed for " << name << std::endl;
    } // if

    return true;
  } // register_legion_task

#if 0
  //--------------------------------------------------------------------------//
  //! Legion backend task registration. For documentation on this
  //! method, please see task__::register_functor_task.
  //--------------------------------------------------------------------------//

  template<
    typename FUNCTOR_TYPE
  >
  static
  bool
  register_functor_task(
    processor_type_t processor,
    launch_t launch,
    std::string name
  )
  {
    using wrapper_t = functor_task_wrapper__<FUNCTOR_TYPE>;

    if(!context_t::instance().register_task(
      typeid(FUNCTOR_TYPE).hash_code(), processor, launch, name,
      wrapper_t::registration_callback)) {
      clog(fatal) << "callback registration failed for " << name << std::endl;
    } // if

    return true;
  } // register_functor_task
#endif

  //--------------------------------------------------------------------------//
  //! Legion backend task registration. For documentation on this
  //! method, please see task__::register_task.
  //--------------------------------------------------------------------------//

  template<
    size_t KEY,
    typename RETURN,
    typename ARG_TUPLE,
    RETURN (*DELEGATE)(ARG_TUPLE)
  >
  static
  bool
  register_task(
    processor_type_t processor,
    launch_t launch,
    std::string name
  )
  {
    using wrapper_t = task_wrapper__<KEY, RETURN, ARG_TUPLE, DELEGATE>;

    if(!context_t::instance().register_task(KEY, processor, launch, name,
      wrapper_t::registration_callback)) {
      clog(fatal) << "callback registration failed for " << name << std::endl;
    } // if

    return true;
  } // register_task

  //--------------------------------------------------------------------------//
  //! Legion backend task execution. For documentation on this
  //! method, please see task__::execute_task.
  //--------------------------------------------------------------------------//

  template<
    size_t KEY,
    typename RETURN,
    typename ARG_TUPLE,
    typename ... ARGS
  >
  static
  decltype(auto)
  execute_task(
    launch_type_t launch,
    size_t parent,
    ARGS && ... args
  )
  {
    using namespace Legion;

    // Make a tuple from the task arguments.
    ARG_TUPLE task_args = std::make_tuple(args ...);

    context_t & context_ = context_t::instance();

    // Get the processor type.
    auto processor_type = context_.processor_type<KEY>();

    // Get the runtime and context from the current task.
#if defined(ENABLE_LEGION_TLS)
    auto legion_runtime = Legion::Runtime::get_runtime();
    auto legion_context = Legion::Runtime::get_context();
#else
    auto legion_runtime = context_.runtime(parent);
    auto legion_context = context_.context(parent);
#endif

    // Handle MPI and Legion invocations separately.
    if(processor_type == processor_type_t::mpi) {
      {
      clog_tag_guard(execution);
      clog(info) << "Executing MPI task: " << KEY << std::endl;
      }

      ArgumentMap arg_map;
      IndexLauncher launcher(
        context_.task_id<KEY>(),
        Legion::Domain::from_rect<1>(context_.all_processes()),
        TaskArgument(&task_args, sizeof(ARG_TUPLE)),
        arg_map
      );

      // Enqueue the MPI task.
      auto future =
        legion_runtime->execute_index_space(legion_context, launcher);
      future.wait_all_results();

      // Handoff to the MPI runtime.
      context_.handoff_to_mpi(legion_context, legion_runtime);

      // Wait for MPI to finish execution (synchronous).
      context_.wait_on_mpi(legion_context, legion_runtime);
      
      // Reset the calling state to false.
      context_.unset_call_mpi(legion_context, legion_runtime);

      return legion_future__<RETURN>(future);
    }
    else {
      // Initialize the arguments to pass through the runtime.
      init_args_t init_args(legion_runtime, legion_context);
      init_args.walk(task_args);

      // Switch on launch type: single or index.
      switch(launch) {

        case launch_type_t::single:
          {
          clog_tag_guard(execution);
          clog(info) << "Executing single task: " << KEY << std::endl;

          // Create a task launcher, passing the task arguments.
          TaskLauncher task_launcher(context_.task_id<KEY>(),
            TaskArgument(&task_args, sizeof(ARG_TUPLE)));

#ifdef MAPPER_COMPACTION
         task_launcher.tag=MAPPER_COMPACTED_STORAGE;
#endif

          for(auto& req : init_args.region_reqs){
            task_launcher.add_region_requirement(req);
          }

          // Enqueue the prolog.
          task_prolog_t
            task_prolog(legion_runtime, legion_context, task_launcher);
          task_prolog.walk(task_args);

          // Enqueue the task.
          clog(trace) << "Execute flecsi/legion task " << KEY << " on rank " <<
              legion_runtime->find_local_MPI_rank() << std::endl;
          auto future = legion_runtime->execute_task(legion_context,
            task_launcher);

          // Enqueue the epilog.
          task_epilog_t
            task_epilog(legion_runtime, legion_context);
          task_epilog.walk(task_args);

          return legion_future__<RETURN>(future);
          } // scope

        case launch_type_t::index:
          {
          clog_tag_guard(execution);
          clog(info) << "Executing index task: " << KEY << std::endl;

          // FIXME:
          // FIXME: This looks incomplete!
          // FIXME:
          //FIXME: get launch domain from partitioning of the data used in
          // the task following launch domeing calculation is temporary:
          LegionRuntime::Arrays::Rect<1> launch_bounds(
            LegionRuntime::Arrays::Point<1>(0),
            LegionRuntime::Arrays::Point<1>(5));
          Domain launch_domain = Domain::from_rect<1>(launch_bounds);

          Legion::ArgumentMap arg_map;
          Legion::IndexLauncher index_launcher(
            context_.task_id<KEY>(), launch_domain,
            TaskArgument(&task_args, sizeof(ARG_TUPLE)), arg_map);

#ifdef MAPPER_COMPACTION
          index_launcher.tag=MAPPER_COMPACTED_STORAGE;
#endif
          // Enqueue the task.
          auto future = legion_runtime->execute_index_space(legion_context,
            index_launcher);

          return legion_future__<RETURN>(future);
          } // scope

        default:
          clog(fatal) << "invalid launch type" << std::endl;
          throw std::runtime_error ("invalid launch type");

      } // switch
    } // if
  } // execute_task

  //--------------------------------------------------------------------------//
  // Function interface.
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  //! Legion backend function registration. For documentation on this
  //! method, please see function__::register_function.
  //--------------------------------------------------------------------------//

template<
    typename RETURN,
    typename ARG_TUPLE,
    RETURN (*FUNCTION)(ARG_TUPLE),
    size_t KEY
  >
  static
  bool
  register_function()
  {
    return context_t::instance().template register_function<
      RETURN, ARG_TUPLE, FUNCTION, KEY>();
  } // register_function
  
  //--------------------------------------------------------------------------//
  //! Legion backend function execution. For documentation on this
  //! method, please see function__::execute_function.
  //--------------------------------------------------------------------------//

  template<
    typename FUNCTION_HANDLE,
    typename ... ARGS
  >
  static
  decltype(auto)
  execute_function(
    FUNCTION_HANDLE & handle,
    ARGS && ... args
  )
  {
    return handle(context_t::instance().function(handle.key()),
      std::forward_as_tuple(args ...));
  } // execute_function

}; // struct legion_execution_policy_t

} // namespace execution 
} // namespace flecsi

#endif // flecsi_execution_legion_execution_policy_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
