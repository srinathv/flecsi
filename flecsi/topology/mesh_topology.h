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

#ifndef flecsi_topology_mesh_topology_h
#define flecsi_topology_mesh_topology_h

/*!
  \file mesh_topology.h
  \authors nickm@lanl.gov, bergen@lanl.gov
  \date Initial file creation: Sep 23, 2015
 */

/*

 Description of major features and terms specific to FleCSI:

 mesh dimension MD - e.g: MD = 2 for a 2d mesh. We currently support 2d and 3d
   meshes.

 topological dimension D - the dimensionality associated with entities,
   e.g: D = 0 is interpreted as a vertex, D = 1 an edge or face for MD = 2,
   D = 2 is a cell for MD = 2

 domain M - a sub-mesh or mesh space that holds entities of various topological
   dimension

 connectivity - a directed connection or adjancy between entities of differing
   topological dimension in the same domain. e.g: D1 -> D2 (edges -> faces)
   for MD = 3. Cell to vertex connectivity is supplied by the user and all
   other connectivies are computed by the topology.

 binding - a type of connectivity that connects entities of potentially
   differing topological dimension across two different domains

 entity - an object associated with a topological dimension, e.g: cell. Each
   entity has an associated integer id.

 mesh topology - the top-level container for domains, entities,
   and connectivities, referred to as the low-level interface

 mesh policy - the top-level class that a specialization creates to
   parameterize the mesh topology to define such things as: mesh dimension,
   number of domains, connectivity and binding pairs of interest, and entity
   classes/types per each domain/topological dimension.

 entity set - contains an iterable set of entities. Support set operations such
   as intersection, union, etc. and functional operations like apply, map,
   reduce, etc. to apply a custom function to the set.

*/

#include <algorithm>
#include <iostream>
#include <array>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <functional>
#include <map>
#include <cstring>
#include <type_traits>
#include <memory>

#include "flecsi/execution/context.h"
#include "flecsi/topology/mesh_storage.h"
#include "flecsi/topology/mesh_types.h"
#include "flecsi/topology/partition.h"
#include "flecsi/utils/common.h"
#include "flecsi/utils/set_intersection.h"
#include "flecsi/utils/static_verify.h"

namespace flecsi {
namespace topology {
namespace verify_mesh {

FLECSI_MEMBER_CHECKER(num_dimensions);
FLECSI_MEMBER_CHECKER(num_domains);
FLECSI_MEMBER_CHECKER(entity_types);
FLECSI_MEMBER_CHECKER(connectivities);
FLECSI_MEMBER_CHECKER(bindings);
FLECSI_MEMBER_CHECKER(create_entity);

} // namespace verify_mesh

/*----------------------------------------------------------------------------*
 * class mesh_topology_t
 *----------------------------------------------------------------------------*/

/*!
  \class mesh_topology_t mesh_topology.h
  \brief mesh_topology_t is parameterized on a class (MT) which gives
    information about its entity types, connectivities and more. the mesh
    topology is responsibly for computing connectivity info between entities
    of different topological dimension, e.g: vertex -> cell,
    cell -> edge, etc. and provides methods for traversing these adjancies.
    It also holds vectors containing the entity instances.
 */
template<
  class MT
>
class mesh_topology_t :
public mesh_topology_base_t<
  mesh_storage_t<MT::num_dimensions, MT::num_domains>
>
{
  // static verification of mesh policy

  static_assert(verify_mesh::has_member_num_dimensions<MT>::value,
                "mesh policy missing num_dimensions size_t");

  static_assert(std::is_convertible<decltype(MT::num_dimensions),
    size_t>::value, "mesh policy num_dimensions must be size_t");



  static_assert(verify_mesh::has_member_num_domains<MT>::value,
                "mesh policy missing num_domains size_t");

  static_assert(std::is_convertible<decltype(MT::num_domains),
    size_t>::value, "mesh policy num_domains must be size_t");



  static_assert(verify_mesh::has_member_entity_types<MT>::value,
                "mesh policy missing entity_types tuple");

  static_assert(utils::is_tuple<typename MT::entity_types>::value,
                "mesh policy entity_types is not a tuple");



  static_assert(verify_mesh::has_member_connectivities<MT>::value,
                "mesh policy missing connectivities tuple");

  static_assert(utils::is_tuple<typename MT::connectivities>::value,
                "mesh policy connectivities is not a tuple");



  static_assert(verify_mesh::has_member_bindings<MT>::value,
                "mesh policy missing bindings tuple");

  static_assert(utils::is_tuple<typename MT::bindings>::value,
                "mesh policy bindings is not a tuple");



  static_assert(verify_mesh::has_member_create_entity<MT>::value,
                "mesh policy missing create_entity()");

public:

  using storage_t = mesh_storage_t<MT::num_dimensions, MT::num_domains>;

  using base_t = 
    mesh_topology_base_t<storage_t>;

  using id_t = utils::id_t;

  using offset_t = utils::offset_t;
  
  // used to find the entity type of topological dimension D and domain M
  template<size_t D, size_t M = 0>
  using entity_type = typename find_entity_<MT, D, M>::type;

  //--------------------------------------------------------------------------//
  // This type definition is needed so that data client handles can be
  // specialized for particular data client types, e.g., mesh topologies vs.
  // tree topologies. It is also useful for detecting illegal usage, such as
  // when a user adds data members.
  //--------------------------------------------------------------------------//

  using type_identifier_t = mesh_topology_t;

  // Don't allow the mesh to be copied or copy constructed

  mesh_topology_t & operator=(const mesh_topology_t &) = delete;

  // Allow move operations
  mesh_topology_t(mesh_topology_t && o) = default;


  //! override default move assignement
  mesh_topology_t & operator=(mesh_topology_t && o) = default;

  //! Constructor
  mesh_topology_t(storage_t * ms = nullptr)
    : base_t(ms)
  {
    if(ms != nullptr) {
      initialize_storage();
    } // if
  } // mesh_topology_t()

  mesh_topology_t(const mesh_topology_t & m)
    : base_t(m.ms_) {}

  // The mesh retains ownership of the entities and deletes them
  // upon mesh destruction
  virtual
  ~mesh_topology_t() {}

  void
  initialize_storage() {

    for (size_t from_domain = 0; from_domain < MT::num_domains; ++from_domain) {
      for (size_t to_domain = 0; to_domain < MT::num_domains; ++to_domain) {
        base_t::ms_->topology[from_domain][to_domain].
          init_(from_domain, to_domain);
      } // for
    } // for

    for (size_t to_domain = 0; to_domain < MT::num_domains; ++to_domain) {
      for (size_t to_dim = 0; to_dim <= MT::num_dimensions; ++to_dim) {
        auto& master = 
          base_t::ms_->index_spaces[to_domain][to_dim];

        for (size_t from_domain = 0; from_domain < MT::num_domains;
             ++from_domain) {
          for (size_t from_dim = 0; from_dim <= MT::num_dimensions;
               ++from_dim) {
            get_connectivity_(from_domain, to_domain, from_dim, to_dim).
              get_index_space().set_master(master);
          } // for
        } // for
      } // for
    } // for
  } // intialize_storage

  // A mesh is constructed by creating cells and vertices and associating
  // vertices with cells as in this method.
  template<
    size_t M,
    class C,
    typename V
  >
  void
  init_cell(
    C * cell,
    V && verts
 )
  {
    init_cell_<M>(cell, std::forward<V>(verts));
  } // init_cell

  template<
    size_t M,
    class C,
    typename V
  >
  void
  init_cell(
    C * cell,
    std::initializer_list<V *> verts
 )
  {
    init_cell_<M>(cell, verts);
  } // init_cell

  // Initialize an entities connectivity with a subset of another
  template<
    size_t M,
    size_t D1,
    size_t D2,
    class E1,
    class E2
  >
  void
  init_entity(
    E1 * super,
    E2 && subs
 )
  {
    init_entity_<M,D1,D2>(super, std::forward<E2>(subs));
  } // init_entity

  template<
    size_t M,
    size_t D1,
    size_t D2,
    class E1,
    class E2
  >
  void
  init_entity(
    E1 * super,
    std::initializer_list<E2*> subs
 )
  {
    init_entity_<M,D1,D2>(super, subs);
  } // init_entity

  // Virtual method of num_entities_()
  size_t
  num_entities(
    size_t dim,
    size_t domain=0
 ) const override
  {
    return num_entities_(dim, domain);
  } // num_entities

  /*!
    The init method builds entities as edges/faces and computes adjacencies
    and bindings.
   */
  template<
    size_t M = 0
  >
  void init()
  {
    // Compute mesh connectivity
    using TP = typename MT::connectivities;
    compute_connectivity_<M, std::tuple_size<TP>::value, TP>::compute(*this);

    using BT = typename MT::bindings;
    compute_bindings_<M, std::tuple_size<BT>::value, BT>::compute(*this);
  } // init

  /*!
    Similar to init(), but only compute bindings. This method should be called
    when a domain is sparse, i.e: missing certain entity types such as cells
    and it is not possible to compute connectivities.
   */
  template<
    size_t M = 0
  >
  void init_bindings()
  {
    using BT = typename MT::bindings;
    compute_bindings_<M, std::tuple_size<BT>::value, BT>::compute(*this);
  } // init

  /*!
   Return the number of entities contained in specified topological dimension
   and domain.
   */
  template<
    size_t D,
    size_t M = 0
    >
  decltype(auto)
  num_entities() const
  {
    return base_t::ms_->index_spaces[M][D].size();
  } // num_entities

  /*!
   Return the number of entities contained in specified topological dimension
   and domain.
   */
  template<
    size_t D,
    size_t M = 0
    >
  decltype(auto)
  num_entities(partition_t partition) const
  {
    return base_t::ms_->partition_index_spaces[partition][M][D].size();
  } // num_entities

  /*!
   Get the connectivity of the specified from/to domain and from/to topological
   dimensions.
   */
  const connectivity_t &
  get_connectivity(
    size_t from_domain,
    size_t to_domain,
    size_t from_dim,
    size_t to_dim) const override
  {
    return get_connectivity_(from_domain, to_domain, from_dim, to_dim);
  } // get_connectivity

  /*!
   Get the connectivity of the specified from/to domain and from/to topological
   dimensions.
   */
  connectivity_t &
  get_connectivity(
    size_t from_domain,
    size_t to_domain,
    size_t from_dim,
    size_t to_dim) override
  {
    return get_connectivity_(from_domain, to_domain, from_dim, to_dim);
  } // get_connectivity

  /*!
   Get the connectivity of the specified domain and from/to topological
   dimensions.
   */
  const connectivity_t &
  get_connectivity(
    size_t domain,
    size_t from_dim,
    size_t to_dim) const override
  {
    return get_connectivity_(domain, domain, from_dim, to_dim);
  } // get_connectivity

  /*!
   Get the connectivity of the specified domain and from/to topological
   dimensions.
   */
  connectivity_t &
  get_connectivity(
    size_t domain,
    size_t from_dim,
    size_t to_dim) override
  {
    return get_connectivity_(domain, domain, from_dim, to_dim);
  } // get_connectivity

  size_t
  topological_dimension() const override
  {
    return MT::num_dimensions;
  }

  template<
    size_t M = 0
  >
  const auto &
  get_index_space_(
    size_t dim
  ) const
  {
    return base_t::ms_->index_spaces[M][dim];
  } // get_entities_

  template<
    size_t M = 0
  >
  auto &
  get_index_space_(
    size_t dim
  )
  {
    return base_t::ms_->index_spaces[M][dim];
  } // get_entities_

  template<
    size_t M = 0
  >
  const auto &
  get_index_space_(
    size_t dim,
    partition_t partition
  ) const
  {
    return base_t::ms_->partition_index_spaces[partition][M][dim];
  } // get_entities_

  template<
    size_t M = 0
  >
  auto &
  get_index_space_(
    size_t dim,
    partition_t partition
  )
  {
    return base_t::ms_->partition_index_spaces[partition][M][dim];
  } // get_entities_

  /*!
    Get an entity in domain M of topological dimension D with specified id.
  */
  template<
    size_t D,
    size_t M = 0
  >
  auto
  get_entity(
    id_t global_id
  ) const
  {
    using etype = entity_type<D, M>;
    return static_cast<etype *>(
      base_t::ms_->index_spaces[M][D][global_id.entity()]);
  } // get_entity

  /*!
    Get an entity in domain M of topological dimension D with specified id.
  */
  template<
    size_t M = 0
  >
  auto
  get_entity(
    size_t dim,
    id_t global_id
  )
  {
    return base_t::ms_->index_spaces[M][dim][global_id.entity()];
  } // get_entity

  /*!
    Get an entity in domain M of topological dimension D with specified id.
  */
  template<
    size_t D,
    size_t M = 0
  >
  auto
  get_entity(
    id_t global_id,
    partition_t partition
  ) const
  {
    using etype = entity_type<D, M>;
    return static_cast<etype *>(
      base_t::ms_->partition_index_spaces[partition][M][D][global_id.entity()]);
  } // get_entity

  /*!
    Get an entity in domain M of topological dimension D with specified id.
  */
  template<
    size_t M = 0
  >
  auto
  get_entity(
    size_t dim,
    id_t global_id,
    partition_t partition
  )
  {
    return base_t::ms_->partition_index_spaces[partition]
      [M][dim][global_id.entity()];
  } // get_entity

  /*!
    Get the entities of topological dimension D connected to another entity
    by specified connectivity from domain FM and to domain TM.
  */
  template<
    size_t D,
    size_t FM,
    size_t TM = FM,
    class E
  >
  const auto
  entities(
    const E * e
  ) const
  {

    const connectivity_t & c = get_connectivity(FM, TM, E::dimension, D);
    assert(!c.empty() && "empty connectivity");

    using etype = entity_type<D, TM>;
    using dtype = domain_entity<TM, etype>;

    return c.get_index_space().slice<dtype>(c.range(e->template id<FM>()));
  } // entities

  /*!
    Get the entities of topological dimension D connected to another entity
    by specified connectivity from domain FM and to domain TM.
  */
  template<
    size_t D,
    size_t FM,
    size_t TM = FM,
    class E
  >
  auto
  entities(
    E * e 
  )
  {
    connectivity_t & c = get_connectivity(FM, TM, E::dimension, D);
    assert(!c.empty() && "empty connectivity");

    using etype = entity_type<D, TM>;
    using dtype = domain_entity<TM, etype>;

    return c.get_index_space().slice<dtype>(c.range(e->template id<FM>()));
  } // entities

  /*!
    Get the entities of topological dimension D connected to another entity
    by specified connectivity from domain FM and to domain TM.
  */
  template<
    size_t D,
    size_t FM = 0,
    size_t TM = FM,
    class E
  >
  decltype(auto)
  entities(
    domain_entity<FM, E> & e
  ) const
  {
    return entities<D, FM, TM>(e.entity());
  } // entities

  /*!
    Get the entities of topological dimension D connected to another entity
    by specified connectivity from domain FM and to domain TM.
  */
  template<
    size_t D,
    size_t FM = 0,
    size_t TM = FM,
    class E
  >
  decltype(auto)
  entities(
    domain_entity<FM, E> & e
  )
  {
    return entities<D, FM, TM>(e.entity());
  } // entities

  /*!
    Get the top-level entities of topological dimension D of the specified
    domain M. e.g: cells of the mesh.
  */
  template<
    size_t D,
    size_t M = 0
  >
  auto
  entities() const
  {
    using etype = entity_type<D, M>;
    using dtype = domain_entity<M, etype>;
    return base_t::ms_->index_spaces[M][D].template slice<dtype>();
  } // entities

  /*!
    Get the top-level entities of topological dimension D of the specified
    domain M. e.g: cells of the mesh.
  */
  template<
    size_t D,
    size_t M = 0
  >
  auto
  entities(partition_t partition) const
  {
    using etype = entity_type<D, M>;
    using dtype = domain_entity<M, etype>;
    return base_t::ms_->partition_index_spaces[partition]
      [M][D].template slice<dtype>();
  } // entities

  /*!
    Get the top-level entity id's of topological dimension D of the specified
    domain M. e.g: cells of the mesh.
  */
  template<
    size_t D,
    size_t M = 0
  >
  auto
  entity_ids() const
  {
    return base_t::ms_->index_spaces[M][D].ids();
  } // entity_ids

  /*!
    Get the top-level entity id's of topological dimension D of the specified
    domain M. e.g: cells of the mesh.
  */
  template<
    size_t D,
    size_t M = 0
  >
  auto
  entity_ids(partition_t partition) const
  {
    return base_t::ms_->partition_index_spaces[partition][M][D].ids();
  } // entity_ids

  /*!
    Get the entity id's of topological dimension D connected to another entity
    by specified connectivity from domain FM and to domain TM.
  */
  template<
    size_t D,
    size_t FM = 0,
    size_t TM = FM,
    class E
  >
  decltype(auto)
  entity_ids(
    domain_entity<FM, E> & e
  )
  {
    return entity_ids<D, FM, TM>(e.entity());
  } // entities

  /*!
    Get the entity id's of topological dimension D connected to another entity
    by specified connectivity from domain FM and to domain TM.
  */
  template<
    size_t D,
    size_t FM = 0,
    size_t TM = FM,
    class E
  >
  auto
  entity_ids(
    const E * e
  ) const
  {
    const connectivity_t & c = get_connectivity(FM, TM, E::dimension, D);
    assert(!c.empty() && "empty connectivity");
    return c.get_index_space().ids(c.range(e->template id<FM>()));
  } // entities

  /*!
    Get the entities of topological dimension D connected to another entity
    by specified connectivity from domain FM and to domain TM.
  */
  template<
    size_t D,
    size_t FM,
    size_t TM = FM,
    class E
  >
  void
  reverse_entities(
    E * e
  )
  {
    auto & c = get_connectivity(FM, TM, E::dimension, D);
    assert(!c.empty() && "empty connectivity");
    c.reverse_entities(e->template id<FM>());
  } // entities

  /*!
    Get the entities of topological dimension D connected to another entity
    by specified connectivity from domain FM and to domain TM.
  */
  template<
    size_t D,
    size_t FM = 0,
    size_t TM = FM,
    class E
  >
  void
  reverse_entities(
    domain_entity<FM, E> & e
  )
  {
    return reverse_entities<D, FM, TM>(e.entity());
  } // entities


  /*!
    Get the entities of topological dimension D connected to another entity
    by specified connectivity from domain FM and to domain TM.
  */
  template<
    size_t D,
    size_t FM,
    size_t TM = FM,
    class E,
    class U
  >
  void
  reorder_entities(
    E * e,
    U && order
  )
  {
    auto & c = get_connectivity(FM, TM, E::dimension, D);
    assert(!c.empty() && "empty connectivity");
    c.reorder_entities(e->template id<FM>(), std::forward<U>(order));
  } // entities

  /*!
    Get the entities of topological dimension D connected to another entity
    by specified connectivity from domain FM and to domain TM.
  */
  template<
    size_t D,
    size_t FM = 0,
    size_t TM = FM,
    class E,
    class U
  >
  void
  reverse_entities(
    domain_entity<FM, E> & e,
    U && order
  )
  {
    return reorder_entities<D, FM, TM>(e.entity(), std::forward<U>(order));
  } // entities

  template<
    typename I
  >
  void
  compute_graph_partition(
    size_t domain,
    size_t dim,
    const std::vector<I>& partition_sizes,
    std::vector<mesh_graph_partition<I>>& partitions
  ){

    using int_t = I;

    partitions.reserve(partition_sizes.size());

    int_t total_size = 0;
    for(auto pi : partition_sizes){
      total_size += pi;
    }

    size_t n = num_entities_(dim, domain);
    size_t pn = n / total_size;

    size_t to_dim;

    if (dim == 0) {
      // vertex -> vertex via shared edge.
      to_dim = 1;
    } else {
      // edge -> edge via shared vertex, cell -> cell via shared edge/face etc.
      to_dim = dim - 1;
    }

    const connectivity_t& c1 = get_connectivity(domain, dim, to_dim);
    assert(!c1.empty() && "empty connectivity c1");
    const auto& o1 = c1.offsets();

    const connectivity_t& c2 = get_connectivity(domain, to_dim, dim);
    assert(!c2.empty() && "empty connectivity c2");
    const auto& o2 = c2.offsets();

    mesh_graph_partition<int_t> cp;
    cp.offset.reserve(pn);

    size_t offset = 0;
    size_t pi = 0;

    std::vector<int_t> partition;
    partition.push_back(0);

    for(size_t from_id = 0; from_id < n; ++from_id){
      auto to_ids = c1.get_index_space().ids(o1.range(from_id));
      cp.offset.push_back(offset);

      for(auto to_id : to_ids){
        auto ret_ids =
          c2.get_index_space().ids(o2.range(to_id.entity()));

        for(auto ret_id : ret_ids){
          if(ret_id.entity() != from_id){
            cp.index.push_back(ret_id.local_id());
            ++offset;
          }
        }
      }

      size_t m = cp.offset.size();

      if(m >= pn * partition_sizes[pi]){
        partitions.emplace_back(std::move(cp));
        partition.push_back(m + partition.back());
        offset = 0;
        ++pi;
      }
    }

    for(auto& pi : partitions){
      pi.partition = partition;
    }
  }

  /*!
    Debug method to dump the connectivity of the mesh over all domains and
    topological dimensions.
  */

  std::ostream &
  dump(
    std::ostream & stream
  )
  {
    for (size_t from_domain = 0; from_domain < MT::num_domains; ++from_domain) {
      stream << "=================== from domain: " << from_domain
                << std::endl;
      for (size_t to_domain = 0; to_domain < MT::num_domains; ++to_domain) {
        stream << "========== to domain: " << to_domain << std::endl;
        base_t::ms_->topology[from_domain][to_domain].dump(stream);
      }
    }
    return stream;
  } // dump

  void dump()
  {
    dump(std::cout);
  } // dump

  template<
    typename A
  >
  void
  save(
    A & archive
  ) const {
    size_t size;
    char* data = serialize_(size);
    archive.saveBinary(&size, sizeof(size));

    archive.saveBinary(data, size);
    delete [] data;
  } // save

  template<
    typename A>
  void
  load(
    A & archive
  ){
    size_t size;
    archive.loadBinary(&size, sizeof(size));

    char* data = new char [size];
    archive.loadBinary(data, size);
    unserialize_(data);
    delete [] data;
  } // load

  char*
  serialize_(
    uint64_t& size
  ) const
  {
    const size_t alloc_size = 1048576;
    size = alloc_size;

    char* buf = new char [alloc_size];
    uint64_t pos = 0;

    uint32_t num_domains = MT::num_domains;
    std::memcpy(buf + pos, &num_domains, sizeof(num_domains));
    pos += sizeof(num_domains);

    uint32_t num_dimensions = MT::num_dimensions;
    std::memcpy(buf + pos, &num_dimensions, sizeof(num_dimensions));
    pos += sizeof(num_dimensions);

    for(size_t domain = 0; domain < MT::num_domains; ++domain){
      for(size_t dimension = 0; dimension <= MT::num_dimensions; ++dimension){
        uint64_t num_entities = 
          base_t::ms_->entities[domain][dimension].size();
        std::memcpy(buf + pos, &num_entities, sizeof(num_entities));
        pos += sizeof(num_entities);
      }
    }

    for(size_t from_domain = 0; from_domain < MT::num_domains; ++from_domain){
      for(size_t to_domain = 0; to_domain < MT::num_domains; ++to_domain){

        auto& dc = base_t::ms_->topology[from_domain][to_domain];

        for(size_t from_dim = 0; from_dim <= MT::num_dimensions; ++from_dim){
          for(size_t to_dim = 0; to_dim <= MT::num_dimensions; ++to_dim){
            const connectivity_t& c = dc.get(from_dim, to_dim);

            auto& tv = c.to_id_storage();
            uint64_t num_to = tv.size();
            std::memcpy(buf + pos, &num_to, sizeof(num_to));
            pos += sizeof(num_to);

            size_t bytes = num_to * sizeof(id_vector_t::value_type);

            if(size - pos < bytes){
              size += bytes + alloc_size;
              buf = (char*)std::realloc(buf, size);
            }

            std::memcpy(buf + pos, tv.data(), bytes);
            pos += bytes;

            uint64_t num_offsets = c.offsets().size();
            std::memcpy(buf + pos, &num_offsets, sizeof(num_offsets));
            pos += sizeof(num_offsets);

            bytes = num_offsets * sizeof(offset_t);

            if(size - pos < bytes){
              size += bytes + alloc_size;
              buf = (char*)std::realloc(buf, size);
            }

            std::memcpy(buf + pos, c.offsets().storage().buffer(), bytes);
            pos += bytes;
          }
        }
      }
    }

    size = pos;

    return buf;
  }

  void
  unserialize_(
    char* buf
  )
  {
    uint64_t pos = 0;

    uint32_t num_domains;
    std::memcpy(&num_domains, buf + pos, sizeof(num_domains));
    pos += sizeof(num_domains);
    assert(num_domains == MT::num_domains && "domain size mismatch");

    uint32_t num_dimensions;
    std::memcpy(&num_dimensions, buf + pos, sizeof(num_dimensions));
    pos += sizeof(num_dimensions);
    assert(num_dimensions == MT::num_dimensions && "dimension size mismatch");

    unserialize_domains_<storage_t,
      MT, MT::num_domains,
      MT::num_dimensions, 0>::unserialize(*this, buf, pos);

    for(size_t from_domain = 0; from_domain < MT::num_domains; ++from_domain){
      for(size_t to_domain = 0; to_domain < MT::num_domains; ++to_domain){

        auto& dc = base_t::ms_->topology[from_domain][to_domain];

        for(size_t from_dim = 0; from_dim <= MT::num_dimensions; ++from_dim){
          for(size_t to_dim = 0; to_dim <= MT::num_dimensions; ++to_dim){
            connectivity_t& c = dc.get(from_dim, to_dim);

            auto& tv = c.to_id_storage();
            uint64_t num_to;
            std::memcpy(&num_to, buf + pos, sizeof(num_to));
            pos += sizeof(num_to);
            auto ta = (id_vector_t::value_type*)(buf + pos);
            tv.resize(num_to);
            tv.assign(ta, ta + num_to);
            pos += num_to * sizeof(id_vector_t::value_type);

            auto offsets_buf = c.offsets().storage().buffer();
            uint64_t num_offsets;
            std::memcpy(&num_offsets, buf + pos, sizeof(num_offsets));
            pos += sizeof(num_offsets);
            std::memcpy(offsets_buf, buf + pos, 
              num_offsets * sizeof(offset_t));
            pos += num_offsets * sizeof(offset_t);
          }
        }
      }
    }
  }

  void
  append_to_index_space_(
    size_t domain,
    size_t dim,
    std::vector<mesh_entity_base_*>& ents,
    std::vector<id_t>& ids) override
  {
    auto& is =  base_t::ms_->index_spaces[domain][dim];
    is.append_(ents, ids);
  }

private:

  template<size_t DM, size_t I, class TS>
  friend class compute_connectivity_;

  template<size_t DM, size_t I, class TS>
  friend class compute_bindings_;

  template<
    size_t M,
    typename V>
  void
  init_cell_(
    entity_type<MT::num_dimensions, M> * cell,
    V && verts
  )
  {
    auto & c = get_connectivity_(M, MT::num_dimensions, 0);

    assert(cell->template id<M>() == c.from_size() && "id mismatch");

    for (entity_type<0, M> * v : std::forward<V>(verts)) {
      c.push(v->template global_id<M>());
    } // for

    c.add_count(verts.size());
  } // init_cell

  template<
    size_t M,
    size_t D1,
    size_t D2,
    class E2
  >
  void
  init_entity_(
    entity_type<D1, M> * super,
    E2 && subs
  )
  {
    auto & c = get_connectivity_(M, D1, D2);

    assert(super->template id<M>() == c.from_size() && "id mismatch");

    for (auto e : std::forward<E2>(subs)) {
      c.push(e->template global_id<M>());
    } // for

    c.add_count(subs.size());
  } // init_entity

  // Get the number of entities in a given domain and topological dimension
  size_t
  num_entities_(
    size_t dim,
    size_t domain=0
  ) const
  {
    return base_t::ms_->index_spaces[domain][dim].size();
  } // num_entities_

  // Get the number of entities in a given domain and topological dimension
  size_t
  num_entities_(
    size_t dim,
    size_t domain,
    partition_t partition
  ) const
  {
    return base_t::ms_->partition_index_spaces[partition][domain][dim].size();
  } // num_entities_

  /*!
    Build connectivity informaiton and add entities to the mesh for the
    given dimension.

    \remark this is the general one that gets instantiated even though
            it may never get called
  */
  template<
    size_t Domain,
    size_t DimensionToBuild,
    size_t UsingDimension>
    typename std::enable_if< (UsingDimension <= 1 ||
      UsingDimension > MT::num_dimensions) >::type
  build_connectivity()
  {
    assert(false && "shouldn't be in here");
  }

  /*!
    Build connectivity informaiton and add entities to the mesh for the
    given dimension.

    \remark This one is enable_if'd so it never gets instantiated in certain
            cases, otherwise we would need create_entities in wedges
            and vertices
   */
  template<
    size_t Domain,
    size_t DimensionToBuild,
    size_t UsingDimension>
    typename std::enable_if< (UsingDimension > 1 &&
      UsingDimension <= MT::num_dimensions) >::type
  build_connectivity()
  {
    // std::cerr << "build: " << DimensionToBuild
    // << " using " << UsingDimension << std::endl;

    // Sanity check
    static_assert(
      DimensionToBuild <= MT::num_dimensions,
      "DimensionToBuild must be <= total number of dimensions"
   );
    static_assert(
      UsingDimension <= MT::num_dimensions,
      "UsingDimension must be <= total number of dimensions"
   );
    static_assert(
      Domain < MT::num_domains,
      "Domain must be < total number of domains"
   );

    // Reference to storage from cells to the entity (to be created here).
    connectivity_t & cell_to_entity =
      get_connectivity_(Domain, UsingDimension, DimensionToBuild);

    // Storage for entity-to-vertex connectivity information.
    connection_vector_t entity_vertex_conn;

    // Helper variables
    size_t max_cell_entity_conns = 1;

    // keep track of the local ids, since they may be added out of order
    std::vector<size_t> entity_ids;

    domain_connectivity<MT::num_dimensions> & dc = 
      base_t::ms_->topology[Domain][Domain];

    // Get connectivity for cells to vertices.
    connectivity_t & cell_to_vertex = dc.template get<UsingDimension>(0);
    assert(!cell_to_vertex.empty());

    const size_t _num_cells = num_entities<UsingDimension, Domain>();

    // Storage for cell-to-entity connectivity information.
    connection_vector_t cell_entity_conn(_num_cells);

    // This map is primarily used to make sure that entities are not
    // created multiple times, i.e., that they are unique.  The
    // emplace method of the map is used to only define a new entity
    // if it does not already exist in the map.
    id_vector_map_t entity_vertices_map;

    // This buffer should be large enough to hold all entities
    // vertices that potentially need to be created
    std::array<id_t, 4096> entity_vertices;

    using cell_type = entity_type<UsingDimension, Domain>;
    using entity_type = entity_type<DimensionToBuild, Domain>;

    auto& is = 
      base_t::ms_->index_spaces[Domain][DimensionToBuild].
      template cast<domain_entity<Domain, entity_type>>();
    
    auto& cis = 
      base_t::ms_->index_spaces[Domain][UsingDimension].
      template cast<domain_entity<Domain, cell_type>>();

    // Lookup the index space for the entity type being created.
    constexpr size_t cell_index_space =
      find_index_space_from_dimension__<
        std::tuple_size<typename MT::entity_types>::value,
        typename MT::entity_types,
        UsingDimension,
        Domain
      >::find();

    // Lookup the index space for the vertices from the mesh
    // specialization.
    constexpr size_t vertex_index_space =
      find_index_space_from_dimension__<
        std::tuple_size<typename MT::entity_types>::value,
        typename MT::entity_types,
        0,
        Domain
      >::find();

    // Lookup the index space for the entity type being created.
    constexpr size_t entity_index_space =
      find_index_space_from_dimension__<
        std::tuple_size<typename MT::entity_types>::value,
        typename MT::entity_types,
        DimensionToBuild,
        Domain
      >::find();

    // get the global to local index space map
    auto & context_ = flecsi::execution::context_t::instance();
    size_t color = context_.color();
    auto& gis_to_cis = context_.reverse_index_map(cell_index_space);

    // Get the reverse map of the intermediate ids. This map takes
    // vertices defining an entity to the entity id in MIS.
    auto & reverse_intermediate_map =
      context_.reverse_intermediate_map(DimensionToBuild, Domain);
    auto has_intermediate_map = !reverse_intermediate_map.empty();

    // Get the index map for the entity.
    auto & entity_index_map = context_.reverse_index_map(entity_index_space);

    // Get the map of the vertex ids. This map takes
    // local compacted vertex ids to mesh index space ids.
    // CIS -> MIS.
    auto & vertex_map =
      context_.index_map(vertex_index_space);

    // a counter for added entityes
    size_t entity_counter{0};

    for(auto& citr : gis_to_cis){
      size_t c = citr.second;

      // Get the cell object
      auto cell = static_cast<cell_type*>(cis[c]);
      id_t cell_id = cell->template global_id<Domain>();

      // Get storage reference.
      id_vector_t & conns = cell_entity_conn[c];

      // Try to optimize storage.
      conns.reserve(max_cell_entity_conns);

      // This call allows the users specialization to create
      // whatever entities are needed to complete the mesh.
      //
      // p.first:   The number of entities per cell.
      // p.second:  A std::vector of id_t containing the ids of the
      //            vertices that define the entity.
      auto sv = cell->template create_entities(cell_id,
        DimensionToBuild, dc, entity_vertices.data());

      size_t n = sv.size();

      // iterate over the newly-defined entities
      for (size_t i = 0; i < n; ++i) {
        size_t m = sv[i];

        // Get the vertices that define this entity by getting
        // a pointer to the vector-of-vector data and then constructing
        // a vector of ids for only this entity.
        id_t * a = &entity_vertices[i * m];
        id_vector_t ev(a, a + m);

        // Sort the ids for the current entity so that they are
        // monotonically increasing. This ensures that entities are
        // created uniquely (using emplace_back below) because the ids
        // will always occur in the same order for the same entity.
        std::sort(ev.begin(), ev.end());

        //
        // The following set of steps use the vertices that define
        // the entity to be created to lookup the id so
        // that the topology creates it at the correct offset.
        // This requires:
        //
        // 1) lookup the MIS vertex ids
        // 2) create a vector of the MIS vertex ids
        // 3) lookup the MIS id of the entity
        // 4) lookup the CIS id of the entity
        //
        // The CIS id of the entity is passed to the create_entity
        // method. The specialization developer must pass this
        // information to 'make' so that the coloring id of the
        // entity is consitent with the id/offset of the entity
        // created by the topology.
        //

        size_t entity_id;
        if ( has_intermediate_map ) {

          std::vector<size_t> vertices_mis;
          vertices_mis.reserve(m);

          // Push the MIS vertex ids onto a vector to search for the
          // associated entity.
          for(id_t * aptr{a}; aptr<(a+m); ++aptr) {
            vertices_mis.push_back(vertex_map[aptr->entity()]);
          } // for

          // Lookup the MIS id of the entity.
          const auto entity_id_mis = reverse_intermediate_map.at(vertices_mis);

          // Lookup the CIS id of the entity.
          entity_id = entity_index_map.at(entity_id_mis);

        }
        else {

          entity_id = entity_counter;

        } // intermediate_map

        id_t id = id_t::make<DimensionToBuild, Domain>(entity_id, color);

        // Emplace the sorted vertices into the entity map
        auto itr = entity_vertices_map.emplace(
            std::move(ev), id_t::make<DimensionToBuild, Domain>(
          entity_id, cell_id.partition()));

        // Add this id to the cell to entity connections
        conns.push_back(itr.first->second);
      
        // If the insertion took place
        if (itr.second) {

          // what does this do?
          id_vector_t ev2 = id_vector_t(a, a + m);
          entity_vertex_conn.emplace_back(std::move(ev2));
          entity_ids.emplace_back( entity_id );

          max_cell_entity_conns =
            std::max(max_cell_entity_conns, conns.size());

          auto ent =
            MT::template create_entity<Domain, DimensionToBuild>(this, m, id);

          ++entity_counter;

        } // if
      } // for
    } // for
  
    // sort the entity connectivity. Entities may have been created out of
    // order.  Sort them using the list of entity ids we kept track of
    if ( has_intermediate_map )
      utils::reorder_destructive(
        entity_ids.begin(),
        entity_ids.end(),
        entity_vertex_conn.begin()
      );

    // Set the connectivity information from the created entities to
    // the vertices.
    connectivity_t & entity_to_vertex = dc.template get<DimensionToBuild>(0);
    entity_to_vertex.init(entity_vertex_conn);
    cell_to_entity.init(cell_entity_conn);
  } // build_connectivity

  /*!
     used internally to compute connectivity information for
     topological dimension
       FD -> TD where FD < TD
   */
  template<
    size_t FM,
    size_t TM,
    size_t FD,
    size_t TD
  >
  void
  transpose()
  {
    //std::cerr << "transpose: " << FD << " -> " << TD << std::endl;

    // The connectivity we will be populating
    auto & out_conn = get_connectivity_(FM, TM, FD, TD);
    if (!out_conn.empty()) {
      return;
    } // if
    
    // get the list of "to" entities
    const auto & to_entities = entities<TD, TM>();

    index_vector_t pos(num_entities_(FD, FM), 0);

    // Count how many connectivities go into each slot
    for (auto to_entity : to_entities) {
      for (id_t from_id : entity_ids<FD, TM, FM>(to_entity)) {
        ++pos[from_id.entity()];
      }
    }

    out_conn.resize(pos);

    std::fill(pos.begin(), pos.end(), 0);
    

    // now do the actual transpose
    for (auto to_entity : to_entities) {
      for (auto from_id : entity_ids<FD, TM, FM>(to_entity)) {
        auto from_lid = from_id.entity();
        out_conn.set(from_lid, to_entity->template global_id<TM>(),
            pos[from_lid]++);
      }
    }

    // now we need to sort the connecvtivity arrays:
    // .. we have to make sure the order of connectivity information apears in
    //    in order of global id

    // we need the context to get the global-to-local mapping
    const auto & context_ = flecsi::execution::context_t::instance();

    // find the from index space and get the mapping from global to local
    constexpr size_t to_index_space =
      find_index_space_from_dimension__<
        std::tuple_size<typename MT::entity_types>::value,
        typename MT::entity_types,
        TD,
        TM
      >::find();

    const auto& to__cis_to_gis = context_.index_map(to_index_space);

    // do the final sort of the connectivity arrays
    for (auto from_id : entity_ids<FD, TM>()) {
      // get the connectivity array
      size_t count;
      auto conn = out_conn.get_entities( from_id.entity(), count );
      // pack it into a list of id and global id pairs
      std::vector< std::pair<size_t, id_t> > gids( count );
      std::transform(
        conn, conn+count, gids.begin(),
        [&](auto id) {
          return std::make_pair( to__cis_to_gis.at(id.entity()), id );
        }
      );
      // sort via global id 
      std::sort(
        gids.begin(),
        gids.end(),
        []( auto a, auto b ) {
          return a.first < b.first;
        }
      );
      // upack the results
      std::transform(
        gids.begin(), gids.end(), conn,
        [](auto id_pair) {
          return id_pair.second;
        }
      );
    }
  } // transpose

  /*!
     Used internally to compute connectivity information for
     topological dimension
       FD -> TD using FD -> D' and D' -> TD
   */
  template<
    size_t FM,
    size_t TM,
    size_t FD,
    size_t TD,
    size_t D
  >
  void
  intersect()
  {
    // std::cerr << "intersect: " << FD << " -> " << TD << std::endl;

    // The connectivity we will be populating
    connectivity_t & out_conn = get_connectivity_(FM, TM, FD, TD);
    if (!out_conn.empty()) {
      return;
    } // if

    // the number of each entity type
    auto num_from_ent = num_entities_(FD, FM);
    auto num_to_ent = num_entities_(TD, FM);

    // Temporary storage for connection id's
    connection_vector_t conns(num_from_ent);

    // Keep track of which to id's we have visited
    using visited_vec = std::vector<bool>;
    visited_vec visited(num_to_ent);

    size_t max_size = 1;

    // Read connectivities
    connectivity_t & c = get_connectivity_(FM, FD, D);
    assert(!c.empty());

    connectivity_t & c2 = get_connectivity_(TM, TD, D);
    assert(!c2.empty());
    
    // Iterate through entities in "from" topological dimension
    for(auto from_entity : entities<FD, FM>()){

      id_t from_id = from_entity->template global_id<FM>();
      id_vector_t & ents = conns[from_id.entity()];
      ents.reserve(max_size);

      size_t count;
      id_t * ep = c.get_entities(from_id.entity(), count);

      // Create a copy of to vertices so they can be sorted
      id_vector_t from_verts(ep, ep+count);
      // sort so we have a unique key for from vertices
      std::sort(from_verts.begin(), from_verts.end());

      // initially set all to id's to unvisited
      for (auto from_ent2 : entities<D, FM>(from_entity)) {
        for (id_t to_id : entity_ids<TD, TM>(from_ent2)) {
          visited[to_id.entity()] = false;
        }
      }

      // Loop through each from entity again
      for (auto from_ent2 : entities<D, FM>(from_entity)) {
        for (id_t to_id : entity_ids<TD, TM>(from_ent2)) {

          // If we have already visited, skip
          if (visited[to_id.entity()]) {
            continue;
          } // if

          visited[to_id.entity()] = true;

          // If the topological dimensions are the same, always add to id
          if (FD == TD) {
            if (from_id != to_id) {
              ents.push_back(to_id);
            } // if
          } else {
            size_t count;
            id_t * ep = c2.get_entities(to_id.entity(), count);

            // Create a copy of to vertices so they can be sorted
            id_vector_t to_verts(ep, ep + count);
            // Sort to verts so we can do an inclusion check
            std::sort(to_verts.begin(), to_verts.end());

            // If from vertices contains the to vertices add to id
            // to this connection set
            if (D < TD) {
              if (std::includes(from_verts.begin(), from_verts.end(),
                                  to_verts.begin(), to_verts.end()))
                ents.emplace_back(to_id);
            }
            // If we are going through a higher level, then set
            // intersection is sufficient. i.e. one set does not need to
            // be a subset of the other
            else {
              if (utils::intersects(from_verts.begin(), from_verts.end(),
                                      to_verts.begin(), to_verts.end()))
                ents.emplace_back(to_id);
            } // if

          } // if
        } // for
      } // for

      max_size = std::max(ents.size(), max_size);
    } // for

    // Finally create the connection from the temporary conns
    out_conn.init(conns);
  } // intersect

  /*!
     Used to compute connectivity information for topological dimension
       D1 -> D2
   */
  template<
    size_t M,
    size_t FD,
    size_t TD
  >
  void
  compute_connectivity()
  {
    // std::cerr << "compute: " << FD << " -> " << TD << std::endl;

    // Get the output connectivity
    connectivity_t & out_conn = get_connectivity_(M, FD, TD);

    // Check if we have already computed it
    if (!out_conn.empty()) {
      return;
    } // if

    // if we don't have cell -> vertex connectivities, then
    // try building cell -> vertex connectivity through the
    // faces (3d) or edges(2d)
    static_assert(MT::num_dimensions <= 3,
                   "this needs to be re-thought for higher dimensions");

    if (get_connectivity_(M, MT::num_dimensions, 0).empty()) {
      assert(!get_connectivity_(M, MT::num_dimensions-1, 0).empty() &&
              " need at least edges(2d)/faces(3) -> vertex connectivity");
      // assume we have cell -> faces, so invert it to get faces -> cells
      transpose<M, M, MT::num_dimensions-1, MT::num_dimensions>();
      // invert faces -> vertices to get vertices -> faces
      transpose<M, M, 0, MT::num_dimensions-1>();
      // build cells -> vertices via intersections with faces
      intersect<M, M, MT::num_dimensions, 0, MT::num_dimensions-1>();
    }

    // Check if we need to build entities, e.g: edges or faces
    if (num_entities_(FD, M) == 0) {
      if (get_connectivity_(M, FD+1, 0).empty())
        build_connectivity<M, FD, MT::num_dimensions>();
      else
        build_connectivity<M, FD, FD+1>();
    } // if

    if (num_entities_(TD, M) == 0) {
      if (get_connectivity_(M, TD+1, 0).empty())
        build_connectivity<M, TD, MT::num_dimensions>();
      else
        build_connectivity<M, TD, TD+1>();
    } // if

    if (num_entities_(FD, M) == 0 && num_entities_(TD, M) == 0) {
      return;
    } // if

    // Depending on the corresponding topological dimensions, call transpose
    // or intersect as need
     if (FD < TD) {
      compute_connectivity<M, TD, FD>();
      transpose<M, M, FD, TD>();
    } else {
       if (FD == 0 && TD == 0) {
         // compute vertex to vertex connectivities through shared cells.
         compute_connectivity<M, FD, MT::num_dimensions>();
         compute_connectivity<M, MT::num_dimensions, TD>();
         intersect<M, M, FD, TD, MT::num_dimensions>();
       } else {
         // computer connectivities through shared vertices.
         compute_connectivity<M, FD, 0>();
         compute_connectivity<M, 0, TD>();
         intersect<M, M, FD, TD, 0>();
       }
    } // if
  } // compute_connectivity

  /*!
    if the to-dimension is larger than the from-dimension, build the bindings
    using the create_bound_entities functionality
  */
  template<
    size_t FM,
    size_t TM,
    size_t FD,
    size_t TD,
    typename std::enable_if< (FM < TM) >::type* = nullptr
  >
  void
  _compute_bindings()
  {

    // if the connectivity for a transpose exists, do it
    if(!get_connectivity_(TM, FM, TD, FD).empty())
      transpose<FM, TM, FD, TD>();

    // otherwise try building the connectivity directly
    else if (num_entities_(TD, TM) == 0)
      build_bindings<FM, TM, TD>();

  } // compute_bindings

  /*!
    if the from-dimension is larger than the to-dimension, we want
    to transpose.  So make sure the opposite connectivity exists first
  */
  template<
    size_t FM,
    size_t TM,
    size_t FD,
    size_t TD,
    typename = typename std::enable_if< (FM > TM) >::type
  >
  void
  _compute_bindings()
  {

    // build the opposite connectivity first
    _compute_bindings< TM, FM, TD, FD>();

    // now apply a transpose to get the requested connectivity
    transpose<FM, TM, FD, TD>();

  } // compute_bindings

  /*!
    in the odd case the from-dimension matches the to-dimension, try and
    build the connectivity between the two
  */
  template<
    size_t FM,
    size_t TM,
    size_t FD,
    size_t TD
  >
  typename std::enable_if< (FM == TM) >::type
  _compute_bindings()
  {

    // compute connectivities through shared vertices at the at the lowest
    // dimension (doesn't matter which one really)
    _compute_bindings<0, TM, 0, FD>();
    _compute_bindings<0, TM, 0, TD>();

    // now try and transpose it
    auto & trans_conn = get_connectivity_(TM, FM, TD, FD);
    if(!trans_conn.empty())
      transpose<FM, TM, FD, TD>();

  } // compute_bindings


  /*!
    Main driver for computing bindings
  */
  template<
    size_t FM,
    size_t TM,
    size_t FD,
    size_t TD>
  void
  compute_bindings()
  {
    // std::cerr << "compute: , dom " << FM << " -> " << TM
    //           <<  ", dim " << FD << " -> " << TD << std::endl;

    // check if requested connectivity is already there, nothing to do
    connectivity_t & out_conn = get_connectivity_(FM, TM, FD, TD);

    if (!out_conn.empty()) return;

    _compute_bindings< FM, TM, FD, TD >();

  }

  /*!
    Build bindings associated with a from/to domain and topological dimension.
    compute_bindings will call this on each binding found in the tuple of
    bindings specified in the mesh type/traits mesh specialization.
   */
  template<
    size_t FM,
    size_t TM,
    size_t TD
  >
  void
  build_bindings()
  {

    // std::cerr << "build bindings: dom " << FM << " -> " << TM
    //           << " dim " << TD << std::endl;

    // Sanity check
    static_assert(TD <= MT::num_dimensions, "invalid dimension");

    // Helper variables
    size_t entity_id = 0;
    const size_t _num_cells = num_entities<MT::num_dimensions, FM>();

    // Storage for cell connectivity information
    connection_vector_t cell_conn(_num_cells);

    // Get cell definitions from domain 0
    auto & cells = base_t::ms_->index_spaces[FM][MT::num_dimensions];

    static constexpr size_t M0 = 0;

    for (size_t i = 0; i < MT::num_dimensions; ++i) {
      get_connectivity_<TM, FM, TD>(i).init();
    }
    for (size_t i = 0; i < TD; ++i) {
      get_connectivity_(TM, TM, TD, i).init();
    }

    // This buffer should be large enough to hold all entities
    // that potentially need to be created
    std::array<id_t, 4096> entity_ids;

    using to_entity_type = entity_type<TD, TM>;

    // Iterate over cells
    for (auto c : cells) {
      // Map used to ensure unique entity creation
      id_vector_map_t entity_ids_map;

      // Get a cell object.
      auto cell = static_cast<entity_type<MT::num_dimensions, M0> *>(c);
      id_t cell_id = cell->template global_id<FM>();

      domain_connectivity<MT::num_dimensions> & primal_conn =
        base_t::ms_->topology[FM][FM];
      domain_connectivity<MT::num_dimensions> & domain_conn =
        base_t::ms_->topology[FM][TM];

      // p.first:   The number of entities per cell.
      // p.second:  A std::vector of id_t containing the ids of the
      //            entities that define the bound entity.

      auto sv = cell->create_bound_entities(
        FM, TM, TD, cell_id, primal_conn, domain_conn, entity_ids.data());

      size_t n = sv.size();

      // Iterate over the newly-defined entities
      id_vector_t & conns = cell_conn[cell_id.entity()];

      conns.reserve(n);

      size_t pos = 0;

      auto& is = base_t::ms_->index_spaces[TM][TD].template cast<
        domain_entity<TM, to_entity_type>>();

      for (size_t i = 0; i < n; ++i) {
        size_t m = sv[i];

        id_t create_id = id_t::make<TD, TM>(entity_id, cell_id.partition());

        // Add this id to the cell entity connections
        conns.push_back(create_id);

        uint32_t dim_flags = 0;
        uint32_t dom_flags = 0;
        size_t num_vertices = 0;

        for (size_t k = 0; k < m; ++k) {
          id_t global_id = entity_ids[pos + k];
          auto dim = global_id.dimension();
          auto dom = global_id.domain();
          get_connectivity_(TM, dom, TD, dim).push(global_id);
          if (dom == FM) {
            dim_flags |= 1U << dim;
            num_vertices += dim == 0 ? 1 : 0;
          }
          else
            dom_flags |= 1U << dim;
        }

        for (size_t i = 0; i < MT::num_dimensions; ++i) {
          if (dim_flags & (1U << i)) {
            get_connectivity_<TM, FM, TD>(i).end_from();
          }
        }

        for (size_t i = 0; i < TD; ++i) {
          if (dom_flags & (1U << i)) {
            get_connectivity_(TM, TM, TD, i).end_from();
          }
        }

        auto ent = MT::template create_entity<TM, TD>(this, num_vertices);

        ++entity_id;

        pos += m;
      } // for
    } // for

    // Reference to storage from cells to the entity (to be created here).
    connectivity_t & cell_out =
      get_connectivity_(FM, TM, MT::num_dimensions, TD);
    cell_out.init(cell_conn);

  } // build_bindings

  /*!
   Implementation of get_connectivity for various get_connectivity convenience
   methods.
   */
  const connectivity_t &
  get_connectivity_(
    size_t from_domain,
    size_t to_domain,
    size_t from_dim,
    size_t to_dim) const
  {
    assert(from_domain < MT::num_domains && "invalid from domain");
    assert(to_domain < MT::num_domains && "invalid to domain");
    return base_t::ms_->topology[from_domain][to_domain].
      get(from_dim, to_dim);
  } // get_connectivity

  /*!
   Implementation of get_connectivity for various get_connectivity convenience
   methods.
   */
  connectivity_t &
  get_connectivity_(
    size_t from_domain,
    size_t to_domain,
    size_t from_dim,
    size_t to_dim)
  {
    assert(from_domain < MT::num_domains && "invalid from domain");
    assert(to_domain < MT::num_domains && "invalid to domain");
    return base_t::ms_->topology[from_domain][to_domain].
      get(from_dim, to_dim);
  } // get_connectivity

  /*!
   Implementation of get_connectivity for various get_connectivity convenience
   methods.
   */
  template<
    size_t FM,
    size_t TM,
    size_t FD
  >
  connectivity_t &
  get_connectivity_(
    size_t to_dim
  )
  {
    return base_t::ms_->topology[FM][TM].template get<FD>(to_dim);
  } // get_connectivity

  /*!
   Implementation of get_connectivity for various get_connectivity convenience
   methods.
   */
  template<
    size_t FM,
    size_t TM,
    size_t FD,
    size_t TD
  >
  connectivity_t &
  get_connectivity_()
  {
    return base_t::ms_->topology[FM][TM].template get<FD, TD>();
  } // get_connectivity

  /*!
   Implementation of get_connectivity for various get_connectivity convenience
   methods.
   */
  const connectivity_t &
  get_connectivity_(
    size_t domain,
    size_t from_dim,
    size_t to_dim) const
  {
    return get_connectivity_(domain, domain, from_dim, to_dim);
  } // get_connectivity

  connectivity_t &
  get_connectivity_(
    size_t domain,
    size_t from_dim,
    size_t to_dim)
  {
    return get_connectivity_(domain, domain, from_dim, to_dim);
  } // get_connectivity

}; // class mesh_topology_t

} // namespace topology
} // namespace flecsi

#endif // flecsi_topology_mesh_topology_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
