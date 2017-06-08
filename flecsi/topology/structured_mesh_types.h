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

#ifndef flecsi_structured_mesh_types_h
#define flecsi_structured_mesh_types_h

#include <array>
#include <unordered_map>
#include <cassert>
#include <iostream>
#include <vector>
#include <cmath>

#include "flecsi/data/data_client.h"
#include "flecsi/topology/mesh_utils.h"
#include "flecsi/topology/structured_index_space.h"


namespace flecsi {
namespace topology {  

template <typename T, T M>
struct typeify {
  static constexpr T value = M;
};

template <size_t M>
using domain_ = typeify<size_t, M>;


class structured_mesh_topology_base_t;

class structured_mesh_entity_base_{
public:
 using id_t = size_t;
};


template <size_t N>
class structured_mesh_entity_base_t : public structured_mesh_entity_base_ 
{
 public:
  virtual ~structured_mesh_entity_base_t() {}

  // get/set methods to get/set id_vector of the entity w.r.t specific domain
 /* template<size_t M> 
  id_vector_t get_id_vector() const
  {
    return idv_[M];
  }
  
  id_vector_t get_id_vector(size_t domain) const
  {
    return idv_[domain];
  }

  template<size_t M>
  void set_id_vector(const id_vector_t &id)
  {
    idv_[M] = id; 
  }
*/
  // get/set methods to get/set offset of the entity w.r.t specific domain
  template<size_t M>
  size_t offset() const
  {
    return offset_[M];
  }

  size_t offset(size_t domain) const
  {
    return offset_[domain];
  }

  template<size_t M>
  void set_offset(const size_t &id)
  {
    offset_[M] = id;
  }


  template <class MT>
  friend class structured_mesh_topology_t;

 private:
  std::array<id_t,N> offset_;
//  std::array<id_vector_t, N> idv_;

}; // class structured_mesh_entity_base_t


template <size_t D, size_t N>
class structured_mesh_entity_t : public structured_mesh_entity_base_t<N>
{
 public:
  static const size_t dimension = D;

  structured_mesh_entity_t() {}
  virtual ~structured_mesh_entity_t() {}
}; // class mesh_entity_t

/*
* D = num_dimensions, NM = num_domains
*/
template <size_t D, size_t NM>
struct structured_mesh_storage_t {

  using index_spaces_t = 
    std::array<structured_index_space<structured_mesh_entity_base_*>, size_t(pow(2,D))>;

  std::array<index_spaces_t, NM> index_spaces;

}; // struct mesh_storage_t


class structured_mesh_topology_base_t : public data::data_client_t
{
public:

  // Default constructor
  structured_mesh_topology_base_t() = default;

  // Don't allow the mesh to be copied or copy constructed
  structured_mesh_topology_base_t(const structured_mesh_topology_base_t &) = delete;
  structured_mesh_topology_base_t & operator=(const structured_mesh_topology_base_t &) = delete;

  /// Allow move operations
  structured_mesh_topology_base_t(structured_mesh_topology_base_t &&) = default;

  //! override default move assignement
  structured_mesh_topology_base_t & operator=(structured_mesh_topology_base_t && o)
  {
    // call base_t move operator
    data::data_client_t::operator=(std::move(o));
    // return a reference to the object
    return *this;
  };


  /*!
    Return the number of entities in for a specific domain and topology dim.
   */
  virtual size_t num_entities(size_t dim, size_t domain) const = 0;

  /*!
    This method should be called to construct an entity rather than
    calling the constructor directly. This way, the ability to have
    extra initialization behavior is reserved.
  */
  template <class T, class... S>
  T * make(S &&... args)
  {
    T * entity = new T(std::forward<S>(args)...);
    return entity;
  }; // make


}; // structured_mesh_topology_base_t


} // namespace topology
} // namespace flecsi

#endif // flecsi_structured_mesh_types_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
