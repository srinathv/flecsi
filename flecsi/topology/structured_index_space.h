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

#ifndef flecsi_structured_index_space_h
#define flecsi_structured_index_space_h

#include <cassert>
#include <algorithm>
#include <type_traits>
#include <vector>

namespace flecsi {
namespace topology {

using id_t = size_t;
using id_vector_t = std::vector<size_t>;  

//E is entity type
template<class E>
class structured_index_space{
public:

  // iterator 
  class iterator_t{
    public:
     iterator_t(id_t id):current{id}{};
    
     iterator_t& operator++(){
       ++current;
       return *this;
     };
  
    bool operator!=(const iterator_t& rhs){
        return (this->current != rhs.current);
    };
    
    id_t& operator*(){
        return current;
    };
   private:
     id_t current;
  };

  //constructors
// structured_index_space(const id_vector_t &lbnds,const id_vector_t &ubnds) 
  void init(const id_vector_t &lbnds, const id_vector_t &ubnds) 
   {
    assert(lbnds.size() == ubnds.size());

    size_t count = 1;
    for (size_t b = 0; b < lbnds.size(); ++b)
    {
      count *= ubnds[b]-lbnds[b]+1;
    }
   
    offset_ = 0;
    size_ = count;
    dim_ = lbnds.size();  
    lower_bnds_ = lbnds;
    upper_bnds_ = ubnds;
    
    std::cout<<std::endl;
    for (size_t b=0; b<lower_bnds_.size();++b){
      std::cout<<"lower_bnds["<<b<<"]="<<lower_bnds_[b]<<std::endl;
      std::cout<<"upper_bnds["<<b<<"]="<<upper_bnds_[b]<<std::endl;
    }
    
   }

  //structured_index_space(size_t offset,size_t size):offset_{offset}, size_{size}{};
  //structured_index_space(size_t size):offset_{0}, size_{size}{};

  structured_index_space(){};

  //provide slicing
  template <class S> 
  structured_index_space(const structured_index_space<S> &is)
  : offset_(is.offset_), size_(is.size_), dim_(is.dim_), lower_bnds_(is.lower_bnds_), upper_bnds_(is.upper_bnds_){};

  ~structured_index_space(){};

 
  auto begin()
  {
    return iterator_t(offset_);
  };

  auto end()
  {
  return iterator_t(offset_+size_);
  };
  
  size_t size() const 
  {
    return size_;
  }

  /*template<class S = E>
  auto slice()
  {
    return structured_index_space<S>(*this);
  }*/

  id_t get_offset_from_indices(id_vector_t &idv) 
  {
    //add range check for idv to make sure it lies within the bounds

    size_t value =0;
    size_t factor;

    for (int i = 0; i <dim_; ++i)
    {
      factor = 1;
      for (int j=i+1; j<dim_; ++j)
      {
        factor *= upper_bnds_[j]-lower_bnds_[j]+1;
      }

      value += idv[dim_-i-1]*factor;
    }

    return value;
  }


  auto get_indices_from_offset(id_t offset)
  {
    id_vector_t id(dim_);
    size_t factor;
    size_t rem = offset, value;


    std::cout<<"entity = "<<offset<<std::endl;

    for (int i=0; i<dim_; ++i)
    {
      factor = 1; 
      for (int j=i; j<dim_-1; ++j)
      {
       factor *= upper_bnds_[j]-lower_bnds_[j] + 1; 
      }
      value = rem/factor;
      id[dim_-i-1] = value;
      rem -= value*factor;
      
      std::cout<<"factor ="<<factor<<", value = "<<value<<", rem = "<<rem<<std::endl;
    }
 
   return id;
  }

 auto get_size_in_direction(size_t d){
  assert(d>=0 && d <dim_);
  return (upper_bnds_[d] - lower_bnds_[d]+1);
 }

template<size_t D>
bool check_index_limits(size_t index){
  return (index >= lower_bnds_[D] && index <= upper_bnds_[D]);
}
 
 private:
   id_t offset_;
   id_t size_;
   id_t dim_; 
   id_vector_t lower_bnds_;
   id_vector_t upper_bnds_;

};
} // namespace topology
} // namespace flecsi

#endif // flecsi_structured_index_space_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
