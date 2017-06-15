/*~--------------------------------------------------------------------------~*
 *  @@@@@@@@  @@           @@@@@@   @@@@@@@@ @@
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

#ifndef flecsi_structured_mesh_topology_h
#define flecsi_structured_mesh_topology_h


#include <cassert>
#include <iostream>
#include <cstring>
#include <array>
#include <vector>
#include <map>
#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <functional>
#include <type_traits>

#include "flecsi/utils/common.h"
#include "flecsi/utils/set_intersection.h"
#include "flecsi/utils/static_verify.h"
#include "flecsi/topology/structured_mesh_types.h"

namespace flecsi {
namespace topology {
namespace verify_mesh {

  FLECSI_MEMBER_CHECKER(num_dimensions);
  FLECSI_MEMBER_CHECKER(num_domains);
  FLECSI_MEMBER_CHECKER(lower_bounds);
  FLECSI_MEMBER_CHECKER(upper_bounds);
//  FLECSI_MEMBER_CHECKER(entity_types);

} // namespace verify_mesh

template<
  class MT
>
class structured_mesh_topology_t : public structured_mesh_topology_base_t
{
 
  /*
  * Verify the existence of following fields in the mesh policy MT
  * num_dimensions
  * num_domains
  * entity_types
  * mesh_bounds (id-based)
  * mesh_resolution
  * domain specific ordering policy
  */
  // static verification of mesh policy

  static_assert(verify_mesh::has_member_num_dimensions<MT>::value,
                "mesh policy missing num_dimensions size_t");
  
  static_assert(std::is_convertible<decltype(MT::num_dimensions),
    size_t>::value, "mesh policy num_dimensions must be size_t");

  static_assert(verify_mesh::has_member_num_domains<MT>::value,
                "mesh policy missing num_domains size_t");
  
  static_assert(std::is_convertible<decltype(MT::num_domains),
    size_t>::value, "mesh policy num_domains must be size_t");

  static_assert(verify_mesh::has_member_lower_bounds<MT>::value,
                "mesh policy missing lower_bounds array");
  
  static_assert(std::is_convertible<decltype(MT::lower_bounds),
    std::array<size_t,MT::num_dimensions>>::value,"mesh policy lower_bounds is not an array");

  static_assert(verify_mesh::has_member_upper_bounds<MT>::value,
                "mesh policy missing upper_bounds array");
  
  static_assert(std::is_convertible<decltype(MT::upper_bounds),
    std::array<size_t,MT::num_dimensions>>::value,"mesh policy upper_bounds is not an array");


  static_assert(MT::lower_bounds.size() == MT::upper_bounds.size(),
     "mesh bounds have inconsistent sizes");

//  static_assert(verify_mesh::has_member_entity_types<MT>::value,
//                "mesh policy missing entity_types tuple");
  
//  static_assert(utils::is_tuple<typename MT::entity_types>::value,
//                "mesh policy entity_types is not a tuple");


public:
  // used to find the entity type of topological dimension D and domain M
  //template<size_t D, size_t M = 0>
  //using entity_type = typename find_entity_<MT, D, M>::type;

  using id_vector_t = std::vector<size_t>;  

  // Don't allow the mesh to be copied or copy constructed
  structured_mesh_topology_t(const structured_mesh_topology_t &) = delete;
  structured_mesh_topology_t & operator=(const structured_mesh_topology_t &) = delete;

  // Allow move operations
  structured_mesh_topology_t(structured_mesh_topology_t && o) = default;
  structured_mesh_topology_t & operator=(structured_mesh_topology_t && o) = default;

  //! Constructor
  structured_mesh_topology_t()
  {
      meshdim_ = MT::num_dimensions;  
      meshbnds_low_.insert(meshbnds_low_.begin(),  MT::lower_bounds.begin(), MT::lower_bounds.end());
      meshbnds_up_.insert(meshbnds_up_.begin(),  MT::upper_bounds.begin(), MT::upper_bounds.end());

      id_vector_t vec, ubnds;      
      for (size_t i = 0; i < pow(2,meshdim_); ++i){
        ubnds.clear();
        vec = get_bnds(meshdim_, i);

        for (size_t j = 0; j < meshdim_; ++j)
            ubnds.push_back(MT::upper_bounds[j] - vec[j]);
            
       // for (int k=0; k<ubnds.size(); k++)
       //   std::cout<<"ubnds["<<k<<"] = "<<ubnds[k]<<std::endl;    
      
        ms_.index_spaces[0][i].init({MT::lower_bounds.begin(), MT::lower_bounds.end()}, ubnds);
        
        std::cout<<"IS-size = "<<ms_.index_spaces[0][i].template size()<<std::endl;
      }
}

  // mesh destructor
  virtual ~structured_mesh_topology_t(){}
  
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
    return ms_.index_spaces[M][D].size();
  } // num_entities

  template<
    size_t D,
    size_t M = 0
  >
  auto get_indices(id_t entity) 
  {
    return ms_.index_spaces[M][D].template get_indices_from_offset(entity);
  }
  
 template<
    size_t D,
    size_t M = 0
  >
  auto get_offset(id_vector_t &idv) 
  {
    return ms_.index_spaces[M][D].template get_offset_from_indices(idv);
  }

  auto lower_bounds()
  {
    return meshbnds_low_;
  }

  auto upper_bounds()
  {
    return meshbnds_up_;
  }

  //Query type 1: Provides traversal over whole index space
  /*!
    Get the top-level entities of topological dimension D of the specified
    domain M. e.g: cells of the mesh.
  */
  template<
    size_t D,
    size_t M = 0
  >
  auto
  entities()
  {
    return ms_.index_spaces[M][D]; 
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
   return ms_.index_spaces[M][D];
  } // entities


  //Query type 2: FEM-type adjacency queries 
  /*!
    Given an entity of a specific dimension FD, find all entities incident 
    on it from dimension TD.
    Supports queries between non-equal dimensions. The queries for entities from the 
    same dimension can be obtained through the more finer level queries for intra-index
    space. 
  */

  template<size_t FD, size_t TD, size_t N>
  auto get_entities(size_t ent)
  {
    assert(FD != TD);
    id_vector_t adj;
    
    if (FD < TD) 
      adj = entities_up_<FD,TD,N>(ent);
    else if (FD > TD)
      adj = entities_down_<FD,TD,N>(ent);
    else
      std::cerr<<("Requesting adjacencies between non-valid dimensions \n");       

    return adj;

  } //entities


/* Query type 3: FD-type stencil queries between entities
 *               of the same dimension. 
 * D : dimension, N: domain
 * xoff, yoff, zoff: offsets w.r.t current entity
*/
template<size_t D, size_t N, std::intmax_t xoff>
auto entities(size_t ent)
{
    assert(xoff != 0); 
    size_t value = ent;
    auto indices = ms_.index_spaces[N][D].template get_indices_from_offset(ent);
    if (ms_.index_spaces[N][D].template check_index_limits<0>(xoff+indices[0]))
      value += xoff;
   // else 
     // value = -1;
    return value;
} //entities

template<size_t D, size_t N, std::intmax_t xoff, std::intmax_t yoff>
auto entities(size_t ent)
{
    assert(!((xoff == 0) && (yoff == 0))); 
    size_t ind = get_index_in_storage(ent,D,N);
    size_t value = ent;
    size_t nx;

    if (ind == 2)
      ent = ent - ms_.index_spaces[N][1].template size();

    auto indices = ms_.index_spaces[N][ind].template get_indices_from_offset(ent);
    
    if((ms_.index_spaces[N][ind].template check_index_limits<0>(xoff+indices[0]))
       && (ms_.index_spaces[N][ind].template check_index_limits<1>(yoff+indices[1])))
    {
      nx = ms_.index_spaces[N][ind].template get_size_in_direction(0);
      value += xoff + nx*yoff;
    }
  return value; 
} //entities

template<size_t D, size_t N, std::intmax_t xoff, std::intmax_t yoff, std::intmax_t zoff>
auto entities(size_t ent)
{
    assert(!((xoff == 0) && (yoff == 0) && (zoff == 0))); 
    size_t ind = get_index_in_storage(ent,D,N);
    size_t value = ent;
    size_t nx, ny;
    
    if (ind == 2) //edge-y
    {
      ent = ent - ms_.index_spaces[N][1].template size();
    }
    else if (ind == 3)//edge-z
    {
      ent = ent - (ms_.index_spaces[N][1].template size()) - (ms_.index_spaces[N][2].template size());
    }
    else if (ind == 5)//face-y
    {
      ent = ent - ms_.index_spaces[N][4].template size();
    }
    else if (ind == 6)//face-z
    {
      ent = ent - (ms_.index_spaces[N][4].template size()) - (ms_.index_spaces[N][5].template size());
    }
    else
      ent = ent;

    auto indices = ms_.index_spaces[N][ind].template get_indices_from_offset(ent);
    if((ms_.index_spaces[N][ind].template check_index_limits<0>(xoff+indices[0]))
       && (ms_.index_spaces[N][ind].template check_index_limits<1>(yoff+indices[1]))
       && (ms_.index_spaces[N][ind].template check_index_limits<2>(zoff+indices[2])))
     {
         nx = ms_.index_spaces[N][ind].template get_size_in_direction(0);
         ny = ms_.index_spaces[N][ind].template get_size_in_direction(1);
         value += xoff + nx*yoff+nx*ny*zoff;
     }
     //else
     //   value = -1;
    return value;
} //entities


private:

  size_t meshdim_; 
  id_vector_t meshbnds_low_;
  id_vector_t meshbnds_up_;
  structured_mesh_storage_t<MT::num_dimensions, MT::num_domains> ms_;


  // Get the number of entities in a given domain and topological dimension
  size_t
  num_entities_(
    size_t dim,
    size_t domain=0
  ) const
  {
    return ms_.index_spaces[domain][dim].size();
  } // num_entities_


  // Up adjacencies
  template<size_t FD, size_t TD, size_t N>
  auto entities_up_(size_t ent)
  {
    id_vector_t adj;
    
    if (meshdim_ == 1)
    {
      adj = entities_up_1D_<FD,TD,N>(ent);
    }
    else if (meshdim_ == 2)
    {
     adj = entities_up_2D_<FD,TD,N>(ent);
    }
    else if (meshdim_ == 3)
    {
      adj = entities_up_3D_<FD,TD,N>(ent);
    }
    else
     std::cerr<<"requesting adjacencies for entities with dimension greater than 3\n";

   return adj;    
  } //entities_up_

  // Down adjacencies
  template<size_t FD, size_t TD, size_t N>
  auto entities_down_(size_t ent)
  {
    id_vector_t adj;
    
    if (meshdim_ == 1)
    {
      adj = entities_down_1D_<FD,TD,N>(ent);
    }
    else if (meshdim_ == 2)
    {
     adj = entities_down_2D_<FD,TD,N>(ent);
    }
    else if (meshdim_ == 3)
    {
      adj = entities_down_3D_<FD,TD,N>(ent);
    }
    else
     std::cerr<<"requesting adjacencies for entities with dimension greater than 3\n";

   return adj;    
  } //entities_down_

   /****************************************************
   *        Lowest level 1D up/down-adjacencies        *
   *****************************************************/
  template<size_t FD, size_t TD, size_t N>
  auto entities_up_1D_(size_t ent)
  {
    id_vector_t adj, indices;
    size_t index = get_index_in_storage(ent, FD, N);
    indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent);
    
    assert(FD == 0);
    auto offset = ms_.index_spaces[N][TD].template get_offset_from_indices(indices);

    adj.push_back(offset);
    adj.push_back(entities<TD,N,-1>(offset));
    return adj;
  } //entities_up_1D_


  template<size_t FD, size_t TD, size_t N>
  auto entities_down_1D_(size_t ent)
  {
    id_vector_t adj, indices;
    size_t index = get_index_in_storage(ent, FD, N);
    indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent);
   
    auto offset = ms_.index_spaces[N][TD].template get_offset_from_indices(indices);

    adj.push_back(offset);
    adj.push_back(entities<TD,N,1>(offset));
    return adj;
  }//entities_down_1D_


  /****************************************************
   *       Lowest level 2D up/down-adjacencies        *
   ****************************************************/
  template<size_t FD, size_t TD, size_t N>
  auto entities_up_2D_(size_t ent)
  {
    id_vector_t adj, indices, ngb_indices;
    size_t index = get_index_in_storage(ent, FD, N);
    //indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent);

    size_t offset, xoffset, yoffset, nx;  

    if ((FD == 0) && (TD == 1)) //V-->E
    {
      indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent);
      nx = ms_.index_spaces[N][TD].template size();
 
      if (indices[0] <= (ms_.index_spaces[N][TD+1].template max<0>()))
      { 
        yoffset = ms_.index_spaces[N][TD+1].template get_offset_from_indices(indices);
        adj.push_back(yoffset+nx);
      }        

      if (indices[1] <= (ms_.index_spaces[N][TD].template max<1>()))
      {
        xoffset = ms_.index_spaces[N][TD].template get_offset_from_indices(indices);
        adj.push_back(xoffset);
      }

      if (indices[0] >= (ms_.index_spaces[N][TD+1].template min<0>()+1))
      { 
        ngb_indices = indices;
        ngb_indices[0] = ngb_indices[0] - 1;
        yoffset = ms_.index_spaces[N][TD+1].template get_offset_from_indices(ngb_indices);
        adj.push_back(yoffset+nx);
      }        

      if (indices[1] >= (ms_.index_spaces[N][TD].template min<1>()+1))
      {
        ngb_indices = indices;
        ngb_indices[1] = ngb_indices[1] - 1;
        xoffset = ms_.index_spaces[N][TD].template get_offset_from_indices(ngb_indices);
        adj.push_back(xoffset); 
      }
    }
    else if ((FD == 0) && (TD == 2)) //V-->F
    {
      indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent);
      if ( (indices[0] <= (ms_.index_spaces[N][TD+1].template max<0>())) && 
           (indices[1] <= (ms_.index_spaces[N][TD+1].template max<1>())) )
      {
        offset = ms_.index_spaces[N][TD+1].template get_offset_from_indices(indices);
        adj.push_back(offset);
      }

     if ( (indices[0] >= (ms_.index_spaces[N][TD+1].template min<0>()+1)) && 
           (indices[1] <= (ms_.index_spaces[N][TD+1].template max<1>())) )
      {
        ngb_indices = indices;
        ngb_indices[0] = ngb_indices[0] - 1;
        offset = ms_.index_spaces[N][TD+1].template get_offset_from_indices(ngb_indices);
        adj.push_back(offset);
      }

      if ( (indices[0] >= (ms_.index_spaces[N][TD+1].template min<0>()+1)) && 
           (indices[1] >= (ms_.index_spaces[N][TD+1].template min<1>()+1)) )
      {
        ngb_indices = indices;
        ngb_indices[0] = ngb_indices[0] - 1;
        ngb_indices[1] = ngb_indices[1] - 1;
        offset = ms_.index_spaces[N][TD+1].template get_offset_from_indices(ngb_indices);
        adj.push_back(offset);
      }

      if ( (indices[0] <= (ms_.index_spaces[N][TD+1].template max<0>())) && 
           (indices[1] >= (ms_.index_spaces[N][TD+1].template min<1>()+1)) )
      {
        ngb_indices = indices;
        ngb_indices[1] = ngb_indices[1] - 1;
        offset = ms_.index_spaces[N][TD+1].template get_offset_from_indices(ngb_indices);
        adj.push_back(offset);
      }
    }
    else if ((FD == 1) && (TD == 2)) //E-->F
    {
      if (ent < ms_.index_spaces[N][FD].template size()) //Ex -->F
      {
        indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent);
        if (indices[0] <= (ms_.index_spaces[N][TD+1].template max<0>()))
        {
          offset = ms_.index_spaces[N][TD+1].template get_offset_from_indices(indices);
          adj.push_back(offset);
        }

        if (indices[0] >= (ms_.index_spaces[N][TD+1].template min<0>()+1))
        {
          indices[0] = indices[0] - 1;
          offset = ms_.index_spaces[N][TD+1].template get_offset_from_indices(indices);
          adj.push_back(offset);
        }
      }
      else //Ey -->F
      {
        nx = ms_.index_spaces[N][FD].template size();
        indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent-nx);
        if (indices[1] <= (ms_.index_spaces[N][TD+1].template max<1>()))
        {
          offset = ms_.index_spaces[N][TD+1].template get_offset_from_indices(indices);
          adj.push_back(offset);
        }

        if (indices[1] >= (ms_.index_spaces[N][TD+1].template min<1>()+1))
        {
          indices[1] = indices[1] - 1;
          offset = ms_.index_spaces[N][TD+1].template get_offset_from_indices(indices);
          adj.push_back(offset);
        }
      } 
    }
    else
      std::cerr<<"errs msg"<<std::endl;

    return adj;
  }//entities_up_2D_

  template<size_t FD, size_t TD, size_t N>
  auto entities_down_2D_(size_t ent)
  {
    id_vector_t adj, indices;
    size_t index = get_index_in_storage(ent, FD, N);  
    size_t offset, xoffset, yoffset, nx;  

    if ((FD == 1) && (TD == 0)) //E-->V
    {
      if (ent < ms_.index_spaces[N][FD].template size()) //Ex-->V
       {
          indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent);
          offset = ms_.index_spaces[N][TD].template get_offset_from_indices(indices);
          adj.push_back(offset);
          adj.push_back(entities<TD,N,0,1>(offset));
       }
       else //Ey-->V
       {
          size_t nx = ms_.index_spaces[N][FD].template size();
          indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent-nx);
          offset = ms_.index_spaces[N][TD].template get_offset_from_indices(indices);   
          adj.push_back(offset);
          adj.push_back(entities<TD,N,1,0>(offset));
       }
    }
    else if ((FD == 2) && (TD == 0)) //F-->V
    {
      indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent);
      offset = ms_.index_spaces[N][TD].template get_offset_from_indices(indices);
      adj.push_back(offset);
      adj.push_back(entities<TD,N,1,0>(offset));
      adj.push_back(entities<TD,N,1,1>(offset));
      adj.push_back(entities<TD,N,0,1>(offset));
    }
    else if ((FD == 2) && (TD == 1)) //F-->E
    {
      indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent);
      nx = ms_.index_spaces[N][TD].template size();
      yoffset = ms_.index_spaces[N][TD+1].template get_offset_from_indices(indices);
      xoffset = ms_.index_spaces[N][TD].template get_offset_from_indices(indices);  
      
      adj.push_back(yoffset+nx);
      adj.push_back(entities<TD,N,1,0>(xoffset)); 
      adj.push_back(entities<TD,N,0,1>(yoffset+nx));     
      adj.push_back(xoffset);
    }

    return adj;
  }//entities_down_2D_



  /****************************************************
   *       Lowest level 3D up/down-adjacencies        *
   ****************************************************/
  template<size_t FD, size_t TD, size_t N>
  auto entities_up_3D_(size_t ent)
  {
    id_vector_t adj, indices, ngb_indices;
    size_t index = get_index_in_storage(ent, FD, N);
   
    size_t offset, xoffset, yoffset, zoffset, nx, ny;  

    if ((FD == 0) && (TD == 1)) //V-->E
    {
      indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent);
      nx = ms_.index_spaces[N][TD].template size();
      ny = ms_.index_spaces[N][TD+1].template size();
 
      if (indices[1] <= (ms_.index_spaces[N][TD].template max<1>()))
      {
        xoffset = ms_.index_spaces[N][TD].template get_offset_from_indices(indices);
        adj.push_back(xoffset);//Ex(i,j,k)
      }
      
     if (indices[1] >= (ms_.index_spaces[N][TD].template min<1>()+1))
      {
        ngb_indices = indices;
        ngb_indices[1] = ngb_indices[1] - 1;
        xoffset = ms_.index_spaces[N][TD].template get_offset_from_indices(ngb_indices);
        adj.push_back(xoffset); //Ex(i,j-1,k)
      }
     
      if (indices[0] <= (ms_.index_spaces[N][TD+1].template max<0>()))
      { 
        yoffset = ms_.index_spaces[N][TD+1].template get_offset_from_indices(indices);
        adj.push_back(yoffset+nx);//Ey(i,j,k)
      }        


      if (indices[0] >= (ms_.index_spaces[N][TD+1].template min<0>()+1))
      { 
        ngb_indices = indices;
        ngb_indices[0] = ngb_indices[0] - 1;
        yoffset = ms_.index_spaces[N][TD+1].template get_offset_from_indices(ngb_indices);
        adj.push_back(yoffset+nx);//Ey(i-1,j,k)
      }        

      
      if (indices[2] <= (ms_.index_spaces[N][TD+2].template max<2>()))
      { 
        zoffset = ms_.index_spaces[N][TD+2].template get_offset_from_indices(indices);
        adj.push_back(zoffset+nx+ny);//Ez(i,j,k)
      }        

       if (indices[2] >= (ms_.index_spaces[N][TD+2].template min<2>()+1))
      { 
        ngb_indices = indices;
        ngb_indices[2] = ngb_indices[2] - 1;
        zoffset = ms_.index_spaces[N][TD+2].template get_offset_from_indices(ngb_indices);
        adj.push_back(zoffset+nx+ny);//Ez(i,j,k-1)
      }  
    }
    else if ((FD == 0) && (TD == 2)) //V-->F
    {
      indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent);
      nx = ms_.index_spaces[N][TD+2].template size();
      ny = ms_.index_spaces[N][TD+3].template size();
      
      //V-->Fx
      if ( (indices[1] <= (ms_.index_spaces[N][TD+2].template max<1>())) && 
           (indices[2] <= (ms_.index_spaces[N][TD+2].template max<2>())) )
      {
        offset = ms_.index_spaces[N][TD+2].template get_offset_from_indices(indices);
        adj.push_back(offset);//Fx(i,j,k)
      }

     if ( (indices[1] >= (ms_.index_spaces[N][TD+2].template min<1>()+1)) && 
           (indices[2] <= (ms_.index_spaces[N][TD+2].template max<2>())) )
      {
        ngb_indices = indices;
        ngb_indices[1] = ngb_indices[1] - 1;
        offset = ms_.index_spaces[N][TD+2].template get_offset_from_indices(ngb_indices);
        adj.push_back(offset);//Fx(i,j-1,k)
      }

      if ( (indices[1] >= (ms_.index_spaces[N][TD+2].template min<1>()+1)) && 
           (indices[2] >= (ms_.index_spaces[N][TD+2].template min<2>()+1)) )
      {
        ngb_indices = indices;
        ngb_indices[1] = ngb_indices[1] - 1;
        ngb_indices[2] = ngb_indices[2] - 1;
        offset = ms_.index_spaces[N][TD+2].template get_offset_from_indices(ngb_indices);
        adj.push_back(offset);//Fx(i,j-1,k-1)
      }

      if ( (indices[1] <= (ms_.index_spaces[N][TD+2].template max<1>())) && 
           (indices[2] >= (ms_.index_spaces[N][TD+2].template min<2>()+1)) )
      {
        ngb_indices = indices;
        ngb_indices[2] = ngb_indices[2] - 1;
        offset = ms_.index_spaces[N][TD+2].template get_offset_from_indices(ngb_indices);
        adj.push_back(offset);//Fx(i,j,k-1)
      }
       
      //V-->Fy
      if ( (indices[0] <= (ms_.index_spaces[N][TD+3].template max<0>())) && 
           (indices[2] <= (ms_.index_spaces[N][TD+3].template max<2>())) )
      {
        offset = ms_.index_spaces[N][TD+3].template get_offset_from_indices(indices);
        adj.push_back(offset+nx);//Fy(i,j,k)
      }

     if ( (indices[0] >= (ms_.index_spaces[N][TD+3].template min<0>()+1)) && 
           (indices[2] <= (ms_.index_spaces[N][TD+3].template max<2>())) )
      {
        ngb_indices = indices;
        ngb_indices[0] = ngb_indices[0] - 1;
        offset = ms_.index_spaces[N][TD+3].template get_offset_from_indices(ngb_indices);
        adj.push_back(offset+nx);//Fy(i-1,j,k)
      }

      if ( (indices[0] >= (ms_.index_spaces[N][TD+3].template min<0>()+1)) && 
           (indices[2] >= (ms_.index_spaces[N][TD+3].template min<2>()+1)) )
      {
        ngb_indices = indices;
        ngb_indices[0] = ngb_indices[0] - 1;
        ngb_indices[2] = ngb_indices[2] - 1;
        offset = ms_.index_spaces[N][TD+3].template get_offset_from_indices(ngb_indices);
        adj.push_back(offset+nx);//Fy(i-1,j,k-1)
      }

      if ( (indices[0] <= (ms_.index_spaces[N][TD+3].template max<0>())) && 
           (indices[2] >= (ms_.index_spaces[N][TD+3].template min<2>()+1)) )
      {
        ngb_indices = indices;
        ngb_indices[2] = ngb_indices[2] - 1;
        offset = ms_.index_spaces[N][TD+3].template get_offset_from_indices(ngb_indices);
        adj.push_back(offset+nx);//Fy(i,j,k-1)
      }
      
      //V-->Fz
      if ( (indices[0] <= (ms_.index_spaces[N][TD+4].template max<0>())) && 
           (indices[1] <= (ms_.index_spaces[N][TD+4].template max<1>())) )
      {
        offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(indices);
        adj.push_back(offset+nx+ny);//Fz(i,j,k)
      }

     if ( (indices[0] >= (ms_.index_spaces[N][TD+4].template min<0>()+1)) && 
           (indices[1] <= (ms_.index_spaces[N][TD+4].template max<1>())) )
      {
        ngb_indices = indices;
        ngb_indices[0] = ngb_indices[0] - 1;
        offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(ngb_indices);
        adj.push_back(offset+nx+ny);//Fz(i-1,j,k)
      }

      if ( (indices[0] >= (ms_.index_spaces[N][TD+4].template min<0>()+1)) && 
           (indices[1] >= (ms_.index_spaces[N][TD+4].template min<1>()+1)) )
      {
        ngb_indices = indices;
        ngb_indices[0] = ngb_indices[0] - 1;
        ngb_indices[1] = ngb_indices[1] - 1;
        offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(ngb_indices);
        adj.push_back(offset+nx+ny);//Fz(i-1,j-1,k)
      }

      if ( (indices[0] <= (ms_.index_spaces[N][TD+4].template max<0>())) && 
           (indices[1] >= (ms_.index_spaces[N][TD+4].template min<1>()+1)) )
      {
        ngb_indices = indices;
        ngb_indices[1] = ngb_indices[1] - 1;
        offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(ngb_indices);
        adj.push_back(offset+nx+ny);//Fz(i,j-1,k)
      }  
    }
    else if ((FD == 0) && (TD == 3)) //V-->C
    {
      indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent);
      
      if ( (indices[0] <= (ms_.index_spaces[N][TD+4].template max<0>())) && 
           (indices[1] <= (ms_.index_spaces[N][TD+4].template max<1>())) &&
           (indices[2] <= (ms_.index_spaces[N][TD+4].template max<2>()))) 
      {
        offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(indices);
        adj.push_back(offset);//C(i,j,k)
      }

     if ( (indices[0] >= (ms_.index_spaces[N][TD+4].template min<0>()+1)) && 
           (indices[1] <= (ms_.index_spaces[N][TD+4].template max<1>()))  &&
           (indices[2] <= (ms_.index_spaces[N][TD+4].template max<2>())))
      {
        ngb_indices = indices;
        ngb_indices[0] = ngb_indices[0] - 1;
        offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(ngb_indices);
        adj.push_back(offset);//C(i-1,j,k)
      }

      if ( (indices[0] >= (ms_.index_spaces[N][TD+4].template min<0>()+1)) && 
           (indices[1] >= (ms_.index_spaces[N][TD+4].template min<1>()+1))  &&
           (indices[2] <= (ms_.index_spaces[N][TD+4].template max<2>())))
      {
        ngb_indices = indices;
        ngb_indices[0] = ngb_indices[0] - 1;
        ngb_indices[1] = ngb_indices[1] - 1;
        offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(ngb_indices);
        adj.push_back(offset);//C(i-1,j-1,k)
      }

      if ( (indices[0] <= (ms_.index_spaces[N][TD+4].template max<0>())) && 
           (indices[1] >= (ms_.index_spaces[N][TD+4].template min<1>()+1)) &&
           (indices[2] <= (ms_.index_spaces[N][TD+4].template max<2>())))
      {
        ngb_indices = indices;
        ngb_indices[1] = ngb_indices[1] - 1;
        offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(ngb_indices);
        adj.push_back(offset);//C(i,j-1,k)
      }
      
      if ( (indices[0] <= (ms_.index_spaces[N][TD+4].template max<0>())) && 
           (indices[1] <= (ms_.index_spaces[N][TD+4].template max<1>())) &&
           (indices[2] >= (ms_.index_spaces[N][TD+4].template max<2>()+1))) 
      {
        ngb_indices = indices;
        ngb_indices[2] = ngb_indices[2] - 1;
        offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(ngb_indices);
        adj.push_back(offset);//C(i,j,k-1)
      }

     if ( (indices[0] >= (ms_.index_spaces[N][TD+4].template min<0>()+1)) && 
           (indices[1] <= (ms_.index_spaces[N][TD+4].template max<1>()))  &&
           (indices[2] >= (ms_.index_spaces[N][TD+4].template max<2>()+1)))
      {
        ngb_indices = indices;
        ngb_indices[0] = ngb_indices[0] - 1;
        ngb_indices[2] = ngb_indices[2] - 1;
        offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(ngb_indices);
        adj.push_back(offset);//C(i-1,j,k-1)
      }

      if ( (indices[0] >= (ms_.index_spaces[N][TD+4].template min<0>()+1)) && 
           (indices[1] >= (ms_.index_spaces[N][TD+4].template min<1>()+1))  &&
           (indices[2] >= (ms_.index_spaces[N][TD+4].template max<2>()+1)))
      {
        ngb_indices = indices;
        ngb_indices[0] = ngb_indices[0] - 1;
        ngb_indices[1] = ngb_indices[1] - 1;
        ngb_indices[2] = ngb_indices[2] - 1;
        offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(ngb_indices);
        adj.push_back(offset);//C(i-1,j-1,k-1)
      }

      if ( (indices[0] <= (ms_.index_spaces[N][TD+4].template max<0>())) && 
           (indices[1] >= (ms_.index_spaces[N][TD+4].template min<1>()+1)) &&
           (indices[2] >= (ms_.index_spaces[N][TD+4].template max<2>()+1)))
      {
        ngb_indices = indices;
        ngb_indices[1] = ngb_indices[1] - 1;
        ngb_indices[2] = ngb_indices[2] - 1;
        offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(ngb_indices);
        adj.push_back(offset);//C(i,j-1,k-1)
      }
    }
    else if ((FD == 1) && (TD == 2)) //E-->F
    {
      if (ent < ms_.index_spaces[N][FD].template size()) //Ex -->F
      {
        indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent);
        nx = ms_.index_spaces[N][TD+2].template size();
        ny = ms_.index_spaces[N][TD+3].template size();
        
        if (indices[0] <= (ms_.index_spaces[N][TD+4].template max<0>()) &&
            indices[2] <= (ms_.index_spaces[N][TD+4].template max<2>()) )
        {
          offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(indices);
          adj.push_back(offset+nx+ny);
        } //Fz(i,j,k)
        
        if (indices[0] <= (ms_.index_spaces[N][TD+2].template max<0>()) &&
            indices[2] <= (ms_.index_spaces[N][TD+2].template max<2>()) )
        {
          offset = ms_.index_spaces[N][TD+2].template get_offset_from_indices(indices);
          adj.push_back(offset);
        }//Fx(i,j,k)
        
         if (indices[0] >= (ms_.index_spaces[N][TD+4].template min<0>()+1) &&
            indices[2] <= (ms_.index_spaces[N][TD+4].template max<2>()) )
        {
          ngb_indices = indices;
          ngb_indices[0] = ngb_indices[0] - 1;
          offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(indices);
          adj.push_back(offset+nx+ny);
        } //Fz(i-1,j,k)
        
        if (indices[0] <= (ms_.index_spaces[N][TD+2].template max<0>()) &&
            indices[2] >= (ms_.index_spaces[N][TD+2].template min<2>()+1) )
        {
          ngb_indices = indices;
          ngb_indices[2] = ngb_indices[2] - 1;
          offset = ms_.index_spaces[N][TD+2].template get_offset_from_indices(indices);
          adj.push_back(offset);
        }//Fx(i,j,k-1)
      }
      else if ((ent >= ms_.index_spaces[N][FD].template size()) && 
       (ent <( ms_.index_spaces[N][FD+1].template size())+
       ( ms_.index_spaces[N][FD].template size())) ) //Ey -->F
      {
        nx = ms_.index_spaces[N][FD].template size();
        indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent-nx);
        nx = ms_.index_spaces[N][TD+2].template size();
        ny = ms_.index_spaces[N][TD+3].template size();
        
        if (indices[1] <= (ms_.index_spaces[N][TD+4].template max<1>()) &&
            indices[2] <= (ms_.index_spaces[N][TD+4].template max<2>()) )
        {
          offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(indices);
          adj.push_back(offset+nx+ny);
        } //Fz(i,j,k)
        
        if (indices[1] <= (ms_.index_spaces[N][TD+3].template max<1>()) &&
            indices[2] <= (ms_.index_spaces[N][TD+3].template max<2>()) )
        {
          offset = ms_.index_spaces[N][TD+3].template get_offset_from_indices(indices);
          adj.push_back(offset+nx);
        }//Fy(i,j,k)
        
         if (indices[1] >= (ms_.index_spaces[N][TD+4].template min<1>()+1) &&
            indices[2] <= (ms_.index_spaces[N][TD+4].template max<2>()) )
        {
          ngb_indices = indices;
          ngb_indices[1] = ngb_indices[1] - 1;
          offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(indices);
          adj.push_back(offset+nx+ny);
        } //Fz(i,j-1,k)
        
        if (indices[1] <= (ms_.index_spaces[N][TD+3].template max<1>()) &&
            indices[2] >= (ms_.index_spaces[N][TD+3].template min<2>()+1) )
        {
          ngb_indices = indices;
          ngb_indices[2] = ngb_indices[2] - 1;
          offset = ms_.index_spaces[N][TD+3].template get_offset_from_indices(indices);
          adj.push_back(offset+nx);
        }//Fy(i,j,k-1) 
      } 
      else //Ez-->F
      { 
       nx = ms_.index_spaces[N][FD].template size();
       ny = ms_.index_spaces[N][FD+1].template size();
       indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent-nx-ny);
       nx = ms_.index_spaces[N][TD+2].template size();
       
       if (indices[0] <= (ms_.index_spaces[N][TD+3].template max<0>()) &&
            indices[1] <= (ms_.index_spaces[N][TD+3].template max<1>()) )
        {
          offset = ms_.index_spaces[N][TD+3].template get_offset_from_indices(indices);
          adj.push_back(offset+nx);
        }//Fy(i,j,k)
        
       if (indices[0] <= (ms_.index_spaces[N][TD+2].template max<0>()) &&
            indices[1] <= (ms_.index_spaces[N][TD+2].template max<1>()) )
        {
          offset = ms_.index_spaces[N][TD+2].template get_offset_from_indices(indices);
          adj.push_back(offset);
        } //Fx(i,j,k)
        
       if (indices[0] >= (ms_.index_spaces[N][TD+3].template min<0>()+1) &&
            indices[1] <= (ms_.index_spaces[N][TD+3].template max<1>()) )
        {
          ngb_indices = indices;
          ngb_indices[0] = ngb_indices[0] - 1;
          offset = ms_.index_spaces[N][TD+3].template get_offset_from_indices(indices);
          adj.push_back(offset+nx);
        }//Fy(i-1,j,k)
      
        if (indices[0] <= (ms_.index_spaces[N][TD+2].template max<1>()) &&
            indices[1] >= (ms_.index_spaces[N][TD+2].template min<2>()+1) )
        {
          ngb_indices = indices;
          ngb_indices[1] = ngb_indices[1] - 1;
          offset = ms_.index_spaces[N][TD+2].template get_offset_from_indices(indices);
          adj.push_back(offset);
        } //Fx(i,j-1,k)  
      }
    }
    else if ((FD == 1) && (TD == 3)) //E-->C
    {
      if (ent < ms_.index_spaces[N][FD].template size()) //Ex -->C
      {
        indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent);
      
        if ((indices[0] <= (ms_.index_spaces[N][TD+4].template max<0>())) && 
           (indices[1] <= (ms_.index_spaces[N][TD+4].template max<1>())) &&
           (indices[2] <= (ms_.index_spaces[N][TD+4].template max<2>()))) 
         {
          offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(indices);
          adj.push_back(offset);
         } //C(i,j,k)
        
        if ( (indices[0] >= (ms_.index_spaces[N][TD+4].template min<0>()+1)) && 
           (indices[1] <= (ms_.index_spaces[N][TD+4].template max<1>()))  &&
           (indices[2] <= (ms_.index_spaces[N][TD+4].template max<2>())))
         {
          ngb_indices = indices;
          ngb_indices[0] = ngb_indices[0] - 1;
          offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(ngb_indices);
          adj.push_back(offset);
         } //C(i-1,j,k)
      
        if ( (indices[0] >= (ms_.index_spaces[N][TD+4].template min<0>()+1)) && 
           (indices[1] <= (ms_.index_spaces[N][TD+4].template max<1>()))  &&
           (indices[2] >= (ms_.index_spaces[N][TD+4].template min<2>()+1)))
        {
         ngb_indices = indices;
         ngb_indices[0] = ngb_indices[0] - 1;
         ngb_indices[2] = ngb_indices[2] - 1;
         offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(ngb_indices);
         adj.push_back(offset);
        }//C(i-1,j,k-1)
      
        if ( (indices[0] <= (ms_.index_spaces[N][TD+4].template max<0>())) && 
           (indices[1] <= (ms_.index_spaces[N][TD+4].template max<1>()))  &&
           (indices[2] >= (ms_.index_spaces[N][TD+4].template min<2>()+1)))
        {
         ngb_indices = indices;
         ngb_indices[2] = ngb_indices[2] - 1;
         offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(ngb_indices);
         adj.push_back(offset);
        }//C(i,j,k-1) 
      }
      else if  ((ent >= ms_.index_spaces[N][FD].template size()) &&
               (ent <( ms_.index_spaces[N][FD+1].template size())+
               (ms_.index_spaces[N][FD].template size())) ) //Ey -->F
      {
        nx = ms_.index_spaces[N][FD].template size();
        indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent-nx);
      
        if ((indices[0] <= (ms_.index_spaces[N][TD+4].template max<0>())) && 
           (indices[1] <= (ms_.index_spaces[N][TD+4].template max<1>())) &&
           (indices[2] <= (ms_.index_spaces[N][TD+4].template max<2>()))) 
        {
         offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(indices);
         adj.push_back(offset);
        } //C(i,j,k)
        
        if ( (indices[0] <= (ms_.index_spaces[N][TD+4].template max<0>())) && 
           (indices[1] >= (ms_.index_spaces[N][TD+4].template min<1>()+1))  &&
           (indices[2] <= (ms_.index_spaces[N][TD+4].template max<2>())))
        {
         ngb_indices = indices;
         ngb_indices[1] = ngb_indices[1] - 1;
         offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(ngb_indices);
         adj.push_back(offset);
        } //C(i,j-1,k)
      
        if ( (indices[0] <= (ms_.index_spaces[N][TD+4].template max<0>())) && 
           (indices[1] >= (ms_.index_spaces[N][TD+4].template min<1>()+1))  &&
           (indices[2] >= (ms_.index_spaces[N][TD+4].template min<2>()+1)))
        {
         ngb_indices = indices;
         ngb_indices[1] = ngb_indices[1] - 1;
         ngb_indices[2] = ngb_indices[2] - 1;
         offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(ngb_indices);
         adj.push_back(offset);
        }//C(i,j-1,k-1)
      
        if ( (indices[0] <= (ms_.index_spaces[N][TD+4].template max<0>())) && 
           (indices[1] <= (ms_.index_spaces[N][TD+4].template max<1>()))  &&
           (indices[2] >= (ms_.index_spaces[N][TD+4].template min<2>()+1)))
        {
         ngb_indices = indices;
         ngb_indices[2] = ngb_indices[2] - 1;
         offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(ngb_indices);
         adj.push_back(offset);
        }//C(i,j,k-1)
      }
      else //Ez-->C
      {
        nx = ms_.index_spaces[N][FD].template size();
        ny = ms_.index_spaces[N][FD+1].template size();
        indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent-nx-ny);
      
        if ((indices[0] <= (ms_.index_spaces[N][TD+4].template max<0>())) && 
           (indices[1] <= (ms_.index_spaces[N][TD+4].template max<1>())) &&
           (indices[2] <= (ms_.index_spaces[N][TD+4].template max<2>()))) 
        {
         offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(indices);
         adj.push_back(offset);
        } //C(i,j,k)
        
        if ( (indices[0] >= (ms_.index_spaces[N][TD+4].template min<0>()+1)) && 
           (indices[1] <= (ms_.index_spaces[N][TD+4].template max<1>()))  &&
           (indices[2] <= (ms_.index_spaces[N][TD+4].template max<2>())))
        {
         ngb_indices = indices;
         ngb_indices[0] = ngb_indices[0] - 1;
         offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(ngb_indices);
         adj.push_back(offset);
        } //C(i-1,j,k)
      
        if ( (indices[0] >= (ms_.index_spaces[N][TD+4].template min<0>()+1)) && 
           (indices[1] >= (ms_.index_spaces[N][TD+4].template min<1>()+1))  &&
           (indices[2] <= (ms_.index_spaces[N][TD+4].template max<2>())))
        {
         ngb_indices = indices;
         ngb_indices[0] = ngb_indices[0] - 1;
         ngb_indices[1] = ngb_indices[1] - 1;
         offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(ngb_indices);
         adj.push_back(offset);
        }//C(i-1,j-1,k)
      
        if ( (indices[0] <= (ms_.index_spaces[N][TD+4].template max<0>())) && 
           (indices[1] >= (ms_.index_spaces[N][TD+4].template min<1>()+1))  &&
           (indices[2] <= (ms_.index_spaces[N][TD+4].template max<2>())))
        {
         ngb_indices = indices;
         ngb_indices[1] = ngb_indices[1] - 1;
         offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(ngb_indices);
         adj.push_back(offset);
        }//C(i,j-1,k)
      }
    }
    else if ((FD == 2) && (TD == 3)) //F-->C
    {
      if (ent < ms_.index_spaces[N][FD+2].template size()) //Fx -->C
      {
        indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent);
      
        if ((indices[0] <= (ms_.index_spaces[N][TD+4].template max<0>())) && 
           (indices[1] <= (ms_.index_spaces[N][TD+4].template max<1>())) &&
           (indices[2] <= (ms_.index_spaces[N][TD+4].template max<2>()))) 
        {
         offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(indices);
         adj.push_back(offset);
        } //C(i,j,k)
        
        if ( (indices[0] >= (ms_.index_spaces[N][TD+4].template min<0>()+1)) && 
           (indices[1] <= (ms_.index_spaces[N][TD+4].template max<1>()))  &&
           (indices[2] <= (ms_.index_spaces[N][TD+4].template max<2>())))
        {
         ngb_indices = indices;
         ngb_indices[0] = ngb_indices[0] - 1;
         offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(ngb_indices);
         adj.push_back(offset);
        } //C(i-1,j,k)
      }
      else if  ((ent >= ms_.index_spaces[N][FD+2].template size()) && 
               (ent <( ms_.index_spaces[N][FD+3].template size())+
               (ms_.index_spaces[N][FD+2].template size())) ) //Fy -->F
      {
        nx = ms_.index_spaces[N][FD+2].template size();
        indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent-nx);
      
        if ((indices[0] <= (ms_.index_spaces[N][TD+4].template max<0>())) && 
           (indices[1] <= (ms_.index_spaces[N][TD+4].template max<1>())) &&
           (indices[2] <= (ms_.index_spaces[N][TD+4].template max<2>()))) 
        {
         offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(indices);
         adj.push_back(offset);
        } //C(i,j,k)
        
         if ( (indices[0] <= (ms_.index_spaces[N][TD+4].template max<0>())) && 
           (indices[1] >= (ms_.index_spaces[N][TD+4].template min<1>()+1))  &&
           (indices[2] <= (ms_.index_spaces[N][TD+4].template max<2>())))
        {
         ngb_indices = indices;
         ngb_indices[1] = ngb_indices[1] - 1;
         offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(ngb_indices);
         adj.push_back(offset);
        }//C(i,j-1,k)
      }
      else //Fz-->C
      {
        nx = ms_.index_spaces[N][FD+2].template size();
        ny = ms_.index_spaces[N][FD+3].template size();
        indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent-nx-ny);
      
        if ((indices[0] <= (ms_.index_spaces[N][TD+4].template max<0>())) && 
           (indices[1] <= (ms_.index_spaces[N][TD+4].template max<1>())) &&
           (indices[2] <= (ms_.index_spaces[N][TD+4].template max<2>()))) 
        {
         offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(indices);
         adj.push_back(offset);
        } //C(i,j,k)
        
         if ( (indices[0] <= (ms_.index_spaces[N][TD+4].template max<0>())) && 
           (indices[1] <= (ms_.index_spaces[N][TD+4].template max<1>()))  &&
           (indices[2] >= (ms_.index_spaces[N][TD+4].template min<2>()+1)))
        {
         ngb_indices = indices;
         ngb_indices[2] = ngb_indices[2] - 1;
         offset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(ngb_indices);
         adj.push_back(offset);
        }//C(i,j,k-1)
      }
    } 
    else 
      std::cerr<<"requesting wrong adjacency query";
  
  return adj;

  }//entities_up_3D_

  
  template<size_t FD, size_t TD, size_t N>
  auto entities_down_3D_(size_t ent)
  {
    id_vector_t adj, indices;
    size_t index = get_index_in_storage(ent, FD, N);  
    size_t offset, xoffset, yoffset, zoffset, nx, ny;  

    if ((FD == 1) && (TD == 0)) //E-->V
      {
        if (ent < ms_.index_spaces[N][FD].template size()) //Ex-->V
         {
            indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent);
            offset = ms_.index_spaces[N][TD].template get_offset_from_indices(indices);
            adj.push_back(offset); //V(i,j,k)
            adj.push_back(entities<TD,N,0,1,0>(offset));//V(i,j+1,k)
         }
       else if ((ent >= ms_.index_spaces[N][FD].template size()) &&
               (ent <( ms_.index_spaces[N][FD+1].template size())+
               (ms_.index_spaces[N][FD].template size())) ) //Ey -->V
         {
            nx = ms_.index_spaces[N][FD].template size();
            indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent-nx);
            offset = ms_.index_spaces[N][TD].template get_offset_from_indices(indices);   
            adj.push_back(offset);//V(i,j,k)
            adj.push_back(entities<TD,N,1,0,0>(offset));//V(i+1,j,k)
         }
       else //Ez-->V
        {
            nx = ms_.index_spaces[N][FD].template size();
            ny = ms_.index_spaces[N][FD+1].template size();
            indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent-nx-ny);
            offset = ms_.index_spaces[N][TD].template get_offset_from_indices(indices);   
            adj.push_back(offset);//V(i,j,k)
            adj.push_back(entities<TD,N,0,0,1>(offset));//V(i,j,k+1)
        
        }   
      }
    else if ((FD == 2) && (TD == 0)) //F-->V
    {
      if (ent < ms_.index_spaces[N][FD+2].template size()) //Fx-->V
      {
        indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent);
        offset = ms_.index_spaces[N][TD].template get_offset_from_indices(indices);
        adj.push_back(offset); //V(i,j,k)
        adj.push_back(entities<TD,N,0,1,0>(offset));//V(i,j+1,k)
        adj.push_back(entities<TD,N,0,1,1>(offset));//V(i,j+1,k+1)
        adj.push_back(entities<TD,N,0,0,1>(offset));//V(i,j,k+1)
      }
      else if ((ent >= ms_.index_spaces[N][FD+2].template size()) && 
              (ent <( ms_.index_spaces[N][FD+3].template size())+
              (ms_.index_spaces[N][FD+2].template size())) ) //Fy -->V
      {
        nx = ms_.index_spaces[N][FD+2].template size();
        indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent-nx);
        offset = ms_.index_spaces[N][TD].template get_offset_from_indices(indices);
        adj.push_back(offset);//V(i,j,k)
        adj.push_back(entities<TD,N,1,0,0>(offset));//V(i+1,j,k)
        adj.push_back(entities<TD,N,1,0,1>(offset));//V(i+1,j,k+1)
        adj.push_back(entities<TD,N,0,0,1>(offset));//V(i,j,k+1)
      }
      else//Fz-->V 
      {
        nx = ms_.index_spaces[N][FD+2].template size();
        ny = ms_.index_spaces[N][FD+3].template size();
        indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent-nx-ny);
        offset = ms_.index_spaces[N][TD].template get_offset_from_indices(indices);
        adj.push_back(offset);//V(i,j,k)
        adj.push_back(entities<TD,N,1,0,0>(offset));//V(i+1,j,k)
        adj.push_back(entities<TD,N,1,1,0>(offset));//V(i+1,j+1,k)
        adj.push_back(entities<TD,N,0,1,0>(offset));//V(i,j+1,k)
      }
    }   
    else if ((FD == 2) && (TD == 1)) //F-->E
    {
      if (ent < ms_.index_spaces[N][FD+2].template size()) //Fx-->E
      {
        indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent);
        nx = ms_.index_spaces[N][TD].template size();
        ny = ms_.index_spaces[N][TD+1].template size();
        zoffset = ms_.index_spaces[N][TD+2].template get_offset_from_indices(indices);
        xoffset = ms_.index_spaces[N][TD].template get_offset_from_indices(indices);  
      
        adj.push_back(xoffset);//Ex(i,j,k)
        adj.push_back(entities<TD,N,0,1,0>(zoffset+nx+ny));//Ez(i,j+1,k) 
        adj.push_back(entities<TD,N,0,0,1>(xoffset));//Ex(i,jk+1)     
        adj.push_back(zoffset+nx+ny);//Ez(i,j,k)
      }
      else if ((ent >= ms_.index_spaces[N][FD+2].template size()) && 
              (ent <( ms_.index_spaces[N][FD+3].template size())+
              (ms_.index_spaces[N][FD+2].template size())) ) //Fy -->E
      {
        nx = ms_.index_spaces[N][FD+2].template size();
        indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent-nx);
        nx = ms_.index_spaces[N][TD].template size();
        ny = ms_.index_spaces[N][TD+1].template size();
        yoffset = ms_.index_spaces[N][TD+1].template get_offset_from_indices(indices);
        zoffset = ms_.index_spaces[N][TD+2].template get_offset_from_indices(indices);  
      
        adj.push_back(yoffset+nx);//Ey(i,j,k)
        adj.push_back(entities<TD,N,1,0,0>(zoffset+nx+ny));//Ez(i+1,j,k) 
        adj.push_back(entities<TD,N,0,0,1>(yoffset+nx));//Ey(i,j,k+1)     
        adj.push_back(zoffset+nx+ny);//Ez(i,j,k)
      }
      else //Fz-->E
      {
        nx = ms_.index_spaces[N][FD+2].template size();
        ny = ms_.index_spaces[N][FD+3].template size();
        indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent-nx-ny);
        nx = ms_.index_spaces[N][TD].template size();
        ny = ms_.index_spaces[N][TD+1].template size();
        xoffset = ms_.index_spaces[N][TD].template get_offset_from_indices(indices);
        yoffset = ms_.index_spaces[N][TD+1].template get_offset_from_indices(indices);  
      
        adj.push_back(yoffset+nx);//Ey(i,j,k)
        adj.push_back(entities<TD,N,1,0,0>(xoffset));//Ex(i+1,j,k) 
        adj.push_back(entities<TD,N,0,1,0>(yoffset+nx));//Ey(i,j+1,k)     
        adj.push_back(xoffset);//Ex(i,j,k)
      }
    }   
    else if ((FD == 3) && (TD == 0)) //C-->V   
    {
      indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent);
      offset = ms_.index_spaces[N][TD].template get_offset_from_indices(indices);
      adj.push_back(offset);//C(i,j,k)
      adj.push_back(entities<TD,N,1,0,0>(offset));//C(i+1,j,k)
      adj.push_back(entities<TD,N,1,1,0>(offset));//C(i+1,j+1,k)
      adj.push_back(entities<TD,N,0,1,0>(offset));//C(i,j+1,k)
      adj.push_back(entities<TD,N,0,0,1>(offset));//C(i,j,k+1)
      adj.push_back(entities<TD,N,1,0,1>(offset));//C(i+1,j,k+1)
      adj.push_back(entities<TD,N,1,1,1>(offset));//C(i+1,j+1,k+1)
      adj.push_back(entities<TD,N,0,1,1>(offset));//C(i,j+1,k+1)
 
    }
   else if ((FD == 3) && (TD == 1)) //C-->E
   {
      indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent);
      nx = ms_.index_spaces[N][TD].template size();
      ny = ms_.index_spaces[N][TD+1].template size();
      xoffset = ms_.index_spaces[N][TD].template get_offset_from_indices(indices);
      yoffset = ms_.index_spaces[N][TD+1].template get_offset_from_indices(indices);
      zoffset = ms_.index_spaces[N][TD+2].template get_offset_from_indices(indices);  
      
      adj.push_back(yoffset+nx); //Ey(i,j,k)
      adj.push_back(entities<TD,N,1,0,0>(xoffset));//Ex(i+1,j,k) 
      adj.push_back(entities<TD,N,0,1,0>(yoffset+nx)); //Ey(i,j+1,k)
      adj.push_back(xoffset);//Ex(i,j,k)
     
      adj.push_back(zoffset+nx+ny); //Ez(i,j,k)
      adj.push_back(entities<TD,N,1,0,0>(zoffset+nx+ny)); //Ez(i+1,j,k)
      adj.push_back(entities<TD,N,1,1,0>(zoffset+nx+ny)); //Ez(i+1,j+1,k)    
      adj.push_back(entities<TD,N,0,1,0>(zoffset+nx+ny)); //Ez(i,j+1,k) 
     
      adj.push_back(entities<TD,N,0,0,1>(yoffset+nx)); //Ey(i,j,k+1)
      adj.push_back(entities<TD,N,1,0,1>(xoffset));//Ex(i+1,j,k+1) 
      adj.push_back(entities<TD,N,0,1,1>(yoffset+nx)); //Ey(i,j+1,k+1)
      adj.push_back(entities<TD,N,0,0,1>(xoffset));//Ex(i,j,k+1)     
   }
   else if ((FD == 3) && (TD == 2)) //C-->F
   {
      indices = ms_.index_spaces[N][index].template get_indices_from_offset(ent);
      nx = ms_.index_spaces[N][TD+2].template size();
      ny = ms_.index_spaces[N][TD+3].template size();
      xoffset = ms_.index_spaces[N][TD+2].template get_offset_from_indices(indices);
      yoffset = ms_.index_spaces[N][TD+3].template get_offset_from_indices(indices);
      zoffset = ms_.index_spaces[N][TD+4].template get_offset_from_indices(indices);  
      
      adj.push_back(yoffset+nx); //Fy(i,j,k)
      adj.push_back(entities<TD,N,1,0,0>(xoffset));//Fx(i+1,j,k) 
      adj.push_back(entities<TD,N,0,1,0>(yoffset+nx)); //Fy(i,j+1,k)
      adj.push_back(xoffset);//Fx(i,j,k)
      adj.push_back(zoffset+nx+ny); //Fz(i,j,k)
      adj.push_back(entities<TD,N,0,0,1>(zoffset+nx+ny)); //Fz(i,j,k+1)
   } 
  else 
    std::cerr<<"requesting wrong adjacency query";
  
  return adj;
      
  } //entities_down_3D_
  

  //Utility methods
  auto get_bnds(size_t mdim, size_t entdim)
  {
   id_vector_t bnds; 
  
   if (mdim == 1){
    switch(entdim){
      case 0:   
         bnds.push_back(0);
         break;
      case 1: 
         bnds.push_back(1);
         break;
      default:
        std::cerr<<"Requesting non-valid entity bounds for dimension 1";
    }
   }
   else if (mdim == 2){
    switch(entdim){
      case 0:   
         bnds.push_back(0);
         bnds.push_back(0);
         break;
      case 1:
         bnds.push_back(0);
         bnds.push_back(1);
         break;
      case 2:
         bnds.push_back(1);
         bnds.push_back(0);
         break;
      case 3:
         bnds.push_back(1);
         bnds.push_back(1);
         break;
      default:
        std::cerr<<"Requesting non-valid entity bounds for dimension 2";
    }
   }
   else if (mdim == 3){
    switch(entdim){
      case 0:
         bnds.push_back(0);
         bnds.push_back(0);
         bnds.push_back(0);
         break;
      case 1:
         bnds.push_back(0);
         bnds.push_back(1);
         bnds.push_back(0);
         break;
      case 2:
         bnds.push_back(1);
         bnds.push_back(0);
         bnds.push_back(0);
         break;
      case 3:
         bnds.push_back(0);
         bnds.push_back(0);
         bnds.push_back(1);
         break;
      case 4:
         bnds.push_back(0);
         bnds.push_back(1);
         bnds.push_back(1);
         break;
      case 5:
         bnds.push_back(1);
         bnds.push_back(0);
         bnds.push_back(1);
         break;
      case 6:
         bnds.push_back(1);
         bnds.push_back(1);
         bnds.push_back(0);
         break;
      case 7:
         bnds.push_back(1);
         bnds.push_back(1);
         bnds.push_back(1);
         break;
      default:
        std::cerr<<"Requesting non-valid entity bounds for dimension 3";
     }
   }
   else 
    std::cerr<<"Requesting entity bounds for dimension greater than 3";
      
   return bnds;
  }//get_bnds

  auto get_index_in_storage(size_t ent, size_t dim, size_t dom)
  {
    size_t index;
    if (meshdim_ == 1)
    {
      switch(dim)
      {
        case 0:
          index = 0;
          break;
        case 1: 
          index = 1;
          break;
        default:
          std::cerr<<"non-valid index request\n";  
       }
    }
    else if (meshdim_ == 2)
    {
      switch(dim)
      {
        case 0:
          index = 0;
          break;
        case 1: 
          index = 1 + (ent >= ms_.index_spaces[dom][dim].template size());
          break;
        case 2:
          index = 3;
          break;
        default:
          std::cerr<<"non-valid index request";  
      }
    }
    else if (meshdim_ == 3)
    {
      switch(dim)
      {
        case 0:
          index = 0;
          break;
        case 1: 
          if (ent < ms_.index_spaces[dom][dim].template size())
              index = 1;
          else if ((ent < (ms_.index_spaces[dom][dim+1].template size())+
                  (ms_.index_spaces[dom][dim].template size())) && 
                  (ent >= ms_.index_spaces[dom][dim].template size()))
              index = 2; 
          else if ((ent < (ms_.index_spaces[dom][dim+2].template size())+
                  (ms_.index_spaces[dom][dim+1].template size())+
                  (ms_.index_spaces[dom][dim].template size())) && 
                  (ent >= (ms_.index_spaces[dom][dim+1].template size())+
                  (ms_.index_spaces[dom][dim].template size()))) 
             index = 3;
          else 
            std::cerr<<"non-valid index request";
          break;
        case 2:
          if (ent < ms_.index_spaces[dom][dim+2].template size())
              index = 4;
          else if ((ent < (ms_.index_spaces[dom][dim+3].template size())+
                  (ms_.index_spaces[dom][dim+2].template size())) && 
                  (ent >= ms_.index_spaces[dom][dim+2].template size()))
              index = 5; 
          else if ((ent < (ms_.index_spaces[dom][dim+4].template size())+
                  (ms_.index_spaces[dom][dim+3].template size())+
                  (ms_.index_spaces[dom+2][dim].template size())) && 
                  (ent >= (ms_.index_spaces[dom][dim+3].template size())+
                  (ms_.index_spaces[dom][dim+2].template size()))) 
             index = 6;
          else 
            std::cerr<<"non-valid index request";
          break;
        case 3:
          index = 7;
          break;
        default:
          std::cerr<<"non-valid index request";  
      }
    }
    else
      std::cerr<<"Requesting index location for dim greater than 3\n";
   return index;    
  }

}; // class structured_mesh_topology_t

} // namespace topology
} // namespace flecsi

#endif // flecsi_structured_mesh_topology_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
