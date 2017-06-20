#include <iostream>
#include <cinchtest.h>
#include "flecsi/topology/structured_mesh_topology.h"

using namespace std;
using namespace flecsi;
using namespace topology;

class Vertex: public structured_mesh_entity_t<0,1>{
public:
  Vertex(){}

  Vertex(structured_mesh_topology_base_t &){}
};

class Cell : public structured_mesh_entity_t<3, 1>{
public:

  Cell(){}

  Cell(structured_mesh_topology_base_t &){}
};


class TestMesh3dType{
public:
  static constexpr size_t num_dimensions = 3;
  static constexpr size_t num_domains = 1;

  static constexpr std::array<size_t,num_dimensions> lower_bounds = {0,0,0};
  static constexpr std::array<size_t,num_dimensions> upper_bounds = {3,2,2};

  using entity_types = std::tuple<
  std::pair<domain_<0>, Vertex>,
  std::pair<domain_<0>, Cell>>;

};

constexpr std::array<size_t,TestMesh3dType::num_dimensions> TestMesh3dType::lower_bounds;
constexpr std::array<size_t,TestMesh3dType::num_dimensions> TestMesh3dType::upper_bounds;

using id_vector_t = std::vector<size_t>;
using TestMesh = structured_mesh_topology_t<TestMesh3dType>;

TEST(structured3d, simple){

  auto mesh = new TestMesh;
 
  size_t nv, nex, ney, nez, nfx, nfy, nfz, nc;
  id_vector_t adj;

  auto lbnd = mesh->lower_bounds();
  auto ubnd = mesh->upper_bounds();

  CINCH_CAPTURE() << "3D Logically structured mesh with bounds: [" <<lbnd[0]<<
  ", "<<lbnd[1]<<", "<<lbnd[2]<<"] - ["<<ubnd[0]<<", "<<ubnd[1]<<", "<<ubnd[2]<<"] \n"<< endl;

  
  nv  = mesh->num_entities(0,0);
  nex = mesh->num_entities(1,0);
  ney = mesh->num_entities(2,0);
  nez = mesh->num_entities(3,0);	
  nfx = mesh->num_entities(4,0);
  nfy = mesh->num_entities(5,0);
  nfz = mesh->num_entities(6,0);
  nc  = mesh->num_entities(7,0);
  
  CINCH_CAPTURE() << "NV   = "<< nv  << endl;
  CINCH_CAPTURE() << "NE_X = "<< nex << endl;
  CINCH_CAPTURE() << "NE_Y = "<< ney << endl;
  CINCH_CAPTURE() << "NE_Z = "<< nez << endl;
  CINCH_CAPTURE() << "NF_X = "<< nfx << endl;
  CINCH_CAPTURE() << "NF_Y = "<< nfy << endl;
  CINCH_CAPTURE() << "NF_Z = "<< nfz << endl;
  CINCH_CAPTURE() << "NC   = "<< nc  << endl;
  CINCH_CAPTURE()<<endl;
 
  //Loop over all vertices and test intra index space queries
  CINCH_CAPTURE()<<"------Vertices------"<<endl;
  for (auto vertex: mesh->entities<0>()){
   CINCH_CAPTURE() << "---- vertex id: " << vertex << endl; 
   CINCH_CAPTURE() << "  -- indices "<< endl; 
   for (auto idv : mesh->get_indices<0>(vertex))
    CINCH_CAPTURE() << "  ---- " <<idv << endl;
   auto id = mesh->get_indices<0>(vertex); 
   auto offset = mesh->get_offset<0>(id);
   CINCH_CAPTURE() << "  ---- offset " <<offset<< endl;
   ASSERT_EQ(vertex,offset);   

   //V-->V
   CINCH_CAPTURE() << "  -- stencil [ 1  0  0] " <<mesh->entities< 0, 0, 1, 0, 0 >(vertex) << endl;
   CINCH_CAPTURE() << "  -- stencil [-1  0  0] " <<mesh->entities< 0, 0,-1, 0, 0 >(vertex) << endl; 
   CINCH_CAPTURE() << "  -- stencil [ 0  1  0] " <<mesh->entities< 0, 0, 0, 1, 0 >(vertex) << endl;
   CINCH_CAPTURE() << "  -- stencil [ 0 -1  0] " <<mesh->entities< 0, 0, 0,-1, 0 >(vertex) << endl;
   CINCH_CAPTURE() << "  -- stencil [ 0  0  1] " <<mesh->entities< 0, 0, 0, 0, 1 >(vertex) << endl;
   CINCH_CAPTURE() << "  -- stencil [ 0  0 -1] " <<mesh->entities< 0, 0, 0, 0,-1 >(vertex) << endl;

  //std::cout<<"v = "<<vertex<<std::endl;
  
   //V-->E
   CINCH_CAPTURE() << "  -- query V-->E "<< endl; 
   for (auto edge : mesh->get_entities<0,1,0>(vertex))
    CINCH_CAPTURE() << "  ---- " <<edge<< endl; 
   
   //V-->F
   CINCH_CAPTURE() << "  -- query V-->F "<< endl; 
   for (auto face : mesh->get_entities<0,2,0>(vertex))
    CINCH_CAPTURE() << "  ---- " <<face<< endl; 
  
   //V-->C
   CINCH_CAPTURE() << "  -- query V-->C "<< endl; 
   for (auto cell : mesh->get_entities<0,3,0>(vertex))
    CINCH_CAPTURE() << "  ---- " <<cell<< endl; 
   
   CINCH_CAPTURE()<<endl;
  }
  
  //Loop over all edges in X-direction and test intra index space queries
  CINCH_CAPTURE()<<"------Edges-X------"<<endl;
  for (auto edgex: mesh->entities<1>()){
   CINCH_CAPTURE() << "---- edgex id: " << edgex << endl;
   CINCH_CAPTURE() << "  -- indices "<< endl; 
   for (auto idv : mesh->get_indices<1>(edgex))
    CINCH_CAPTURE() << "  ---- " <<idv << endl; 
   auto id = mesh->get_indices<1>(edgex); 
   auto offset = mesh->get_offset<1>(id);
   CINCH_CAPTURE() << "  ---- offset " <<offset<< endl;
   ASSERT_EQ(edgex,offset);   

   //E-->E
   CINCH_CAPTURE() << "  -- stencil [ 1  0  0] " <<mesh->entities< 1, 0, 1, 0, 0 >(edgex) << endl;
   CINCH_CAPTURE() << "  -- stencil [-1  0  0] " <<mesh->entities< 1, 0,-1, 0, 0 >(edgex) << endl; 
   CINCH_CAPTURE() << "  -- stencil [ 0  1  0] " <<mesh->entities< 1, 0, 0, 1, 0 >(edgex) << endl;
   CINCH_CAPTURE() << "  -- stencil [ 0 -1  0] " <<mesh->entities< 1, 0, 0,-1, 0 >(edgex) << endl;
   CINCH_CAPTURE() << "  -- stencil [ 0  0  1] " <<mesh->entities< 1, 0, 0, 0, 1 >(edgex) << endl;
   CINCH_CAPTURE() << "  -- stencil [ 0  0 -1] " <<mesh->entities< 1, 0, 0, 0,-1 >(edgex) << endl;
   
   //E-->V
   CINCH_CAPTURE() << "  -- query E-->V "<< endl; 
   for (auto vert : mesh->get_entities<1,0,0>(edgex))
    CINCH_CAPTURE() << "  ---- " <<vert<< endl; 
   
   //E-->F
   CINCH_CAPTURE() << "  -- query E-->F "<< endl; 
   for (auto face : mesh->get_entities<1,2,0>(edgex))
    CINCH_CAPTURE() << "  ---- " <<face<< endl; 
   
   //E-->C
   CINCH_CAPTURE() << "  -- query E-->C "<< endl; 
   for (auto cell : mesh->get_entities<1,3,0>(edgex))
    CINCH_CAPTURE() << "  ---- " <<cell<< endl;
 
   CINCH_CAPTURE()<<endl;
  }

  
  //Loop over all edges in Y-directions and test intra index space queries
  CINCH_CAPTURE()<<"------Edges-Y------"<<endl;
  for (auto edgey: mesh->entities<2>()){
   CINCH_CAPTURE() << "---- edgey id: " << edgey << endl;
   CINCH_CAPTURE() << "  -- indices "<< endl; 
   for (auto idv : mesh->get_indices<2>(edgey))
    CINCH_CAPTURE() << "  ---- " <<idv << endl; 
   auto id = mesh->get_indices<2>(edgey); 
   auto offset = mesh->get_offset<2>(id);
   CINCH_CAPTURE() << "  ---- offset " <<offset<< endl;
   ASSERT_EQ(edgey,offset);   
  
   //E-->E
   CINCH_CAPTURE() << "  -- stencil [ 1  0  0] " <<mesh->entities< 1, 0, 1, 0, 0 >(edgey+nex) << endl;
   CINCH_CAPTURE() << "  -- stencil [-1  0  0] " <<mesh->entities< 1, 0,-1, 0, 0 >(edgey+nex) << endl; 
   CINCH_CAPTURE() << "  -- stencil [ 0  1  0] " <<mesh->entities< 1, 0, 0, 1, 0 >(edgey+nex) << endl;
   CINCH_CAPTURE() << "  -- stencil [ 0 -1  0] " <<mesh->entities< 1, 0, 0,-1, 0 >(edgey+nex) << endl;
   CINCH_CAPTURE() << "  -- stencil [ 0  0  1] " <<mesh->entities< 1, 0, 0, 0, 1 >(edgey+nex) << endl;
   CINCH_CAPTURE() << "  -- stencil [ 0  0 -1] " <<mesh->entities< 1, 0, 0, 0,-1 >(edgey+nex) << endl;
   
   //E-->V
   CINCH_CAPTURE() << "  -- query E-->V "<< endl; 
   for (auto vert : mesh->get_entities<1,0,0>(edgey+nex))
    CINCH_CAPTURE() << "  ---- " <<vert<< endl; 
   
   //E-->F
   CINCH_CAPTURE() << "  -- query E-->F "<< endl; 
   for (auto face : mesh->get_entities<1,2,0>(edgey+nex))
    CINCH_CAPTURE() << "  ---- " <<face<< endl; 
   
   //E-->C
   CINCH_CAPTURE() << "  -- query E-->C "<< endl; 
   for (auto cell : mesh->get_entities<1,3,0>(edgey+nex))
    CINCH_CAPTURE() << "  ---- " <<cell<< endl;

   CINCH_CAPTURE()<<endl;
  }
 

  //Loop over all edges in Z-directions and test intra index space queries
  CINCH_CAPTURE()<<"------Edges-Z------"<<endl;
  for (auto edgez: mesh->entities<3>()){
   CINCH_CAPTURE() << "---- edgez id: " << edgez << endl;
   CINCH_CAPTURE() << "  -- indices "<< endl; 
   for (auto idv : mesh->get_indices<3>(edgez))
    CINCH_CAPTURE() << "  ---- " <<idv << endl; 
   auto id = mesh->get_indices<3>(edgez); 
   auto offset = mesh->get_offset<3>(id);
   CINCH_CAPTURE() << "  ---- offset " <<offset<< endl;
   ASSERT_EQ(edgez,offset);   
  
   //E-->E
   CINCH_CAPTURE() << "  -- stencil [ 1  0  0] " <<mesh->entities< 1, 0, 1, 0, 0 >(edgez+nex+ney) << endl;
   CINCH_CAPTURE() << "  -- stencil [-1  0  0] " <<mesh->entities< 1, 0,-1, 0, 0 >(edgez+nex+ney) << endl; 
   CINCH_CAPTURE() << "  -- stencil [ 0  1  0] " <<mesh->entities< 1, 0, 0, 1, 0 >(edgez+nex+ney) << endl;
   CINCH_CAPTURE() << "  -- stencil [ 0 -1  0] " <<mesh->entities< 1, 0, 0,-1, 0 >(edgez+nex+ney) << endl;
   CINCH_CAPTURE() << "  -- stencil [ 0  0  1] " <<mesh->entities< 1, 0, 0, 0, 1 >(edgez+nex+ney) << endl;
   CINCH_CAPTURE() << "  -- stencil [ 0  0 -1] " <<mesh->entities< 1, 0, 0, 0,-1 >(edgez+nex+ney) << endl;
   
   //E-->V
   CINCH_CAPTURE() << "  -- query E-->V "<< endl; 
   for (auto vert : mesh->get_entities<1,0,0>(edgez+nex+ney))
    CINCH_CAPTURE() << "  ---- " <<vert<< endl; 
   
   //E-->F
   CINCH_CAPTURE() << "  -- query E-->F "<< endl; 
   for (auto face : mesh->get_entities<1,2,0>(edgez+nex+ney))
    CINCH_CAPTURE() << "  ---- " <<face<< endl; 
   
   //E-->C
   CINCH_CAPTURE() << "  -- query E-->C "<< endl; 
   for (auto cell : mesh->get_entities<1,3,0>(edgez+nex+ney))
    CINCH_CAPTURE() << "  ---- " <<cell<< endl;

   CINCH_CAPTURE()<<endl;
  }

  //Loop over all faces in X-direction and test intra index space queries
  CINCH_CAPTURE()<<"------Faces-X------"<<endl;
  for (auto facex: mesh->entities<4>()){
   CINCH_CAPTURE() << "---- facex id: " << facex << endl;
   CINCH_CAPTURE() << "  -- indices "<< endl; 
   for (auto idv : mesh->get_indices<4>(facex))
    CINCH_CAPTURE() << "  ---- " <<idv << endl; 
   auto id = mesh->get_indices<4>(facex); 
   auto offset = mesh->get_offset<4>(id);
   CINCH_CAPTURE() << "  ---- offset " <<offset<< endl;
   ASSERT_EQ(facex,offset);   
  
   //F-->F
   CINCH_CAPTURE() << "  -- stencil [ 1  0  0] " <<mesh->entities< 2, 0, 1, 0, 0 >(facex) << endl;
   CINCH_CAPTURE() << "  -- stencil [-1  0  0] " <<mesh->entities< 2, 0,-1, 0, 0 >(facex) << endl; 
   CINCH_CAPTURE() << "  -- stencil [ 0  1  0] " <<mesh->entities< 2, 0, 0, 1, 0 >(facex) << endl;
   CINCH_CAPTURE() << "  -- stencil [ 0 -1  0] " <<mesh->entities< 2, 0, 0,-1, 0 >(facex) << endl;
   CINCH_CAPTURE() << "  -- stencil [ 0  0  1] " <<mesh->entities< 2, 0, 0, 0, 1 >(facex) << endl;
   CINCH_CAPTURE() << "  -- stencil [ 0  0 -1] " <<mesh->entities< 2, 0, 0, 0,-1 >(facex) << endl;
   
   //F-->V
   CINCH_CAPTURE() << "  -- query F-->V "<< endl; 
   for (auto vert : mesh->get_entities<2,0,0>(facex))
    CINCH_CAPTURE() << "  ---- " <<vert<< endl; 
   
   //F-->E
   CINCH_CAPTURE() << "  -- query F-->E "<< endl; 
   for (auto edge : mesh->get_entities<2,1,0>(facex))
    CINCH_CAPTURE() << "  ---- " <<edge<< endl; 
   
   //F-->C
   CINCH_CAPTURE() << "  -- query F-->C "<< endl; 
   for (auto cell : mesh->get_entities<2,3,0>(facex))
    CINCH_CAPTURE() << "  ---- " <<cell<< endl;

   CINCH_CAPTURE()<<endl;
  }

  //Loop over all faces in Y-direction and test intra index space queries
  CINCH_CAPTURE()<<"------Faces-Y------"<<endl;
  for (auto facey: mesh->entities<5>()){
   CINCH_CAPTURE() << "---- facey id: " << facey << endl;
   CINCH_CAPTURE() << "  -- indices "<< endl; 
   for (auto idv : mesh->get_indices<5>(facey))
    CINCH_CAPTURE() << "  ---- " <<idv << endl; 
   auto id = mesh->get_indices<5>(facey); 
   auto offset = mesh->get_offset<5>(id);
   CINCH_CAPTURE() << "  ---- offset " <<offset<< endl;
   ASSERT_EQ(facey,offset);   
  
   //F-->F
   CINCH_CAPTURE() << "  -- stencil [ 1  0  0] " <<mesh->entities< 2, 0, 1, 0, 0 >(facey+nfx) << endl;
   CINCH_CAPTURE() << "  -- stencil [-1  0  0] " <<mesh->entities< 2, 0,-1, 0, 0 >(facey+nfx) << endl; 
   CINCH_CAPTURE() << "  -- stencil [ 0  1  0] " <<mesh->entities< 2, 0, 0, 1, 0 >(facey+nfx) << endl;
   CINCH_CAPTURE() << "  -- stencil [ 0 -1  0] " <<mesh->entities< 2, 0, 0,-1, 0 >(facey+nfx) << endl;
   CINCH_CAPTURE() << "  -- stencil [ 0  0  1] " <<mesh->entities< 2, 0, 0, 0, 1 >(facey+nfx) << endl;
   CINCH_CAPTURE() << "  -- stencil [ 0  0 -1] " <<mesh->entities< 2, 0, 0, 0,-1 >(facey+nfx) << endl;
   
   //std::cout<<"fy = "<<facey<<std::endl;
   //F-->V
   CINCH_CAPTURE() << "  -- query F-->V "<< endl; 
   for (auto vert : mesh->get_entities<2,0,0>(facey+nfx))
    CINCH_CAPTURE() << "  ---- " <<vert<< endl; 
   
   //std::cout<<"fy = "<<facey<<std::endl;
   //F-->E
   CINCH_CAPTURE() << "  -- query F-->E "<< endl; 
   for (auto edge : mesh->get_entities<2,1,0>(facey+nfx))
    CINCH_CAPTURE() << "  ---- " <<edge<< endl; 
   
//   std::cout<<"fy = "<<facey<<std::endl;
   //F-->C
   CINCH_CAPTURE() << "  -- query F-->C "<< endl; 
   for (auto cell : mesh->get_entities<2,3,0>(facey+nfx))
    CINCH_CAPTURE() << "  ---- " <<cell<< endl;

   CINCH_CAPTURE()<<endl;
  }
  
  //Loop over all faces in Z-direction and test intra index space queries
  CINCH_CAPTURE()<<"------Faces-Z------"<<endl;
  for (auto facez: mesh->entities<6>()){
   CINCH_CAPTURE() << "---- facez id: " << facez << endl;
   CINCH_CAPTURE() << "  -- indices "<< endl; 
   for (auto idv : mesh->get_indices<6>(facez))
    CINCH_CAPTURE() << "  ---- " <<idv << endl; 
   auto id = mesh->get_indices<6>(facez); 
   auto offset = mesh->get_offset<6>(id);
   CINCH_CAPTURE() << "  ---- offset " <<offset<< endl;
   ASSERT_EQ(facez,offset);   
  
   //F-->F
   CINCH_CAPTURE() << "  -- stencil [ 1  0  0] " <<mesh->entities< 2, 0, 1, 0, 0 >(facez+nfx+nfy) << endl;
   CINCH_CAPTURE() << "  -- stencil [-1  0  0] " <<mesh->entities< 2, 0,-1, 0, 0 >(facez+nfx+nfy) << endl; 
   CINCH_CAPTURE() << "  -- stencil [ 0  1  0] " <<mesh->entities< 2, 0, 0, 1, 0 >(facez+nfx+nfy) << endl;
   CINCH_CAPTURE() << "  -- stencil [ 0 -1  0] " <<mesh->entities< 2, 0, 0,-1, 0 >(facez+nfx+nfy) << endl;
   CINCH_CAPTURE() << "  -- stencil [ 0  0  1] " <<mesh->entities< 2, 0, 0, 0, 1 >(facez+nfx+nfy) << endl;
   CINCH_CAPTURE() << "  -- stencil [ 0  0 -1] " <<mesh->entities< 2, 0, 0, 0,-1 >(facez+nfx+nfy) << endl;
  
 
   //std::cout<<"fz = "<<facez<<std::endl;
   //F-->V
   CINCH_CAPTURE() << "  -- query F-->V "<< endl; 
   for (auto vert : mesh->get_entities<2,0,0>(facez+nfx+nfy))
    CINCH_CAPTURE() << "  ---- " <<vert<< endl; 
   
   //F-->E
   CINCH_CAPTURE() << "  -- query F-->E "<< endl; 
   for (auto edge : mesh->get_entities<2,1,0>(facez+nfx+nfy))
    CINCH_CAPTURE() << "  ---- " <<edge<< endl; 
   
   //F-->C
   CINCH_CAPTURE() << "  -- query F-->C "<< endl; 
   for (auto cell : mesh->get_entities<2,3,0>(facez+nfx+nfy))
    CINCH_CAPTURE() << "  ---- " <<cell<< endl;

   CINCH_CAPTURE()<<endl;
  }
  
  //Loop over all cells and test intra index space queries
  CINCH_CAPTURE()<<"------Cells------"<<endl;
  for (auto cell: mesh->entities<7>()){
   CINCH_CAPTURE() << "---- cell id: " << cell << endl;
   CINCH_CAPTURE() << "  -- indices "<< endl; 
   for (auto idv : mesh->get_indices<7>(cell))
    CINCH_CAPTURE() << "  ---- " <<idv << endl; 
   auto id = mesh->get_indices<7>(cell); 
   auto offset = mesh->get_offset<7>(id);
   CINCH_CAPTURE() << "  ---- offset " <<offset<< endl;
   ASSERT_EQ(cell,offset);   
  
   //C-->C
   CINCH_CAPTURE() << "  -- stencil [ 1  0  0] " <<mesh->entities< 3, 0, 1, 0, 0 >(cell) << endl;
   CINCH_CAPTURE() << "  -- stencil [-1  0  0] " <<mesh->entities< 3, 0,-1, 0, 0 >(cell) << endl; 
   CINCH_CAPTURE() << "  -- stencil [ 0  1  0] " <<mesh->entities< 3, 0, 0, 1, 0 >(cell) << endl;
   CINCH_CAPTURE() << "  -- stencil [ 0 -1  0] " <<mesh->entities< 3, 0, 0,-1, 0 >(cell) << endl;
   CINCH_CAPTURE() << "  -- stencil [ 0  0  1] " <<mesh->entities< 3, 0, 0, 0, 1 >(cell) << endl;
   CINCH_CAPTURE() << "  -- stencil [ 0  0 -1] " <<mesh->entities< 3, 0, 0, 0,-1 >(cell) << endl;
   
   //C-->V
   CINCH_CAPTURE() << "  -- query C-->V "<< endl; 
   for (auto vert : mesh->get_entities<3,0,0>(cell))
    CINCH_CAPTURE() << "  ---- " <<vert<< endl; 
   
   //C-->E
   CINCH_CAPTURE() << "  -- query C-->E "<< endl; 
   for (auto edg : mesh->get_entities<3,1,0>(cell))
    CINCH_CAPTURE() << "  ---- " <<edg<< endl; 
   
   //C-->F
   CINCH_CAPTURE() << "  -- query C-->F "<< endl; 
   for (auto face : mesh->get_entities<3,2,0>(cell))
    CINCH_CAPTURE() << "  ---- " <<face<< endl;

   CINCH_CAPTURE()<<endl;
  }
  

  ASSERT_TRUE(CINCH_EQUAL_BLESSED("structured3d.blessed"));
  //CINCH_WRITE("structured3d.blessed");

} // TEST
