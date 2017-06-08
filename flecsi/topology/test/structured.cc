#include <iostream>
#include <cinchtest.h>
#include "flecsi/topology/structured_mesh_topology.h"

using namespace std;
using namespace flecsi;
using namespace topology;


class Vertex : public structured_mesh_entity_t<0, 1>{
public:
  Vertex(){}

  Vertex(structured_mesh_topology_base_t &){}
};

class Face : public structured_mesh_entity_t<2, 1>{
public:

  Face(){}

  Face(structured_mesh_topology_base_t &){}
};


class TestMesh2dType{
public:
  static constexpr size_t num_dimensions = 2;
  static constexpr size_t num_domains = 1;

  static constexpr std::array<size_t,num_dimensions> lower_bounds = {0,0};
  static constexpr std::array<size_t,num_dimensions> upper_bounds = {3,2};

  using entity_types = std::tuple<
  std::pair<domain_<0>, Vertex>,
  std::pair<domain_<0>, Face>>;

};

constexpr std::array<size_t,TestMesh2dType::num_dimensions> TestMesh2dType::lower_bounds;
constexpr std::array<size_t,TestMesh2dType::num_dimensions> TestMesh2dType::upper_bounds;

using id_vector_t = std::vector<size_t>;
using TestMesh = structured_mesh_topology_t<TestMesh2dType>;

TEST(structured, simple){

  auto mesh = new TestMesh; 
  size_t nv, nex, ney, nf;
  id_vector_t adj;

  nv = mesh->num_entities(0,0);
  nex = mesh->num_entities(1,0);
  ney = mesh->num_entities(2,0);
  nf = mesh->num_entities(3,0);
  
  CINCH_CAPTURE() << "NV = " << nv << endl;
  CINCH_CAPTURE() << "NE_X = "<< nex << endl;
  CINCH_CAPTURE() << "NE_Y = "<< ney << endl;
  CINCH_CAPTURE() << "NF = " << nf << endl;
 
  //Loop over all vertices and test intra index space queries
  for (auto vertex: mesh->entities<0>()){
   CINCH_CAPTURE() << "---- vertex id: " << vertex << endl;

   //V-->V
   CINCH_CAPTURE() << "  -- stencil [1 0] " <<mesh->entities<0,0,1,0>(vertex) << endl; std::cout<<std::endl;
   CINCH_CAPTURE() << "  -- stencil [0 1] " <<mesh->entities<0,0,0,1>(vertex) << endl;
   CINCH_CAPTURE() << "  -- stencil [-1 0] " <<mesh->entities<0,0,-1,0>(vertex) << endl; 
   CINCH_CAPTURE() << "  -- stencil [0 -1] " <<mesh->entities<0,0,0,-1>(vertex) << endl;
  
   //V-->E
   adj.clear();
   adj = mesh->get_entities<0,1,0>(vertex);

   CINCH_CAPTURE() << "  -- query V-->E "<< endl; 
   for (size_t k=0; k< adj.size(); ++k)
    CINCH_CAPTURE() << "  ---- " <<adj[k] << endl; 
   
   //V-->F
   adj.clear();
   adj = mesh->get_entities<0,2,0>(vertex);

   CINCH_CAPTURE() << "  -- query V-->F "<< endl; 
   for (size_t k=0; k< adj.size(); ++k)
    CINCH_CAPTURE() << "  ---- " <<adj[k] << endl; 
  }

  
  //Loop over all edges in X-direction and test intra index space queries
  for (auto edgex: mesh->entities<1>()){
   CINCH_CAPTURE() << "---- edgex id: " << edgex << endl;

   //E-->E
   CINCH_CAPTURE() << "  -- stencil [1 0] " <<mesh->entities<1,0,1,0>(edgex) << endl;
   CINCH_CAPTURE() << "  -- stencil [0 1] " <<mesh->entities<1,0,0,1>(edgex) << endl;
   CINCH_CAPTURE() << "  -- stencil [-1 0] " <<mesh->entities<1,0,-1,0>(edgex) << endl; 
   CINCH_CAPTURE() << "  -- stencil [0 -1] " <<mesh->entities<1,0,0,-1>(edgex) << endl;
   
   //E-->V
   adj.clear();
   adj = mesh->get_entities<1,0,0>(edgex);

   CINCH_CAPTURE() << "  -- query E-->V "<< endl; 
   for (size_t k=0; k< adj.size(); ++k)
    CINCH_CAPTURE() << "  ---- " <<adj[k] << endl; 
   
   //E-->F
   adj.clear();
   adj = mesh->get_entities<1,2,0>(edgex);

   CINCH_CAPTURE() << "  -- query E-->F "<< endl; 
   for (size_t k=0; k< adj.size(); ++k)
    CINCH_CAPTURE() << "  ---- " <<adj[k] << endl; 
  }

  std::cout<<std::endl;
  std::cout<<"Starting edge y"<<std::endl;
  //Loop over all edges in Y-directions and test intra index space queries
  for (auto edgey: mesh->entities<2>()){
   CINCH_CAPTURE() << "---- edgey id: " << edgey + nex << endl;
  
   //E-->E
   CINCH_CAPTURE() << "  -- stencil [1 0] " <<mesh->entities<1,0,1,0>(edgey+nex) << endl;
   CINCH_CAPTURE() << "  -- stencil [0 1] " <<mesh->entities<1,0,0,1>(edgey+nex) << endl;
   CINCH_CAPTURE() << "  -- stencil [-1 0] " <<mesh->entities<1,0,-1,0>(edgey+nex) << endl; 
   CINCH_CAPTURE() << "  -- stencil [0 -1] " <<mesh->entities<1,0,0,-1>(edgey+nex) << endl;
   
   //E-->V
   adj.clear();
   adj = mesh->get_entities<1,0,0>(edgey+nex);

   CINCH_CAPTURE() << "  -- query E-->V "<< endl; 
   for (size_t k=0; k< adj.size(); ++k)
    CINCH_CAPTURE() << "  ---- " <<adj[k] << endl; 
   
   //E-->F
   adj.clear();
   adj = mesh->get_entities<1,2,0>(edgey+nex);

   CINCH_CAPTURE() << "  -- query E-->F "<< endl; 
   for (size_t k=0; k< adj.size(); ++k)
    CINCH_CAPTURE() << "  ---- " <<adj[k] << endl; 
  }

   //Loop over all face in Y-directions and test intra index space queries
  for (auto face: mesh->entities<3>()){
   CINCH_CAPTURE() << "---- face id: " << face << endl;
  
   //F-->F
   CINCH_CAPTURE() << "  -- stencil [1 0] " <<mesh->entities<2,0,1,0>(face) << endl;
   CINCH_CAPTURE() << "  -- stencil [0 1] " <<mesh->entities<2,0,0,1>(face) << endl;
   CINCH_CAPTURE() << "  -- stencil [-1 0] " <<mesh->entities<2,0,-1,0>(face) << endl; 
   CINCH_CAPTURE() << "  -- stencil [0 -1] " <<mesh->entities<2,0,0,-1>(face) << endl;
   
   //F-->V
   adj.clear();
   adj = mesh->get_entities<2,0,0>(face);

   CINCH_CAPTURE() << "  -- query F-->V "<< endl; 
   for (size_t k=0; k< adj.size(); ++k)
    CINCH_CAPTURE() << "  ---- " <<adj[k] << endl; 
   
   //F-->E
   adj.clear();
   adj = mesh->get_entities<2,1,0>(face);

   CINCH_CAPTURE() << "  -- query F-->E "<< endl; 
   for (size_t k=0; k< adj.size(); ++k)
    CINCH_CAPTURE() << "  ---- " <<adj[k] << endl; 
  }


/*  for(auto cell : mesh->entities<2>()) {
    CINCH_CAPTURE() << "------- cell id: " << cell.id() << endl;
    for(auto corner : mesh->entities<0, 0, 1>(cell)) {
      CINCH_CAPTURE() << "--- corner id: " << corner.id() << endl;
    }
  }

  for(auto vertex : mesh->entities<0>()) {
    CINCH_CAPTURE() << "------- vertex id: " << vertex.id() << endl;
    for(auto corner : mesh->entities<0, 0, 1>(vertex)) {
      CINCH_CAPTURE() << "--- corner id: " << corner.id() << endl;
    }
  }

*/
  CINCH_WRITE("structured.blessed");

} // TEST
