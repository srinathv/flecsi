#include <iostream>
#include <cinchtest.h>
#include "flecsi/topology/structured_mesh_topology.h"
#include <sys/time.h>


using namespace std;
using namespace flecsi;
using namespace topology;


double wtime()
{
  double y = -1;
  struct timeval tm;
  gettimeofday(&tm,NULL);
  y = (double)(tm.tv_sec)+(double)(tm.tv_usec)*1.e-6;
  return y;
};

struct smi_perf
{
  double meshgen;
  double v2v, v2e, v2f, v2c;
  double e2v, e2e, e2f, e2c;
  double f2v, f2e, f2f, f2c;
  double c2v, c2e, c2f, c2c;
};

class Vertex : public structured_mesh_entity_t<0, 1>{
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
  static constexpr std::array<size_t,num_dimensions> upper_bounds = {160,160,160};

  using entity_types = std::tuple<
  std::pair<domain_<0>, Vertex>,
  std::pair<domain_<0>, Cell>>;

};

constexpr std::array<size_t,TestMesh3dType::num_dimensions> TestMesh3dType::lower_bounds;
constexpr std::array<size_t,TestMesh3dType::num_dimensions> TestMesh3dType::upper_bounds;

using id_vector_t = std::vector<size_t>;
using TestMesh = structured_mesh_topology_t<TestMesh3dType>;

TEST(structured3dtime, simple){
  
  smi_perf *sp = new smi_perf;
  double time_start, time_end;
  

  time_start = wtime();
  auto mesh = new TestMesh;
  sp->meshgen = wtime() - time_start;
 
  size_t nv, nex, ney, nez, nfx, nfy, nfz, nc, val;
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
  
  CINCH_CAPTURE()<<"|V| = "<<nv<<", |E| = "<<nex+ney+nez<<", |F| = "<<nfx+nfy+nfz<<", |C| = "<<nc<<endl;
  
  //Loop over all vertices and test intra index space queries
  //V-->V
  time_start = wtime();
  for (auto vertex: mesh->entities<0>()){
    val = mesh->entities< 0, 0, 1, 0, 0 >(vertex);
    val = mesh->entities< 0, 0,-1, 0, 0 >(vertex); 
    val = mesh->entities< 0, 0, 0, 1, 0 >(vertex);
    val = mesh->entities< 0, 0, 0,-1, 0 >(vertex);
    val = mesh->entities< 0, 0, 0, 0, 1 >(vertex);
    val = mesh->entities< 0, 0, 0, 0,-1 >(vertex);
  }
  sp->v2v = (wtime()-time_start)/(double)(6*nv);
 
  //V-->E
  time_start = wtime();
  for (auto vertex: mesh->entities<0>())
  {
    adj.clear();
    adj = mesh->get_entities<0,1,0>(vertex);
  }
  sp->v2e = (wtime()-time_start)/(double)(nv);
   
  //V-->F
  time_start = wtime();
  for (auto vertex: mesh->entities<0>())
  {
    adj.clear(); 
    adj = mesh->get_entities<0,2,0>(vertex);
  }
  sp->v2f = (wtime()-time_start)/(double)(nv);

  //V-->C
  time_start = wtime();
  for (auto vertex: mesh->entities<0>())
  {
    adj.clear();
    adj = mesh->get_entities<0,3,0>(vertex);
  }
  sp->v2c = (wtime()-time_start)/(double)(nv);
  
  
  //Loop over all edges in X-direction and test intra index space queries
  //E-->V
  time_start = wtime();
  for (auto edgex: mesh->entities<1>())
   {
     adj.clear();
     adj = mesh->get_entities<1,0,0>(edgex);
   }

   for (auto edgey: mesh->entities<2>())
   {
     adj.clear();
     adj = mesh->get_entities<1,0,0>(edgey+nex);
   }

   for (auto edgez: mesh->entities<3>())
   {
     adj.clear();
     adj = mesh->get_entities<1,0,0>(edgez+nex+ney);
   }
  sp->e2v = (wtime()-time_start)/(double)(nex+ney+nez); 
  
  //E-->E
  time_start = wtime();
  for (auto edgex: mesh->entities<1>())
  {
    val = mesh->entities< 1, 0, 1, 0, 0 >(edgex);
    val = mesh->entities< 1, 0,-1, 0, 0 >(edgex); 
    val = mesh->entities< 1, 0, 0, 1, 0 >(edgex);
    val = mesh->entities< 1, 0, 0,-1, 0 >(edgex);
    val = mesh->entities< 1, 0, 0, 0, 1 >(edgex);
    val = mesh->entities< 1, 0, 0, 0,-1 >(edgex);
  }

  for (auto edgey: mesh->entities<2>())
  {
    val = mesh->entities< 1, 0, 1, 0, 0 >(edgey+nex);
    val = mesh->entities< 1, 0,-1, 0, 0 >(edgey+nex); 
    val = mesh->entities< 1, 0, 0, 1, 0 >(edgey+nex);
    val = mesh->entities< 1, 0, 0,-1, 0 >(edgey+nex);
    val = mesh->entities< 1, 0, 0, 0, 1 >(edgey+nex);
    val = mesh->entities< 1, 0, 0, 0,-1 >(edgey+nex);
  }

  for (auto edgez: mesh->entities<3>())
  {
    val = mesh->entities< 1, 0, 1, 0, 0 >(edgez+nex+ney);
    val = mesh->entities< 1, 0,-1, 0, 0 >(edgez+nex+ney); 
    val = mesh->entities< 1, 0, 0, 1, 0 >(edgez+nex+ney);
    val = mesh->entities< 1, 0, 0,-1, 0 >(edgez+nex+ney);
    val = mesh->entities< 1, 0, 0, 0, 1 >(edgez+nex+ney);
    val = mesh->entities< 1, 0, 0, 0,-1 >(edgez+nex+ney);
  }
  sp->e2e = (wtime()-time_start)/(double)(6*(nex+ney+nez));

  //E-->F
  time_start = wtime();
  for (auto edgex: mesh->entities<1>())
   {
     adj.clear();
     adj = mesh->get_entities<1,2,0>(edgex);
   }

  for (auto edgey: mesh->entities<2>())
   {
     adj.clear();
     adj = mesh->get_entities<1,2,0>(edgey+nex);
   }

  for (auto edgez: mesh->entities<3>())
   {
     adj.clear();
     adj = mesh->get_entities<1,2,0>(edgez+nex+ney);
   }
   sp->e2f = (wtime()-time_start)/(double)(nex+ney+nez);

  //E-->C
   time_start = wtime();
  for (auto edgex: mesh->entities<1>())
   {
     adj.clear();
     adj = mesh->get_entities<1,3,0>(edgex);
   }

  for (auto edgey: mesh->entities<2>())
   {
     adj.clear();
     adj = mesh->get_entities<1,3,0>(edgey+nex);
   }

  for (auto edgez: mesh->entities<3>())
   {
     adj.clear();
     adj = mesh->get_entities<1,3,0>(edgez+nex+ney);
   }
   sp->e2c = (wtime()-time_start)/(double)(nex+ney+nez);
  
  //F-->V
   time_start = wtime();
  for (auto facex: mesh->entities<4>())
   {
     adj.clear();
     adj = mesh->get_entities<2,0,0>(facex);
   }

   for (auto facey: mesh->entities<5>())
   {
     adj.clear();
     adj = mesh->get_entities<2,0,0>(facey+nfx);
   }

   for (auto facez: mesh->entities<6>())
   {
     adj.clear();
     adj = mesh->get_entities<2,0,0>(facez+nfx+nfy);
   }
   sp->f2v = (wtime()-time_start)/(double)(nfx+nfy+nfz);

  //F-->E
   time_start = wtime();
    for (auto facex: mesh->entities<4>())
   {
     adj.clear();
     adj = mesh->get_entities<2,1,0>(facex);
   }

   for (auto facey: mesh->entities<5>())
   {
     adj.clear();
     adj = mesh->get_entities<2,1,0>(facey+nfx);
   }

   for (auto facez: mesh->entities<6>())
   {
     adj.clear();
     adj = mesh->get_entities<2,1,0>(facez+nfx+nfy);
   }
   sp->f2e = (wtime()-time_start)/(double)(nfx+nfy+nfz);

   //F-->F
   time_start = wtime();
   for (auto facex: mesh->entities<4>())
  {
    val = mesh->entities< 2, 0, 1, 0, 0 >(facex);
    val = mesh->entities< 2, 0,-1, 0, 0 >(facex); 
    val = mesh->entities< 2, 0, 0, 1, 0 >(facex);
    val = mesh->entities< 2, 0, 0,-1, 0 >(facex);
    val = mesh->entities< 2, 0, 0, 0, 1 >(facex);
    val = mesh->entities< 2, 0, 0, 0,-1 >(facex);
  }

  for (auto facey: mesh->entities<5>())
  {
    val = mesh->entities< 2, 0, 1, 0, 0 >(facey+nfx);
    val = mesh->entities< 2, 0,-1, 0, 0 >(facey+nfx); 
    val = mesh->entities< 2, 0, 0, 1, 0 >(facey+nfx);
    val = mesh->entities< 2, 0, 0,-1, 0 >(facey+nfx);
    val = mesh->entities< 2, 0, 0, 0, 1 >(facey+nfx);
    val = mesh->entities< 2, 0, 0, 0,-1 >(facey+nfx);
  }

  for (auto facez: mesh->entities<6>())
  {
    val = mesh->entities< 2, 0, 1, 0, 0 >(facez+nfx+nfy);
    val = mesh->entities< 2, 0,-1, 0, 0 >(facez+nfx+nfy); 
    val = mesh->entities< 2, 0, 0, 1, 0 >(facez+nfx+nfy);
    val = mesh->entities< 2, 0, 0,-1, 0 >(facez+nfx+nfy);
    val = mesh->entities< 2, 0, 0, 0, 1 >(facez+nfx+nfy);
    val = mesh->entities< 2, 0, 0, 0,-1 >(facez+nfx+nfy);
  }
  sp->f2f = (wtime()-time_start)/(double)(6*(nfx+nfy+nfz));

  //F-->C
  time_start = wtime();
    for (auto facex: mesh->entities<4>())
   {
     adj.clear();
     adj = mesh->get_entities<2,3,0>(facex);
   }

   for (auto facey: mesh->entities<5>())
   {
     adj.clear();
     adj = mesh->get_entities<2,3,0>(facey+nfx);
   }

   for (auto facez: mesh->entities<6>())
   {
     adj.clear();
     adj = mesh->get_entities<2,3,0>(facez+nfx+nfy);
   }
   sp->f2c = (wtime()-time_start)/(double)(nfx+nfy+nfz);

  //C-->V
   time_start = wtime();
   for (auto cell: mesh->entities<7>())
   {
     adj.clear();
     adj = mesh->get_entities<3,0,0>(cell);
   }
   sp->c2v = (wtime()-time_start)/(double)(nc);

  //C-->E
   time_start = wtime();
  for (auto cell: mesh->entities<7>())
   {
     adj.clear();
     adj = mesh->get_entities<3,1,0>(cell);
   }
   sp->c2e = (wtime()-time_start)/(double)(nc);

  //C-->F
   time_start = wtime();
  for (auto cell: mesh->entities<7>())
   {
     adj.clear();
     adj = mesh->get_entities<3,2,0>(cell);
   }
   sp->c2f = (wtime()-time_start)/(double)(nc);
 
  //C-->C
   time_start = wtime();
  for (auto cell: mesh->entities<7>())
  {
    val = mesh->entities< 3, 0, 1, 0, 0 >(cell);
    val = mesh->entities< 3, 0,-1, 0, 0 >(cell); 
    val = mesh->entities< 3, 0, 0, 1, 0 >(cell);
    val = mesh->entities< 3, 0, 0,-1, 0 >(cell);
    val = mesh->entities< 3, 0, 0, 0, 1 >(cell);
    val = mesh->entities< 3, 0, 0, 0,-1 >(cell);
  }
  sp->c2c = (wtime()-time_start)/(double)(6*nc);

  CINCH_CAPTURE()<<"Average times in seconds for mesh traversal\n"<<endl;
  CINCH_CAPTURE()<<"Generating structured mesh :: "<<sp->meshgen<<endl;
  CINCH_CAPTURE()<<endl;
  CINCH_CAPTURE()<<"QUERY: Vertex --> Vertex :: "<<sp->v2v<<endl;
  CINCH_CAPTURE()<<"QUERY: Vertex --> Edge :: "<<sp->v2e<<endl;
  CINCH_CAPTURE()<<"QUERY: Vertex --> Face :: "<<sp->v2f<<endl;
  CINCH_CAPTURE()<<"QUERY: Vertex --> Cell :: "<<sp->v2c<<endl;
  CINCH_CAPTURE()<<endl;
  CINCH_CAPTURE()<<"QUERY: Edge --> Vertex :: "<<sp->e2v<<endl;
  CINCH_CAPTURE()<<"QUERY: Edge --> Edge :: "<<sp->e2e<<endl;
  CINCH_CAPTURE()<<"QUERY: Edge --> Face :: "<<sp->e2f<<endl;
  CINCH_CAPTURE()<<"QUERY: Edge --> Cell :: "<<sp->e2c<<endl;
  CINCH_CAPTURE()<<endl;
  CINCH_CAPTURE()<<"QUERY: Face --> Vertex :: "<<sp->f2v<<endl;
  CINCH_CAPTURE()<<"QUERY: Face --> Edge :: "<<sp->f2e<<endl;
  CINCH_CAPTURE()<<"QUERY: Face --> Face :: "<<sp->f2f<<endl;
  CINCH_CAPTURE()<<"QUERY: Face --> Cell :: "<<sp->f2c<<endl;
  CINCH_CAPTURE()<<endl;
  CINCH_CAPTURE()<<"QUERY: Cell --> Vertex :: "<<sp->c2v<<endl;
  CINCH_CAPTURE()<<"QUERY: Cell --> Edge :: "<<sp->c2e<<endl;
  CINCH_CAPTURE()<<"QUERY: Cell --> Face :: "<<sp->c2f<<endl;
  CINCH_CAPTURE()<<"QUERY: Cell --> Cell :: "<<sp->c2c<<endl;

  CINCH_WRITE("structured3dtime.blessed");

} // TEST
