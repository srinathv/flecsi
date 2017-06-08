#include <iostream>
#include <cinchtest.h>
#include "flecsi/topology/structured_index_space.h"

using namespace std;
using namespace flecsi;
using namespace topology;

using id_vector_t = std::vector<size_t>;

TEST(structured, simple){

  using sid_t = structured_index_space<int>;
  id_vector_t indices;	
  size_t offset;

  //1D 
  CINCH_CAPTURE() << "1D Test "<< endl;
  sid_t cobj_1D;
  cobj_1D.init({0}, {2});
  ASSERT_EQ(cobj_1D.size(), 3);
  CINCH_CAPTURE() << " cobj_1D.size() = " << cobj_1D.size() << endl;

  for (auto obj : cobj_1D){
    CINCH_CAPTURE() << endl;
    CINCH_CAPTURE() << " obj --> " << obj << endl;
    indices = cobj_1D.get_indices_from_offset(obj);

    CINCH_CAPTURE() << "indices = [ ";
    for (auto i: indices)
     CINCH_CAPTURE() << i;
    CINCH_CAPTURE() << " ] " <<endl;

    offset = cobj_1D.get_offset_from_indices(indices);
    CINCH_CAPTURE() << "offset = " <<offset<<endl;
  }

  CINCH_CAPTURE() << endl;
  
  //2D
  CINCH_CAPTURE() << "2D Test "<< endl;
  sid_t cobj_2D;
  cobj_2D.init({0,0},{2,2});	
  ASSERT_EQ(cobj_2D.size(), 9);
  CINCH_CAPTURE() << " cobj.size() = " << cobj_2D.size() << endl;

  for (auto obj : cobj_2D){
    CINCH_CAPTURE() << endl;
    CINCH_CAPTURE() << " obj --> " << obj << endl;
    indices = cobj_2D.get_indices_from_offset(obj);

    CINCH_CAPTURE() << "indices = [ ";
    for (auto i: indices)
     CINCH_CAPTURE() << i<<" ";
    CINCH_CAPTURE() << "] " <<endl;

    offset = cobj_2D.get_offset_from_indices(indices);
    CINCH_CAPTURE() << "offset = " <<offset<<endl;

  }

  CINCH_CAPTURE() << endl;

  //3D
  CINCH_CAPTURE() << "3D Test "<< endl;
  sid_t cobj_3D;
  cobj_3D.init({0,0,0},{2,2,2});
  ASSERT_EQ(cobj_3D.size(), 27);
  CINCH_CAPTURE() << " cobj.size() = " << cobj_3D.size() << endl;

  for (auto obj : cobj_3D){
    CINCH_CAPTURE() << endl;
    CINCH_CAPTURE() << " obj --> " << obj << endl;
    indices = cobj_3D.get_indices_from_offset(obj);

    CINCH_CAPTURE() << "indices = [ ";
    for (auto i: indices)
     CINCH_CAPTURE() << i<<" ";
    CINCH_CAPTURE() << "] " <<endl;

    offset = cobj_3D.get_offset_from_indices(indices);
    CINCH_CAPTURE() << "offset = " <<offset<<endl;
  }

  CINCH_CAPTURE() << endl;
  
  //Output result
  CINCH_WRITE("structured_IS.blessed");
} // TEST
