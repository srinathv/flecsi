/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#ifndef flecsi_default_driver_h
#define flecsi_default_driver_h

#include <iostream>

#include "flecsi/utils/common.h"
#include "flecsi/execution/context.h"
#include "flecsi/execution/execution.h"
#include "flecsi/data/data.h"

/*!
 * \file default_driver.h
 * \authors bergen
 * \date Initial file creation: Jul 24, 2016
 */

namespace flecsi {
namespace execution {

/*----------------------------------------------------------------------------*
 * Fake mesh definition.
 *----------------------------------------------------------------------------*/

enum mesh_index_spaces_t : size_t {
  vertices,
  edges,
  faces,
  cells
}; // enum mesh_index_spaces_t

enum class privileges : size_t {
  read,
  read_write,
  write_discard
}; // enum data_access_type_t

struct mesh_t : public data::data_client_t {

  size_t indices(size_t index_space_id) override {

    switch(index_space_id) {
      case cells:
        return 10;
        break;
      default:
        // FIXME: lookup user-defined index space
        assert(false && "unknown index space");
    } // switch
  }

}; // struct mesh_t

// FIXME: Need to try to hide this
using namespace flecsi::data;
template<typename T>
using dense_field_t = storage_t::st_t<dense>::handle_t<double>;

/*----------------------------------------------------------------------------*
 * Task registration.
 *----------------------------------------------------------------------------*/

void task1(double dval, int ival) {
  std::cout << "Executing task1" << std::endl;
  std::cout << "Value(double): " << dval << std::endl;
  std::cout << "Value(int): " << ival << std::endl;
} // task1

register_task(task1, loc, void, double, int);

double task2(double x, double y, dense_field_t<double> p) {
  std::cout << "Executing task2" << std::endl;
  std::cout << "(x,y): (" << x << "," << y << ")" << std::endl;
  std::cout << "Return: " << x*y << std::endl;
//  return x*y;
} // task2

register_task(task2, loc, void, double, double, dense_field_t<double>);

/*----------------------------------------------------------------------------*
 * Kernel registration.
 *----------------------------------------------------------------------------*/

double eos_copper(double r, double e) {
  std::cout << "Executing function1" << std::endl;
  std::cout << "(r,e): (" << r << "," << e << ")" << std::endl;
  return r*e;
} // function1

register_function(eos_copper, double, double, double);

double eos_steel(double r, double e) {
  std::cout << "Executing function2" << std::endl;
  std::cout << "(r,e): (" << r << "," << e << ")" << std::endl;
  return 2*r*e;
} // function1

register_function(eos_steel, double, double, double);

/*----------------------------------------------------------------------------*
 * User type.
 *----------------------------------------------------------------------------*/

define_function_type(eos_function_t, double, double, double);

struct material_t {
  eos_function_t eos_function;
  double r;
  double e;

  double eos() {
    return execute_function(eos_function, r, e);
  } // eos

  material_t(eos_function_t eos_function_, double r_, double e_)
    : eos_function(eos_function_), r(r_), e(e_) {}
}; // struct material_t

struct copper_t : material_t {
  copper_t(double r_, double e_)
    : material_t(function_handle(eos_copper), r_, e_) {}
};

struct steel_t : material_t {
  steel_t(double r_, double e_)
    : material_t(function_handle(eos_steel), r_, e_) {}
};

/*----------------------------------------------------------------------------*
 * Driver.
 *----------------------------------------------------------------------------*/

void driver(int argc, char ** argv) {

  mesh_t m;

  // FIXME: need this to come from get_handle
  auto p = register_data(m, "pressure", 1, double, dense, cells);
  double alpha(10.0);

  execute_task(task1, loc, alpha, 5);
  execute_task(task2, loc,  alpha, 5.0, p);

  register_data(m, "materials", 1, material_t, dense, cells);

  auto mats = get_accessor(m, "materials", 0, material_t, dense);

  for(size_t i(0); i<5; ++i) {
    mats[i] = copper_t(2.0, 2.0);
  } // for

  for(size_t i(5); i<10; ++i) {
    mats[i] = steel_t(2.0, 2.0);
  } // for

  for(size_t i(0); i<10; ++i) {
    std::cout << mats[i].eos() << std::endl;
  } // for

} // driver

} // namespace execution
} // namespace flecsi

#endif // flecsi_default_driver_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/