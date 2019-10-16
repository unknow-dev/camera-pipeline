// g++ test -g -I ../../include -L ../../bin -lHalide -o test -std=c++11
// DYLD_LIBRARY_PATH=../bin ./test
#include <stdio.h>
#include "Halide.h"
using namespace Halide;
int main(int argc, char **argv) {
  // test 1d example
  std::cout<<"Are we printing anything at all?"<<std::endl;
  Var x;
  Func v;
  v(x) = 0;
  v(0) = 1;
  v(x) = 10;
  v.trace_stores();
  Buffer<uint8_t> b = v.realize(20);


  // Func gradient("gradient");
  // gradient(x, y) = x + y;

  // // And tell Halide that we'd like to be notified of all
  // // evaluations.
  // gradient.trace_stores();

  // // Realize the function over an 8x8 region.
  // printf("Evaluating gradient\n");
  // Buffer<int> output = gradient.realize(8, 8);
}

