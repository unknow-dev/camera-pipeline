#include<stdio.h>
#include<string>
#include "Halide.h"
#include "halide_image_io.h"
#include <math.h>
#include "clock.h"

#define R 0
#define G 1
#define B 2
using namespace std;
using namespace Halide::Tools;
using namespace Halide::ConciseCasts;
using namespace Halide;



int main(int argc, char **argv) {
  Buffer<uint8_t> in_b = Tools::load_image("../images/noise_test/in_n220.jpg");
  Func in = BoundaryConditions::repeat_edge(in_b);

  RDom r(-5, 6, -5, 6);

  Var x, y;
  
  Func box;
  box(x, y) = f32(0);
  box(x, y) = u8(sum(in(x + r.x, x + r.y) / 100.f));
  box.trace_stores();

  double t1 = current_time();
  Buffer<uint8_t> output = box.realize(in_b.width(), in_b.width(), in_.channels());
  double t2 = current_time();

  cout<<"Time: "<<t2-t1<<endl;
  save_image(output, "../images/noise_test/out_n220.jpg");
}
