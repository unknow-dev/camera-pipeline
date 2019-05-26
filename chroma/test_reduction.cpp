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

  Var x, y, c;
  Func box;
  box(x, y, c) = u8(0);
  box(x, y, c) += u8((in(x + r.x, y + r.y, c) / 121.f));

  Func box_;
  box_(x, y, c) = u8(0);
  box_(x, y, c) = u8(sum(in(x + r.x, y + r.y, c))/ 11.f);

  Func diff;
  diff(x, y, c) = box(x, y, c) - box_(x, y, c);
  diff.trace_stores();
  double t1 = current_time();
  Buffer<uint8_t> output = box.realize(in_b.width(), in_b.width(), in_b.channels());
  Buffer<uint8_t> output_ = box_.realize(in_b.width(), in_b.width(), in_b.channels());
  double t2 = current_time();
  cout<<"Time: "<<t2-t1<<endl;
  save_image(output, "../images/noise_test/box.jpg");
  save_image(output_, "../images/noise_test/box_.jpg");
}
