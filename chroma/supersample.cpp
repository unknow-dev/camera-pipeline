// TODO
//fix the clamp bug
//optimize p_1D_f (v. slow), b (slow) and p_1D_b (v. v. slow)

#include<stdio.h>
#include<string>
#include "Halide.h"
#include "halide_image_io.h"
#include <math.h>
#include "clock.h"
using namespace std;
using namespace Halide::Tools;
using namespace Halide::ConciseCasts;
using namespace Halide;
int main(int argc, char **argv) {
  Halide::Buffer<uint8_t> in = load_image("../images/supersample_test/in256.jpg");

  Var x, y, c;
  
  Func in_b = BoundaryConditions::repeat_edge(in);

  Var j_x, j_y;
  
  Expr loc_x = x + j_x;
  Expr loc_y = y + j_y;
  
  Func f;
  f(x, y, c, j_x, j_y) = in_b(loc_x, loc_y, c);
 
  Expr t_x = x % 3;
  Expr t_x_w = t_x * 0.3333333f;

  Expr t_n_1_x = (-t_x_w + 2 * pow(t_x_w, 2));
  Expr t_0_x = (2 - 5 * pow(t_x_w, 2) + 3 * pow(t_x_w, 3));
  Expr t_1_x = (t_x_w + 4 * pow(t_x_w, 2) - 3 * pow(t_x_w, 3));
  Expr t_2_x = (-pow(t_x_w, 2) + pow(t_x_w, 3));
  
  Var i;
  Func b;
  b(x, y, c, i) = (f(x, y/3, c, -1, i) * t_n_1_x + f(x, y/3, c, 0, i) * t_0_x + f(x, y/3, c, 1, i) * t_1_x + f(x, y/3, c, 2, i) * t_2_x) / 2.0f;
 

  Expr t_y = y % 3;
  Expr t_y_w = t_y * 0.333333f;

  Expr t_n_1_y = -t_y_w + 2 * pow(t_y_w, 2);
  Expr t_0_y = 2 - 5 * pow(t_y_w, 2) + 3 * pow(t_y_w, 3);
  Expr t_1_y = t_y_w + 4 * pow(t_y_w, 2) - 3 * pow(t_y_w, 3);
  Expr t_2_y = -pow(t_y_w, 2) + pow(t_y_w, 3);
  
  Func p;
  p(x, y, c) = u8(min((b(x/3, y, c, -1) * t_n_1_y + b(x/3, y, c, 0) * t_0_y + b(x/3, y, c, 1) * t_1_y + b(x/3, y, c, 2) * t_2_y) / 2.0f, 255));
  double t1 = current_time(); 
  Buffer<uint8_t> output = p.realize(3 * in.width(), 3 * in.height(), in.channels());
  
  double t2 = current_time();
  cout<<"Time: "<< t2 - t1 <<endl;
  save_image(output, "../images/supersample_test/out256.jpg");
} 
