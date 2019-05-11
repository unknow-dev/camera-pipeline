#include<stdio.h>
#include<string>
#include "Halide.h"
#include "halide_image_io.h"
#include <math.h>
#include "clock.h"
#define PI 3.141592f
using namespace std;
using namespace Halide::Tools;
using namespace Halide::ConciseCasts;
using namespace Halide;

int main(int argc, char **argv) {
  Halide::Buffer<uint8_t> in = load_image("../images/inputs/lena.jpg");
  float sigma = 2.0f;
  Expr K = 1.6f;

  Var x, y, c;
  RDom rk_1(-2, 3, -2, 3);
  Expr coef_1 = 1.0f/(2.0f * PI * sigma * sigma);
  
  Func diff_gauss_k_1;
  diff_gauss_k_1(x, y) = coef_1 * 1/(exp((x * x + y * y)/(2 * sigma * sigma)));

  RDom rk_2(-2, 3, -2, 3);
  Expr coef_2 = 1.0f/(2 * PI * K * K * sigma * sigma);
  Func diff_gauss_k_2;
  diff_gauss_k_2(x, y) = coef_2 * 1/(exp((x * x + y * y)/(2 * K * K * sigma * sigma)));
  
  Func in_b = BoundaryConditions::repeat_edge(in);
  Func conv_1;
  conv_1(x, y) += u8(in_b(x - rk_1.x, y - rk_1.y) * diff_gauss_k_1(rk_1.x, rk_1.y));

  Func conv_2;
  conv_2(x, y) += u8(in_b(x - rk_2.x, y - rk_2.y) * diff_gauss_k_2(rk_2.x, rk_2.y));

  
  Func img;
  img(x, y) = conv_1(x, y) - conv_2(x, y);

  double t1 = current_time(); 
  Buffer<uint8_t> output = conv_1.realize(in.width(), in.height());
  double t2 = current_time();
  cout<<"Time: "<< t2 - t1 <<endl;
  save_image(output, "../images/outputs/2_diff_gauss_green.jpg");   
} 
