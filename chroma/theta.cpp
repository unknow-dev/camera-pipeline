// TODO
// include sigma and lambda in the optimization
// optimize a_i
// fix supersampling

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
  Halide::Buffer<uint8_t> in = load_image("./images/inputs/2.jpg");

  Var x, y, c;
  Func in_bounded = BoundaryConditions::repeat_edge(in);

  Expr x_r = pow(in.width() - x, 2);
  Expr y_r = pow(in.height() - y, 2);
    

  float w_x = 1; //pixel width
  float w_y = 1; //pixel height
  Func cont_img_conv_h_box;
  RDom h_box(-2 / w_x, 3 / w_x, -2 / w_y, 3 / w_y);
  cont_img_conv_h_box(x, y, c) = sum(in_bounded(x + h_box.x, y, c) + in_bounded(x, y + h_box.y, c)) / 25.0f;
  
  Func img_d;
  img_d(x, y, c) = u8(cont_img_conv_h_box(3 * x, 3 * y,  c));

  double t1 = current_time(); 
  Buffer<uint8_t> output = img_d.realize(in.width()/3, in.height()/3, in.channels());
  double t2 = current_time();
  cout<<"Time: "<< t2 - t1 <<endl;
  save_image(output, "./images/outputs/conv_flower.jpg");
  
} 
