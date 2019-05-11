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
const int n_channels = 3;

int main(int argc, char **argv) {
  Halide::Buffer<uint8_t> in = load_image("./images/inputs/flower_noisy.jpg");
 
  Var x, y, c, i;
  Func a0;
  a0(i) = i;
  Func a1;
  a1(i) = i;
  Func a2;
  a2(i) = i;
  Func a3;
  a3(i) = i;
  
  Func r;
  r(x, y, c) = sqrt((in.width() / 2 - x) * (in.width() / 2 - x) + (in.height() / 2 - y) * (in.height() / 2 - y));
  
  Func mag_inv;
  mag_inv(x, y, c) = in(x, y, c) / (a0(0) + a1(1) * r(x, y, c) + a2(2) * r(x, y, c) + a3(3) * r(x, y, c));
  Func h1;
  h1(x, y, c) = 0;
  h1(1, 1, c) = 1;

  //rect_2(x, y, c) =
  //float w_x, w_y;
   //h_box(x, y, w_x, w_y, c) = rect_2(x / w_x, y / w_y, c);


  float sigma = 1.5f;
  Func kernel;
    kernel(x) = exp(-x*x/(2*sigma*sigma)) / (sqrtf(2*M_PI)*sigma);

    Func in_bounded = BoundaryConditions::repeat_edge(in);

    Func blur_y;
    blur_y(x, y, c) = (kernel(0) * in_bounded(x, y, c) +
		       kernel(1) * (in_bounded(x, y-1, c) +
				    in_bounded(x, y+1, c)) +
		       kernel(2) * (in_bounded(x, y-2, c) +
				    in_bounded(x, y+2, c)) +
		       kernel(3) * (in_bounded(x, y-3, c) +
				    in_bounded(x, y+3, c)));
    Func h_sigma;
    h_sigma(x, y, c) = (kernel(0) * blur_y(x, y, c) +
		       kernel(1) * (blur_y(x-1, y, c) +
				    blur_y(x+1, y, c)) +
		       kernel(2) * (blur_y(x-2, y, c) +
				    blur_y(x+2, y, c)) +
		       kernel(3) * (blur_y(x-3, y, c) +
				    blur_y(x+3, y, c)));

    float lambda = 0.0005;
    RDom k(0, 2, 0, 2);
    Expr first_order_term = lambda * (h1(x + k.x, y + k.y, c) - h_sigma(x + k.x, y + k.y, c));
    Expr second_order_term = lambda * lambda * h1(x + k.x, y + k.y, c) - h_sigma(x + k.x, y + k.y, c) * h1(x + k.x, y + k.y, c) - h_sigma(x + k.x, y + k.y, c);
    
    Func h_sharp_inv;
    h_sharp_inv(x, y, c) = f32(0);
    h_sharp_inv(x, y, c) = h1(k.x, k.y, c) - first_order_term + second_order_term;
    
    Func h_sharp_inv_u8;
    h_sharp_inv_u8(x, y, c) = u8(h_sharp_inv(x, y, c));

    Func h_sharp;
    h_sharp(x, y, c) = f32(0);
    h_sharp(x, y, c) = h1(k.x, k.y, c) + first_order_term;

    Func h_sharp_u8;
    h_sharp_u8(x, y, c) = u8(h_sharp(x, y, c));
    
    double t1 = current_time(); 
    Buffer<uint8_t> output = h_sharp_inv_u8.realize(in.width(), in.height(), in.channels());
    double t2 = current_time();
    cout<<"Time: "<< t2 - t1 <<endl;
    
    save_image(output, "./images/outputs/h_inv_flower00005.jpg");
} 
