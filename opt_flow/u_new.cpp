// g++ pipeline.cpp -g -I ../../include -I ../../tools -L ../../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -o pipeline -std=c++11
// DYLD_LIBRARY_PATH=../../bin ./pipeline
#include "Halide.h"
#include "halide_image_io.h"
#include <stdio.h>
#include "clock.h"
using namespace Halide;
using namespace Halide::Tools;
using namespace Halide::ConciseCasts;
int height;
int width;
int t0;

int main(int argc, char **argv) {

  double start_time = current_time();
  ImageParam img0(type_of<int16_t>(), 2);
  ImageParam img1(type_of<int16_t>(), 2);

  ImageParam prev_u(type_of<double>(), 2);
  ImageParam prev_v(type_of<double>(), 2);
  
  ImageParam next_u(type_of<double>(), 2);
  ImageParam next_v(type_of<double>(), 2);

  ImageParam u(type_of<double>(), 2);
  ImageParam v(type_of<double>(), 2);

  Var x, y, t;
  Func fx;

  fx(x, y) = (img1(x - 3, y) + img1(x - 2, y) * -9 + img1(x - 1, y) * 45 + img1(x + 1, y) * -45 + img1(x + 2, y) * 9 + img1(x + 3, y) * -1) / (60 * img1.width());

  Func fy;
  fy(x, y) = (img1(x, y - 3) + img1(x, y - 2) * -9 + img1(x, y - 1) * 45 + img1(x, y + 1) * -45 + img1(x, y + 2) * 9 + img1(x, y + 3) * -1) / (60 * img1.height());

  Func ft;
  ft(x, y) = img1(x, y) - img0(x, y);

  Func J[3][3];
  J[0][0](x, y) = pow(fx(x, y), 2);
  J[0][1](x, y) = fx(x, y) * fy(x, y);
  J[0][2](x, y) = fx(x, y) * ft(x, y);
  J[1][0](x, y) = fx(x, y) * fy(x, y);
  J[1][1](x, y) = pow(fy(x, y), 2);
  J[1][2](x, y) = fy(x, y) * ft(x, y);
  J[2][0](x, y) = fx(x, y) * ft(x, y);
  J[2][1](x, y) = fx(x, y) * ft(x, y);
  J[2][2](x, y) = pow(ft(x, y), 2);

  Expr omega = 1.9f;
  
  Func u_new;
  u_new(x, y) = f64(0);
  Func v_new;
  v_new(x, y) = f64(0);
  
  // role of alpha from pg 215
  // the smoothness weight α > 0 serves as regularisation parameter: Larger values for α result in a stronger penalisation of large flow gradients and lead to smoother flow fields.
  Expr alpha = 100000.0f;
  
  Expr W = img0.width();
  Expr H = img0.height();
  RDom img_r(0, W, 0, H);
  // update
  Expr u_neg_sum = u_new(img_r.x - 1, img_r.y) + u_new(img_r.x, img_r.y - 1) + prev_u(img_r.x, img_r.y); // prev frame and prev pixels
  Expr u_pos_sum = u(img_r.x + 1, img_r.y) + u(img_r.x, img_r.y + 1) + next_u(img_r.x, img_r.y); // next frame and next pixels

 
 
  u_new(img_r.x, img_r.y) = (1 - omega) * u(img_r.x, img_r.y) + omega * (u_neg_sum + u_pos_sum - (W * H / alpha) * (J[0][1](img_r.x, img_r.y) * v(img_r.x, img_r.y) + J[0][2](img_r.x, img_r.y))) / (6 + H * W / alpha * J[0][0](img_r.x, img_r.y)); // 6 because there are 6 neighboring values that we are computing the new value
  u_new.trace_stores();
  
  Expr v_neg_sum = v_new(img_r.x - 1, img_r.y) + v_new(img_r.x, img_r.y - 1) + prev_v(img_r.x, img_r.y); // prev frame and prev pixels
  Expr v_pos_sum = v(img_r.x + 1, img_r.y) + v(img_r.x, img_r.y + 1) + next_v(img_r.x, img_r.y); // next frame and next pixels
  // update
  v_new(img_r.x, img_r.y) = (1 - omega) * v(img_r.x, img_r.y) + omega * (v_neg_sum + v_pos_sum - (W * H / alpha) * (J[1][0](img_r.x, img_r.y) * u_new(img_r.x, img_r.y) + J[1][2](img_r.x, img_r.y))) / (6 + (W * H) / 6 * J[1][1](img_r.x, img_r.y));

  v_new.trace_stores();
  
  u_new.compile_to_static_library("u_new", {img0, img1, prev_u, next_u, u, v}, "u_pipe");

  std::cout<<(current_time() - start_time) / 100 <<" ms"<<std::endl;
}


