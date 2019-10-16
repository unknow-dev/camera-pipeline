#include "Halide.h"
#include "halide_image_io.h"
#include <stdio.h>
#include "clock.h"
using namespace Halide;
using namespace Halide::Tools;
using namespace Halide::ConciseCasts;
bool have_opencl_or_metal();
int main(int argc, char **argv) {
  
  double start_time = current_time();
  ImageParam img0_orig(type_of<uint8_t>(), 2);
  ImageParam img1_orig(type_of<uint8_t>(), 2);
  
  Expr W = img0_orig.width();
  Expr H = img0_orig.height();
  Var x, y;
  Expr clamped_x = clamp(x, 3, img1_orig.width() - 5);
  Expr clamped_y = clamp(y, 3, img1_orig.height() - 5);
  
  Func img0;
  img0(x, y) = img0_orig(clamped_x, clamped_y);
  Func img1;
  img1(x, y) = img1_orig(clamped_x, clamped_y);

  ImageParam prev_u(type_of<double>(), 2);
  ImageParam prev_v(type_of<double>(), 2);
  
  ImageParam next_u(type_of<double>(), 2);
  ImageParam next_v(type_of<double>(), 2);

  ImageParam u(type_of<double>(), 2);
  ImageParam v(type_of<double>(), 2);

  
  Func fx;

  fx(x, y) = (i16(img1(x - 3, y)) + i16(img1(x - 2, y)) * -9 + i16(img1(x - 1, y)) * 45 + i16(img1(x + 1, y)) * -45 + i16(img1(x + 2, y)) * 9 + i16(img1(x + 3, y)) * -1) / (60 * W);

  Func fy;
  fy(x, y) = (i16(img1(x, y - 3)) + i16(img1(x, y - 2)) * -9 + i16(img1(x, y - 1)) * 45 + i16(img1(x, y + 1)) * -45 + i16(img1(x, y + 2)) * 9 + i16(img1(x, y + 3)) * -1) / (60 * H);

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

  Var xo, yo, xi, yi;
  Var block, thread;
  Expr omega = 1.9f;
  
  Func u_new;
  u_new(x, y) = f64(0);
  Func v_new;
  v_new(x, y) = f64(0);
  
  // role of alpha from pg 215
  // the smoothness weight α > 0 serves as regularisation parameter: Larger values for α result in a stronger penalisation of large flow gradients and lead to smoother flow fields.
  Expr alpha = 100000.0f;
  
  RDom img_r(3, W - 4, 3, H - 4);
  // update
  Expr u_neg_sum = u_new(img_r.x - 1, img_r.y) + u_new(img_r.x, img_r.y - 1) + prev_u(img_r.x, img_r.y); // prev frame and prev pixels
  Expr u_pos_sum = u(img_r.x + 1, img_r.y) + u(img_r.x, img_r.y + 1) + next_u(img_r.x, img_r.y); // next frame and next pixels
 
  u_new(img_r.x, img_r.y) = (1 - omega) * u(img_r.x, img_r.y) + omega * (u_neg_sum + u_pos_sum - (W * H / alpha) * (J[0][1](img_r.x, img_r.y) * v(img_r.x, img_r.y) + J[0][2](img_r.x, img_r.y))) / (6 + H * W / alpha * J[0][0](img_r.x, img_r.y)); // 6 because there are 6 neighboring values that we are computing the new value
  //u_new.trace_stores();
  
  Expr v_neg_sum = v_new(img_r.x - 1, img_r.y) + v_new(img_r.x, img_r.y - 1) + prev_v(img_r.x, img_r.y); // prev frame and prev pixels
  Expr v_pos_sum = v(img_r.x + 1, img_r.y) + v(img_r.x, img_r.y + 1) + next_v(img_r.x, img_r.y); // next frame and next pixels
  // update
  v_new(img_r.x, img_r.y) = (1 - omega) * v(img_r.x, img_r.y) + omega * (v_neg_sum + v_pos_sum - (W * H / alpha) * (J[1][0](img_r.x, img_r.y) * u_new(img_r.x, img_r.y) + J[1][2](img_r.x, img_r.y))) / (6 + (W * H) / 6 * J[1][1](img_r.x, img_r.y));

  //v_new.trace_stores();
  
  Target target = get_host_target();
  
  if (target.os == Target::OSX) {
    target.set_feature(Target::Metal);
  } else {
    target.set_feature(Target::OpenCL);
  }
  
  
  u_new.compile_to_static_library("u_new", {img0_orig, img1_orig, prev_u, next_u, u, v}, "u_new");
  v_new.compile_to_static_library("v_new", {img0_orig, img1_orig, prev_u, prev_v, next_u, next_v, u, v}, "v_new");

  if(!have_opencl_or_metal()) {
    printf("opencl or metal not found\n");
  } else {
    printf("have opencl or metal\n");
    fx.gpu_tile(x, y, xo, yo, xi, yi, 8, 8);
    fx.compute_root();
    fy.gpu_tile(x, y, xo, yo, xi, yi, 8, 8);
    fy.compute_root();
    ft.gpu_tile(x, y, xo, yo, xi, yi, 8, 8);
    ft.compute_root();
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
	J[i][j].compute_root();
	J[i][j].gpu_tile(x, y, xo, yo, xi, yi, 8, 8);
	J[i][j].trace_stores();
      }
    }

    u_new.gpu_tile(x, y, xo, yo, xi, yi, 8, 8);
    v_new.gpu_tile(x, y, xo, yo, xi, yi, 8, 8);
    
  }
  

  std::cout<<(current_time() - start_time) / 100 <<" ms"<<std::endl;
}


// A helper function to check if OpenCL seems to exist on this machine.

#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

bool have_opencl_or_metal() {
#ifdef _WIN32
    return LoadLibrary("OpenCL.dll") != NULL;
#elif __APPLE__
    return dlopen("/System/Library/Frameworks/Metal.framework/Versions/Current/Metal", RTLD_LAZY) != NULL;
#else
    return dlopen("libOpenCL.so", RTLD_LAZY) != NULL;
#endif
}

