// g++ main.cpp u_new.a v_new.a -std=c++11 -I ../../include -I../../tools -L ../../bin `libpng-config --cflags --ldflags` -ljpeg -lpthread -ldl -o main
#include "HalideBuffer.h"
#include "halide_image_io.h"
#include <stdio.h>
#include "u_new.h"
#include "v_new.h"
#include "clock.h"
#include <iostream>
using namespace Halide;
using namespace Halide::Tools;

// bool converge() {
//   double diff_x, diff_y;
//   for (int x = 0; x < width; x++) {
//     for (int y = 0; x < height; y++) {
//       diff_x = v_new_buff(x, y, t) - v_buff(x, y, t);
//       diff_y = u_new_buff(x, y, t) - u_buff(x, y, t);
//     }
//   }
//   if (diff_x < threshold && diff_y < threshold) return true;
//   return false;
// }

std::string trim(std::string s) {
  std::string prepend = "";
  for (int i = 0; i < 4 - s.size(); i++) prepend += "0";

  return prepend + s;
  
}

void print_buffer(Halide::Runtime::Buffer<double> b) {
  for (int i = 0; i < b.height(); i++) {
    for (int j = 0; j < b.width(); j++) {
      printf("%f ", b(i, j));
    }
    printf("\n");
  }
}
int main(int argc, char **argv) {

  int n = 5; // we converge in chunks of n 

  int start_time = current_time();
  Halide::Runtime::Buffer<uint8_t> img[n];
  img[0] = load_image("test1/small/0001.png"); // to determine image width and height
  int W = img[0].width();
  int H = img[0].height();
  
  Halide::Runtime::Buffer<double> u_prev_f(W, H);
  Halide::Runtime::Buffer<double> prev_u(W, H);
  Halide::Runtime::Buffer<double> prev_v(W, H);
  
  int convergence = 10;
  
  Halide::Runtime::Buffer<double> u[n + 2];
    
  Halide::Runtime::Buffer<double> v[n + 2];

  for (int i = 0; i < n + 2; i++) {
    Halide::Runtime::Buffer<double> b1(W, H);
    Halide::Runtime::Buffer<double> b2(W, H);
    u[i] = b1;
    v[i] = b2;
  }

  int num_frames = 5;
  int start_frame = 40;
  for (int f = start_frame; f < start_frame + num_frames; f += n) {
    for (int i = 1; i < n; i++) {
      // std::cout<<trim(std::to_string(f + i - 1))<<std::endl;
      img[i - 1] = load_image("test1/small/" + trim(std::to_string(f + i - 1)) + ".png");
      img[i] = load_image("test1/small/" + trim(std::to_string(f + i)) + ".png");
    }
    
    for (int c = 0; c < convergence; c++) { // for now convergence happens at a constant rate

      u[0] = prev_u;
      v[0] = prev_v;
      
      for (int i = 1; i < n; i++) {
  	int err_u = u_new(img[i - 1].sliced(2, 0), img[i].sliced(2, 0), u[i - 1], u[i + 1], u[i], v[i], u[i]); // the result is stored in the last parameter
	// printf("frame %d u: \n", f + i);
	// print_buffer(u[i]);
	// printf("frame %d v: \n", f + i);
	// print_buffer(v[i]);
  	if (err_u) printf("Error at u_new pipline\n");
  	int err_v = v_new(img[i - 1].sliced(2, 0), img[i].sliced(2,0), u[i - 1], v[i - 1], u[i + 1], v[i + 1], u[i], v[i], v[i]);
  	if (err_v) printf("Error at v_new pipline\n");
      }
    }
    prev_u = u[n];
    prev_v = v[n];
  }

  printf("optical flow 15 frames time: %f ms", (current_time() - start_time) / 100);
}
