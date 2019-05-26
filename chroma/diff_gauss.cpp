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
  Halide::Buffer<uint8_t> in = load_image("../images/chroma_test/scenery.jpg");
  float sigma_1 = 0.25f;
  float sigma_2 = 0.75f;
  
  Var x, y, c;

  Func diff_gauss_1_k;
  diff_gauss_1_k(x) = exp((- x * x) / (2.0f * sigma_1))/(sigma_1 * sqrtf(2 * M_PI));
  Func diff_gauss_2_k;
  diff_gauss_2_k(x) = exp((- x * x) / (2.0f * sigma_1))/(sigma_2 * sqrtf(2 * M_PI));
  
  Func in_bounded = BoundaryConditions::repeat_edge(in);
  Func conv1_y;
  conv1_y(x, y, c) = diff_gauss_1_k(0) * in_bounded(x, y, c) +
    diff_gauss_1_k(1) * (in_bounded(x, y - 1, c) + in_bounded(x, y + 1, c)) +
    diff_gauss_1_k(1) * (in_bounded(x, y - 2, c) + in_bounded(x, y + 2, c)) +
    diff_gauss_1_k(1) * (in_bounded(x, y - 3, c) + in_bounded(x, y + 3, c));

  Func conv1;
  conv1(x, y, c) = diff_gauss_1_k(0) * conv1_y(x, y, c) +
    diff_gauss_1_k(1) * (conv1_y(x - 1, y, c) + conv1_y(x + 1, y, c)) +
    diff_gauss_1_k(1) * (conv1_y(x - 2, y, c) + conv1_y(x + 2, y, c)) +
    diff_gauss_1_k(1) * (conv1_y(x - 3, y, c) + conv1_y(x + 3, y, c));

  Func conv2_y;
  conv2_y(x, y, c) = diff_gauss_2_k(0) * in_bounded(x, y, c) +
    diff_gauss_2_k(1) * (in_bounded(x, y - 1, c) + in_bounded(x, y + 1, c)) +
    diff_gauss_2_k(1) * (in_bounded(x, y - 2, c) + in_bounded(x, y + 2, c)) +
    diff_gauss_2_k(1) * (in_bounded(x, y - 3, c) + in_bounded(x, y + 3, c));

  Func conv2;
  conv2(x, y, c) = diff_gauss_2_k(0) * conv2_y(x, y, c) +
    diff_gauss_2_k(1) * (conv1_y(x - 1, y, c) + conv2_y(x + 1, y, c)) +
    diff_gauss_2_k(1) * (conv1_y(x - 2, y, c) + conv2_y(x + 2, y, c)) +
    diff_gauss_2_k(1) * (conv1_y(x - 3, y, c) + conv2_y(x + 3, y, c));

  
  Func img;
  img(x, y, c) = u8(conv1(x, y, c) - conv2(x, y, c));

  img.trace_stores();
  Func r_e;
  r_e(x, y) = img(x, y, R);
  r_e.trace_stores();
  
  Func g_e;
  g_e(x, y) = img(x, y, G);
  Func b_e;
  b_e(x, y) = img(x, y, B);

  // double t1 = current_time(); 
  // Buffer<uint8_t> r_e_b = r_e.realize(in.width(), in.height());
  //Buffer<uint8_t> g_e_b = g_e.realize(in.width(), in.height());
  //Buffer<uint8_t> b_e_b = b_e.realize(in.width(), in.height());
  //double t2 = current_time();
  //cout<<"Time: "<< t2 - t1 <<endl;
  // save_image(r_e_b, "../images/chroma_test/red.jpg");
  // save_image(g_e_b, "../images/chroma_test/green.jpg");
  // save_image(b_e_b, "../images/chroma_test/blue.jpg");
  
  ImageParam sigma(Float(64), 1); //optimize this parameter
  ImageParam lambda(Float(64), 1); //optimize this parameter
  
  Buffer<double> sigma_b(1);
  Buffer<double> lambda_b(1);

  sigma.set(sigma_b);
  lambda.set(lambda_b);
  
  
  Expr gauss = exp(-x * x / (2.f * sigma(0) * sigma(0))) / (sqrtf(2.f * M_PI) * sigma(0));
  Expr gauss_0 = 1.0f / (sqrtf(2.f * M_PI) * sigma(0));
 
  Func h_sharp_inv_k;
  h_sharp_inv_k(x) = - (lambda(0) * -gauss) + (lambda(0) * lambda(0) * gauss * gauss);
  h_sharp_inv_k(0) = 1.0f - lambda(0) * (1.0f - gauss_0) + lambda(0) * lambda(0) * (1.0f - gauss_0) * (1.0f - gauss_0);

  
  Func r_e_y;
  r_e_y(x, y) = f32(h_sharp_inv_k(0) * r_e(x, y) +
		 h_sharp_inv_k(1) * (r_e(x, y - 1) + r_e(x, y + 1)) +
		 h_sharp_inv_k(2) * (r_e(x, y - 2) + r_e(x, y + 2)) +
		 h_sharp_inv_k(3) * (r_e(x, y - 3) + r_e(x, y + 3)));
  r_e_y.trace_stores();
  Func r_e_;
  r_e_(x, y) = f32(h_sharp_inv_k(0) * r_e_y(x, y) +
		  h_sharp_inv_k(1) * (r_e_y(x - 1, y) + r_e_y(x + 1, y)) +
		  h_sharp_inv_k(2) * (r_e_y(x - 2, y) + r_e_y(x + 2, y)) +
		  h_sharp_inv_k(3) * (r_e_y(x - 3, y) + r_e_y(x + 3, y)));

  r_e_.compute_root();
  //r_e_.trace_stores();
  //Buffer<uint8_t> r_e__b = r_e_.realize(in.width(), in.height());
  //save_image(r_e__b, "../images/chroma_test/diff_gauss_hsharp.jpg");


 
  Expr xi = 9;
  Expr e_ = 2;
  
  RDom left(-(xi + e_), -e_);
  RDom right(e_ + 1, xi + e_ + 1);
  
  Func sum_left;
  sum_left(x, y) += r_e_(x + left, y);

  sum_left.trace_stores();
  Buffer<float> b1 = sum_left.realize(in.width(), in.height());


  Func sum_right;
  sum_right(x, y) = sum(r_e_(x + right, y));

  // sum_right.trace_stores();
  // Buffer<float> b2 = sum_right.realize(in.width(), in.height());


  Func mean_left;
  mean_left(x, y) = sum_left(x, y) / 9.0f;

  // mean_left.trace_stores();
  // Buffer<float> b3 = mean_left.realize(in.width(), in.height());

  Func mean_right;
  mean_right(x, y) = sum_right(x, y) / 9.0f;

  // mean_right.trace_stores();
  // Buffer<float> b4 = mean_right.realize(in.width(), in.height());

  Func sigma_left;
  sigma_left(x, y) = sqrt(sum(pow(r_e_(x + left, y) - mean_left(x, y), 2)) / 9.0f);
  
  // sigma_left.trace_stores();
  // Buffer<float> b5 = sigma_left.realize(in.width(), in.height());
  
  Func sigma_right;
  sigma_right(x, y) = sqrt(sum(pow(r_e_(x + right, y) - mean_right(x, y), 2)) / 9.0f);

  // sigma_right.trace_stores();
  // Buffer<float> b6 = sigma_right.realize(in.width(), in.height());
  
  Func E_ring; //E_ring is the loss function that we are minimizing
  E_ring(x, y) = (sigma_left(x, y) + sigma_right(x, y)) / abs(mean_left(x, y) - mean_right(x, y));

  //   E_ring.trace_stores();

  // Buffer<float> b = E_ring.realize(in.width(), in.height());
  
  RDom r(in);
  Func loss;
  loss() = f64(0.f);
  loss() += E_ring(r.x, r.y);

  loss.trace_stores();

  auto d_loss_d = Halide::propagate_adjoints(loss);

  Param<double> learning_rate;
  learning_rate.set(0.00001);
  Func new_lambda;
  new_lambda(x) = lambda(x) - learning_rate * d_loss_d(lambda)(x); // must only compute at 0
  Func new_sigma;
  new_sigma(x) = sigma(x) - learning_rate * d_loss_d(sigma)(x); // must only compute at 0

  new_lambda.trace_stores();
  new_sigma.trace_stores();

  // // gradient descent loop
  // Pipeline p({loss, new_lambda, new_sigma});


  // sigma_b.fill(0);
  // lambda_b.fill(0);
  // sigma_b(0) = 2.2f;
  // lambda_b(0) = 0.1f;
 

  // auto e = Buffer<double>::make_scalar();
  // const int steps = 100;
  // double initial_error = 0.0;
  // for (int i = 0; i <= steps; i++) {
  //   bool should_print = (i == 0 || i == steps / 2 || i == steps);
  //   if (should_print) {
  //     printf("Iteration %d\n"
  // 	     "sigma: ", i);
  //     printf("%g ", sigma_b(0));
  //     printf("lambda: ");
  //     printf("%g ", lambda_b(0));
  //     printf("\n");
  //   }

  //   //loss.realize(e);
  //   //p.realize({e, lambda_b, sigma_b});
    
  //   if (should_print) {
  //     printf("Err: %g\n", e());
  //   }
    
  //   if (i == 0) {
  //     initial_error = e();
  //   }
  // }
  
  // double final_error = e();
  // if (final_error <= 1e-10 && final_error < initial_error) {
  //   printf("[fit_function] Success!\n");
  //   return 0;
  // } else {
  //   printf("Did not converge\n");
  //   return -1;
  // }
  
} 
