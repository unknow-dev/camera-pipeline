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
  Halide::Buffer<uint8_t> in = load_image("../images/chroma_test/in256.jpg");
  Func m_i = BoundaryConditions::repeat_edge(in);
  
  Var x, y, c;
  Expr sigma = 2.2f; //optimize this parameter
  Expr lambda = 0.35f; //optimize this parameter
  Expr gauss = exp(-x*x/(2*sigma*sigma)) / (sqrtf(2*M_PI)*sigma);
  Expr gauss_0 = 1.0f / (sqrtf(2*M_PI)*sigma);
  
  Func h_sharp_inv_k;
  h_sharp_inv_k(x) = - (lambda * -gauss) + (lambda * lambda * gauss * gauss);
  h_sharp_inv_k(0) = 1.0f - lambda * (1.0f - gauss_0) + lambda * lambda * (1.0f - gauss_0) * (1.0f - gauss_0);
  

  Func m_i_y;
  m_i_y(x, y, c) = (h_sharp_inv_k(0) * m_i(x, y, c) +
  		    h_sharp_inv_k(1) * (m_i(x, y - 1, c) + m_i(x, y + 1, c)) +
  		    h_sharp_inv_k(2) * (m_i(x, y - 2, c) + m_i(x, y + 2, c)) +
  		    h_sharp_inv_k(3) * (m_i(x, y - 3, c) + m_i(x, y + 3, c)));
  
  Func m_i_;
  m_i_(x, y, c) = u8(h_sharp_inv_k(0) * m_i_y(x, y, c) +
  		h_sharp_inv_k(1) * (m_i_y(x - 1, y, c) + m_i_y(x + 1, y, c)) +
  		h_sharp_inv_k(2) * (m_i_y(x - 2, y, c) + m_i_y(x + 2, y, c)) +
  		h_sharp_inv_k(3) * (m_i_y(x - 3, y, c) + m_i_y(x + 3, y, c)));

  
  // Buffer<uint8_t> mi_b = m_i_.realize(in.width(), in.height(), in.channels());
  // save_image(mi_b, "../images/chroma_test/mi_b.jpg");

  // Buffer<uint8_t> mib = m_i.realize(in.width(), in.height(), in.channels());
  // save_image(mib, "../images/chroma_test/mib.jpg");

  Var j_x, j_y;
  
  Func f;
  f(x, y, c, j_x, j_y) = m_i_(x + j_x, y + j_y, c);
 
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
  
  Func m_s;
  m_s(x, y, c) = u8(min((b(x/3, y, c, -1) * t_n_1_y + b(x/3, y, c, 0) * t_0_y + b(x/3, y, c, 1) * t_1_y + b(x/3, y, c, 2) * t_2_y) / 2.0f, 255.f));
  
  Func img;
  img(x, y, c) = u8(0);
  img(x, y, G) = m_i(x, y, G);
  
  // Buffer<uint8_t> m_s_b = m_s.realize(3 * in.width(), 3 * in.height(), in.channels());
  //m_s.trace_stores();
  m_s.compute_root();
  
  
  ImageParam h_inv_ca_rg_k(Float(64), 2);
  
  int kernel_size = 7;
  
  Param<double> learning_rate;
  
  RDom r_rg(0, kernel_size, 0, kernel_size);
  
  Func r_s_;
  r_s_(x, y) += f64(h_inv_ca_rg_k(r_rg.x, r_rg.y) * m_s(x + r_rg.x - (kernel_size/2), y + r_rg.y - (kernel_size/2), R) / (kernel_size * kernel_size));
  //r_s_.trace_stores();
  float w_x = 1; //pixel width;
  float w_y = 1; //pixel height;
  
  RDom h_box(-2 / w_x, 3 / w_x, -2 / w_y, 3 / w_y);
  
  Func r_d_unsamp;
  r_d_unsamp(x, y) += u8((r_s_(x + h_box.x, y + h_box.y))/5.0f); // box
  
  //r_d_unsamp.trace_stores();
  
  Func r_d_;
  r_d_(x, y) = u8(r_d_unsamp(x/3, y/3)); // sample
  r_d_.compute_root();
  
  Func h_sharp_k;
  h_sharp_k(x) = lambda * -gauss;
  h_sharp_k(0) = 1.0f + lambda * (1.0f - gauss_0);
  

  Func r_i__y;
  r_i__y(x, y) = (h_sharp_k(0) * r_d_(x, y) +
  		  h_sharp_k(1) * (r_d_(x, y - 1) + r_d_(x, y + 1)) +
  		  h_sharp_k(2) * (r_d_(x, y - 2) + r_d_(x, y + 2)) +
  		  h_sharp_k(3) * (r_d_(x, y - 3) + r_d_(x, y + 3)));					        		
  Func r_i__;
  r_i__(x, y) = u8(h_sharp_k(0) * r_i__y(x, y) +
  		 h_sharp_k(1) * (r_i__y(x, y - 1) + r_i__y(x, y + 1)) +
  		 h_sharp_k(2) * (r_i__y(x, y - 2) + r_i__y(x, y + 2)) +
  		 h_sharp_k(3) * (r_i__y(x, y - 3) + r_i__y(x, y + 3)));

  // Auto Diff
  RDom r(in);
  Func loss_r;
  loss_r() = f64(0.f);
  Expr diff_r = r_s_(r.x, r.y) - m_s(r.x, r.y, G);
  loss_r() += diff_r * diff_r;
  auto d_loss_r_d = propagate_adjoints(loss_r);
  
  Func new_r_k;

  new_r_k(x, y) = h_inv_ca_rg_k(x, y) - learning_rate * d_loss_r_d(h_inv_ca_rg_k)(x, y);

  
  Buffer<double> h_inv_ca_rg_k_b(kernel_size, kernel_size);
  h_inv_ca_rg_k_b.fill(0);
  
  for (int i = 0; i < kernel_size; i++) {
    for (int j = 0; j < kernel_size; j++) {
      h_inv_ca_rg_k_b(i, j) = (double) rand() / RAND_MAX;
    }
  }

  h_inv_ca_rg_k.set(h_inv_ca_rg_k_b);

  // Buffer<uint8_t> r_s_b = r_s_.realize(in.width() * 3, in.height() * 3);
  // save_image(r_s_b, "../images/chroma_test/512rsb.jpg");

  // Func msr;
  // msr(x, y) = m_s(x, y, R);
  // Buffer<uint8_t> m_s_b = msr.realize(in.width() * 3, in.height() * 3);
  // save_image(m_s_b, "../images/chroma_test/512msb.jpg");
  
  // Buffer<double> r_s_b = r_s_.realize(in.width(), in.height());
  // Buffer<float> m_s_b = m_s.realize(in.width(), in.height(), in.channels());

  // for (int i = 0; i < in.width(); i++) {
  //   for (int j = 0; j < in.height(); j++) {
  //     printf("r_s_b %d %d %f\n", i, j, r_s_b(i, j));
  //     printf("m_s_b %d %d %f\n", i, j, m_s_b(i, j));
  //   }
  // }

  auto e = Buffer<double>::make_scalar();
  Pipeline p({loss_r, new_r_k});
  
  const int steps = 1;
  double initial_error = 0.0;
  
  for (int i = 0; i <= steps; i++) {
    
    bool should_print = (i == 0|| i == steps/2 || i == steps);
    if(should_print) {
      printf("Iteration %d\n"
  	     "Coefficients: ", i);		
      for (int j = 0; j < kernel_size; j++) {
  	for(int k = 0; k < kernel_size; k++) {
  	  printf("%g ", h_inv_ca_rg_k_b(j, k));
  	}
  	printf("\n");
      }
    }
    //double t1 = current_time(); 
    
    p.realize({e, h_inv_ca_rg_k_b});
    
    //double t2 = current_time();
    //cout<<"Time: "<< t2 - t1 <<endl;
    
    if (should_print) {
      printf("Err: %g\n", e());
    }
    
    if (i == 0) {
      initial_error = e();
    } 
  }
  
  double final_error = e();
  
  if (final_error <= 1e-10 && final_error < initial_error) {
    printf("[fit_function] Success!\n");
    //return 0;
  } else {
    printf("Did not converge\n");
    //return -1;
  }

  Buffer<uint8_t> r_d_b = r_d_.realize(in.width(), in.height());
  save_image(r_d_b, "../images/chroma_test/rd.jpg");
  
  Buffer<uint8_t> r_i__b = r_i__.realize(in.width(), in.height());
  save_image(r_i__b, "../images/chroma_test/out256.jpg");

  Func r_s_u8;
  r_s_u8(x, y) = u8(r_s_(x, y));
  
  Buffer<uint8_t> r_s_b = r_s_u8.realize(in.width(), in.height());


  Func g_s;
  g_s(x, y) = u8(m_s(x, y, G));
  Buffer<uint8_t> tgt = g_s.realize(in.width(), in.height());

  save_image(r_s_b, "../images/chroma_test/r_s_.jpg");

  save_image(tgt, "../images/chroma_test/tgt.jpg");
  // Func grad;
  // grad(x, y) = f64(d_loss_r_d(h_inv_ca_rg_k)(x, y));
  // grad.trace_stores();
  // Buffer<double> grad_b = grad.realize(kernel_size, kernel_size);
  //save_image(grad_b, "../images/chroma_test/grad.jpg");
}


