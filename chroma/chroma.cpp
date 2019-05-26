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
  m_s(x, y, c) = u8((b(x/3, y, c, -1) * t_n_1_y + b(x/3, y, c, 0) * t_0_y + b(x/3, y, c, 1) * t_1_y + b(x/3, y, c, 2) * t_2_y) / 2.0f);
  
  Func img;
  img(x, y, c) = u8(0);
  img(x, y, G) = m_i(x, y, G);

  // Buffer<uint8_t> m_s_b = m_s.realize(3 * in.width(), 3 * in.height(), in.channels());
  //m_s.trace_stores();
  m_s.compute_root();

  // Func h_inv_ca_rg_k_y;
  // h_inv_ca_rg_k_y(x) = 1;

  // Func h_inv_ca_rg_k_x;
  // h_inv_ca_rg_k_x(x) = 1;

  Func h_inv_ca_rg_k;
  h_inv_ca_rg_k(x, y) = 1;

  RDom r_rg(-3, 4, -3, 4);
  
  Func r_s_;
  r_s_(x, y) += h_inv_ca_rg_k(r_rg.x, r_rg.y) * m_s(x + r_rg.x, y + r_rg.y, R);
  
  
  // Func r_s_y;
  // r_s_y(x, y) = u8(h_inv_ca_rg_k_y(0) * m_s(x, y, R) +
  // 		   h_inv_ca_rg_k_y(1) * (m_s(x, y - 1, R) + m_s(x, y + 1, R)) +
  // 		   h_inv_ca_rg_k_y(2) * (m_s(x, y - 2, R) + m_s(x, y + 2, R)) +
  // 		   h_inv_ca_rg_k_y(3) * (m_s(x, y - 3, R) + m_s(x, y + 3, R)));

  
  // Func r_s_;
  // r_s_(x, y) = u8(h_inv_ca_rg_k_x(0) * r_s_y(x, y) +
  // 		 h_inv_ca_rg_k_x(1) * (r_s_y(x - 1, y) + r_s_y(x + 1, y)) +
  // 		 h_inv_ca_rg_k_x(2) * (r_s_y(x - 2, y) + r_s_y(x + 2, y)) +
  // 		 h_inv_ca_rg_k_x(3) * (r_s_y(x - 3, y) + r_s_y(x + 3, y)));
  //r_s_(x, y) = m_s(x, y, R);
  //r_s_.trace_stores();
  float w_x = 1; //pixel width;
  float w_y = 1; //pixel height;
  
  RDom h_box(-2 / w_x, 3 / w_x, -2 / w_y, 3 / w_y);

  Func r_d_unsamp;
  r_d_unsamp(x, y) = u8((r_s_(x - 2, y - 2) + r_s_(x - 1, y - 2) + r_s_(x, y - 2) + r_s_(x + 1, y - 2) + r_s_(x + 2, y - 2) +
			r_s_(x - 2, y - 1) + r_s_(x - 1, y - 1) + r_s_(x, y - 1) + r_s_(x + 1, y - 1) + r_s_(x + 2, y - 1) +
			r_s_(x - 2, y) + r_s_(x - 1, y) + r_s_(x, y) + r_s_(x + 1, y) + r_s_(x + 2, y) +
			r_s_(x - 2, y + 1) + r_s_(x - 1, y + 1) + r_s_(x, y + 1) + r_s_(x + 1, y + 1) + r_s_(x + 2, y + 1) +
			r_s_(x - 2, y + 2) + r_s_(x - 1, y + 2) + r_s_(x, y + 2) + r_s_(x + 1, y + 2) + r_s_(x + 2, y + 2))/5.0f); // box
  
  r_d_unsamp.trace_stores();

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
  loss_r() = 0.f;
  Expr diff_r = r_s_(r.x, r.y) - m_s(r.x, r.y, G);
  loss_r() += diff_r * diff_r;
 
  auto d_loss_r_d = propagate_adjoints(loss_r);
  Func d_loss_r_d_h_inv_ca_rg_k = d_loss_r_d(h_inv_ca_rg_k);
  h_inv_ca_rg_k.trace_stores();

  img(x, y, R) = r_i__(x, y);

  ////Blue Channel////

  Func h_inv_ca_bg_k_y;
  h_inv_ca_bg_k_y(x) = 1;

  Func h_inv_ca_bg_k;
  h_inv_ca_bg_k(x, y) = 1;

  RDom r_bg(-3, 4, -3, 4);
  Func b_s_;
  b_s_(x, y) += h_inv_ca_bg_k(r_bg.x, r_bg.y) * m_s(x + r_bg.x, y + r_bg.y, B);

  
  
  
  // Func b_s_y;
  // b_s_y(x, y) = u8(h_inv_ca_bg_k_y(0) * m_s(x, y, B) +
  // 		   h_inv_ca_bg_k_y(1) * (m_s(x, y - 1, B) + m_s(x, y + 1, B)) +
  // 		   h_inv_ca_bg_k_y(2) * (m_s(x, y - 2, B) + m_s(x, y + 2, B)) +
  // 		   h_inv_ca_bg_k_y(3) * (m_s(x, y - 3, B) + m_s(x, y + 3, B)));

  // Func b_s_;
  // b_s_(x, y) = u8(h_inv_ca_bg_k_x(0) * b_s_y(x, y) +
  // 		 h_inv_ca_bg_k_x(1) * (b_s_y(x - 1, y) + b_s_y(x + 1, y)) +
  // 		 h_inv_ca_bg_k_x(2) * (b_s_y(x - 2, y) + b_s_y(x + 2, y)) +
  // 		 h_inv_ca_bg_k_x(3) * (b_s_y(x - 3, y) + b_s_y(x + 3, y)));

  Func b_d_unsamp;
  b_d_unsamp(x, y) = sum(b_s_(x + h_box.x, y) + b_s_(x, y + h_box.y)) / 25.0f;

  Func b_d_;
  b_d_(x, y) = b_d_unsamp(x/3, y/3); // box and sample
 
  Func b_i__y;
  b_i__y(x, y) = (h_sharp_k(0) * b_d_(x, y) +
  		  h_sharp_k(1) * (b_d_(x, y-1) +
  				  b_d_(x, y+1)) +
  		  h_sharp_k(2) * (b_d_(x, y-2) +
  				  b_d_(x, y+2)) +
  		  h_sharp_k(3) * (b_d_(x, y-3) +
  				  b_d_(x, y+3)));
  Func b_i__;
  b_i__(x, y) = u8(h_sharp_k(0) * b_i__y(x, y) +
  		  h_sharp_k(1) * (b_i__y(x, y-1) +
  				  b_i__y(x, y+1)) +
  		 h_sharp_k(2) * (b_i__y(x, y-2) +
  				  b_i__y(x, y+2)) +
  		 h_sharp_k(3) * (b_i__y(x, y-3) +
  				  b_i__y(x, y+3))); 

  // Func loss_b;
  // loss_b() = 0.f;
  // Expr diff_b = b_s_(r.x, r.y) - m_s(r.x, r.y, G);
  // loss_b() += diff_b * diff_b;
  
  
  // auto d_loss_b_d = propagate_adjoints(loss_b);
  // Func d_loss_b_d_h_inv_ca_bg_k = d_loss_b_d(h_inv_ca_bg_k);
  img(x, y, B) = b_i__(x, y);

  r_i__.trace_stores();
  b_i__.trace_stores();
  double t1 = current_time(); 
  Buffer<uint8_t> out = img.realize(in.width(), in.height(), in.channels()); 
  double t2 = current_time();
  cout<<"Time: "<< t2 - t1 <<endl;
  save_image(out, "../images/chroma_test/out_.jpg");  

  
  
} 
