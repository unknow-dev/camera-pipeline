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
  Halide::Buffer<uint8_t> in = load_image("../images/inputs/lena.jpg");
  Func in_bounded = BoundaryConditions::repeat_edge(in);
  
  Expr sigma = 2.2f; //optimize this parameter
  Var x, y, c;
  
  Expr lambda = 0.35f; //optimize this parameter
  Expr gauss = exp(-x*x/(2*sigma*sigma)) / (sqrtf(2*M_PI)*sigma);
  Expr gauss_0 = 1.0f / (sqrtf(2*M_PI)*sigma);
 
  Func h_sharp_inv_k;
  h_sharp_inv_k(x) = - (lambda * -gauss) + (lambda * lambda * gauss * gauss);
  h_sharp_inv_k(0) = 1.0f - lambda * (1.0f - gauss_0) + lambda * lambda * (1.0f - gauss_0) * (1.0f - gauss_0);
												
  Func r_i__y;
  
  // h_sharp_inv_y(x, y, c) = (h_sharp_inv_k(0) * h_sharp(x, y, c) +
  // 		     h_sharp_inv_k(1) * (h_sharp(x, y-1, c) +
  // 				     h_sharp(x, y+1, c)) +
  // 		     h_sharp_inv_k(2) * (h_sharp(x, y-2, c) +
  // 				     h_sharp(x, y+2, c)) +
  // 		      h_sharp_inv_k(3) * (h_sharp(x, y-3, c) +
  // 				     h_sharp(x, y+3, c)));
  r_i__y(x, y) = (h_sharp_inv_k(0) * in_bounded(x, y) +
			    h_sharp_inv_k(1) * (in_bounded(x, y-1) +
						in_bounded(x, y+1)) +
			    h_sharp_inv_k(2) * (in_bounded(x, y-2) +
						in_bounded(x, y+2)) +
			    h_sharp_inv_k(3) * (in_bounded(x, y-3) +
						in_bounded(x, y+3)));
  
  Func r_i_;
  r_i_(x, y) = (h_sharp_inv_k(0) * r_i__y(x, y) +
  		      h_sharp_inv_k(1) * (r_i__y(x-1, y) +
  				      r_i__y(x+1, y)) +
  		      h_sharp_inv_k(2) * (r_i__y(x-2, y) +
  				      r_i__y(x+2, y)) +
  		      h_sharp_inv_k(3) * (r_i__y(x-3, y) +
  				      r_i__y(x+3, y)));
  

  Func r_s; //apply supersampling
  r_s(x, y) = r_i_(x, y);

  Func h_inv_ca_rg_k;
  h_inv_ca_rg_k(x) = 0;
  Func h_inv_ca_bg_k;
  h_inv_ca_bg_k(x) = 0;
  
  Func r_s_y;
  r_s_y(x, y) = u8(h_inv_ca_rg_k(0) * r_s(x, y) +
		   h_inv_ca_rg_k(1) * (r_s(x, y - 1) + r_s(x, y + 1)) +
		   h_inv_ca_rg_k(2) * (r_s(x, y - 2) + r_s(x, y + 2)) +
		   h_inv_ca_rg_k(3) * (r_s(x, y - 3) + r_s(x, y + 3)));
  
  Func r_s_;
  r_s_(x, y) = u8(h_inv_ca_rg_k(0) * r_s_y(x, y) +
		 h_inv_ca_rg_k(1) * (r_s_y(x - 1, y) + r_s_y(x + 1, y)) +
		 h_inv_ca_rg_k(2) * (r_s_y(x - 2, y) + r_s_y(x + 2, y)) +
		 h_inv_ca_rg_k(3) * (r_s_y(x - 3, y) + r_s_y(x + 3, y)));

  float w_x = 1; //pixel width;
  float w_y = 1; //pixel height;
  
  RDom h_box(-2 / w_x, 3 / w_x, -2 / w_y, 3 / w_y);
  r_s_(x, y) = sum(r_s_(x + h_box.x, y, c) + r_s_(x, y + h_box.y, c)) / 25.0f;

  Func r_d_;
  r_d_(x, y) = r_s_(x/3, y/3);

  Expr sigma = 2.2f; //optimize this parameter  
  Expr lambda = 0.35f; //optimize this parameter
  Expr gauss = exp(-x*x/(2*sigma*sigma)) / (sqrtf(2*M_PI)*sigma);
  Expr gauss_0 = 1.0f / (sqrtf(2*M_PI)*sigma);
  Func h_sharp_k;
  h_sharp_k(x) = lambda * -gauss;
  h_sharp_k(0) = 1.0f + lambda * (1.0f - gauss_0);
  
  Func r_i__y;
  r_i__y(x, y) = (h_sharp_k(0) * r_d_(x, y) +
		  h_sharp_k(1) * (r_d_(x, y-1) +
				  r_d_(x, y+1)) +
		  h_sharp_k(2) * (r_d_(x, y-2) +
				  r_d_(x, y+2)) +
		  h_sharp_k(3) * (r_d_(x, y-3) +
				  r_d_(x, y+3)));
  Func r_i__;
  r_i__(x, y) = (h_sharp_k(0) * r_i__y(x, y) +
		  h_sharp_k(1) * (r_i__y(x, y-1) +
				  r_i__y(x, y+1)) +
		 h_sharp_k(2) * (r_i__y(x, y-2) +
				  r_i__y(x, y+2)) +
		 h_sharp_k(3) * (r_i__y(x, y-3) +
				  r_i__y(x, y+3))); 
  

  Func g_i__y;
  g_i__y(x, y) = (h_sharp_inv_k(0) * in_bounded(x, y) +
			    h_sharp_inv_k(1) * (in_bounded(x, y-1) +
						in_bounded(x, y+1)) +
			    h_sharp_inv_k(2) * (in_bounded(x, y-2) +
						in_bounded(x, y+2)) +
			    h_sharp_inv_k(3) * (in_bounded(x, y-3) +
						in_bounded(x, y+3)));
  
  
    
  g_i_(x, y) = (h_sharp_inv_k(0) * g_i__y(x, y) +
  		      h_sharp_inv_k(1) * (g_i__y(x-1, y) +
  				      g_i__y(x+1, y)) +
  		      h_sharp_inv_k(2) * (g_i__y(x-2, y) +
  				      g_i__y(x+2, y)) +
  		      h_sharp_inv_k(3) * (g_i__y(x-3, y) +
  				      g_i__y(x+3, y)));
  
  Func g_s_;
  g_s_(x, y) = g_i_(x, y);
  Buffer<uint8_t> g_s_b = g_s_.realize(0, in.width(), 0, in.height());

  Func b_i_;
  b_i_(x, y) = 
  b_i_(x, y) = (h_sharp_inv_k(0) * b_i__y(x, y) +
  		      h_sharp_inv_k(1) * (b_i__y(x-1, y) +
  				      b_i__y(x+1, y)) +
  		      h_sharp_inv_k(2) * (b_i__y(x-2, y) +
  				      b_i__y(x+2, y)) +
  		      h_sharp_inv_k(3) * (b_i__y(x-3, y) +
  				      b_i__y(x+3, y)));

  Func b_s;
  b_s(x, y) = b_i_(x, y);

  Func b_s_;
  b_s_(x, y) = u8(h_inv_ca_rg_k(0) * b_s_y(x, y) +
		 h_inv_ca_rg_k(1) * (b_s_y(x - 1, y) + b_s_y(x + 1, y)) +
		 h_inv_ca_rg_k(2) * (b_s_y(x - 2, y) + b_s_y(x + 2, y)) +
		 h_inv_ca_rg_k(3) * (b_s_y(x - 3, y) + b_s_y(x + 3, y)));

  Func b_d_;
  b_d_(x, y) = b_s_(x/3, y/3);
 
  Func b_i__y;
  b_i__y(x, y) = (h_sharp_k(0) * b_d_(x, y) +
		  h_sharp_k(1) * (b_d_(x, y-1) +
				  b_d_(x, y+1)) +
		  h_sharp_k(2) * (b_d_(x, y-2) +
				  b_d_(x, y+2)) +
		  h_sharp_k(3) * (b_d_(x, y-3) +
				  b_d_(x, y+3)));
  Func b_i__;
  b_i__(x, y) = (h_sharp_k(0) * b_i__y(x, y) +
		  h_sharp_k(1) * (b_i__y(x, y-1) +
				  b_i__y(x, y+1)) +
		 h_sharp_k(2) * (b_i__y(x, y-2) +
				  b_i__y(x, y+2)) +
		 h_sharp_k(3) * (b_i__y(x, y-3) +
				  b_i__y(x, y+3))); 

    
  // auto d = propogate_adjoints(r_s_, adjoints);
  // Func d_in = d(in);
  // RDom r(g_s_b);
  // Func loss_r;
  // loss_r() = 0.f;
  // Expr diff_r = r_s_(r.x, r.y) - g_s_b(r.x, r.y);
  // loss_r() += diff_r * diff_r;

  // Func loss_b;
  // loss_b() = 0.f;
  // Expr diff_b = b_s_(r.x, r.y) - g_s_b(r.x, r.y);
  // loss_b() += diff_b * diff_b;
  
  // auto d_loss_r_d = propagate_adjoints(loss_r);
  // Func d_loss_r_d_h_inv_ca_rg = d_loss_r_d(h_inv_ca_rg);

  // auto d_loss_b_d = propagate_adjoingts(loss_b);
  // Func d_loss_b_d_h_inv_ca_bg = d_loss_b_d(h_inv_ca_bg);
  
  Func img;
  img(x, y, c) = u8(0);
  img(x, y, R) = r_i__(x, y);
  img(x, y, G) = in_bounded(x, y, G);
  img(x, y, B) = b_i__(x, y);
  Buffer<uint8_t> output = img.realize(in.width(), in.height()); 
  save_image(output, "../images/outputs/chroma_corrected.jpg");

  // Func p; // this is the supersampling function 
  // p(x, y, c) = h_sharp_inv(x, y, 0);

  // Func p_red;
  // p_red(x, y) = p(x, y, R);

  // Func p_green;
  // p_green(x, y) = p(x, y, G);
  
  // Func p_blue;
  // p_blue(x, y) = p(x, y, B);
  
  

  // Func p_red_conv_h_inv_y;
  // p_red_conv_h_inv_y(x, y) = (h_inv_ca_rg_k(0) * p_red(x, y) +
  // 			      h_inv_ca_rg_k(1) * (p_red(x, y-1) +
  // 					 p_red(x, y+1)) +
  // 			      h_inv_ca_rg_k(2) * (p_red(x, y-2) +
  // 					 p_red(x, y+2)) +
  // 			      h_inv_ca_rg_k(3) * (p_red(x, y-3) +
  // 					  p_red(x, y+3)));
  // Func p_red_conv_h_inv;
  // p_red_conv_h_inv(x, y) = (h_sharp_inv_k(0) * p_red_conv_h_inv_y(x, y) +
  // 		      h_sharp_inv_k(1) * (p_red_conv_h_inv_y(x - 1, y) +
  // 					  p_red_conv_h_inv_y(x + 1, y)) +
  // 		      h_sharp_inv_k(2) * (p_red_conv_h_inv_y(x - 2, y) +
  // 					  p_red_conv_h_inv_y(x + 2, y)) +
  // 		      h_sharp_inv_k(3) * (p_red_conv_h_inv_y(x - 3, y) +
  // 					  p_red_conv_h_inv_y(x + 3, y)));

  // Func loss_red;
  // loss_red(x, y) =  pow(abs(p_red_conv_h_inv(x, y) - p_green(x, y)), 2);

  // Func p_blue_conv_h_inv_y;
  // p_blue_conv_h_inv_y(x, y) = (h_inv_ca_bg_k(0) * p_blue(x, y) +
  // 			 h_inv_ca_bg_k(1) * (p_blue(x, y-1) +
  // 					     p_blue(x, y+1)) +
  // 			 h_inv_ca_bg_k(2) * (p_blue(x, y-2) +
  // 					     p_blue(x, y+2)) +
  // 			 h_inv_ca_bg_k(3) * (p_blue(x, y-3) + p_blue(x, y+3))); 

  // Func p_blue_conv_h_inv;
  // p_blue_conv_h_inv(x, y) = (h_inv_ca_bg_k(0) * p_blue(x, y) +
  // 			 h_inv_ca_bg_k(1) * (p_blue(x - 1, y) +
  // 					     p_blue(x + 1, y)) +
  // 			 h_inv_ca_bg_k(2) * (p_blue(x - 2, y) +
  // 					     p_blue(x + 2, y)) +
  // 			 h_inv_ca_bg_k(3) * (p_blue(x - 3, y) + p_blue(x + 3, y))); 

  // Func loss_blue;
  // loss_blue(x, y) = pow(abs(p_blue_conv_h_inv(x, y) - p_green(x, y)), 2);
  
  // // Func loss;
  // // loss() = 0.f;
  // // loss() += loss_red(x, y) + loss_blue(x, y);
  // // auto d_loss_d = Halide::propagate_adjoints(loss);
  // // Func d_loss_d_h_inv_ca_rg_k = d_loss_d(h_inv_ca_rg_k);
  // // Func d_loss_d_h_inv_ca_bg_k = d_loss_d(h_inv_ca_bg_k);

  // Expr x_r = pow(in.width() - x, 2);
  // Expr y_r = pow(in.height() - y, 2);
    

  // float w_x = 1; //pixel width
  // float w_y = 1; //pixel height
  // Func cont_img_conv_h_box;
  // RDom h_box(-2 / w_x, 3 / w_x, -2 / w_y, 3 / w_y);
  // cont_img_conv_h_box(x, y, c) = sum(in_bounded(x + h_box.x, y, c) + in_bounded(x, y + h_box.y, c)) / 25.0f;
  
  // Func img_d;
  // img_d(x, y, c) = u8(cont_img_conv_h_box(3 * x, 3 * y,  c));

  
  // Func img;
  // img(x, y, c) = f32(0);
  // img(x, y, R) = (h_sharp_k(0) * img_d(x, y, R) +
  // 		  h_sharp_k(1) * (img_d(x-1, y, R) +
  // 				  img_d(x+1, y, R)) +
  // 		      h_sharp_k(2) * (img_d(x-2, y, R) +
  // 				      img_d(x+2, y, R)) +
  // 		  h_sharp_k(3) * (img_d(x-2, y, R) +
  // 				  img_d(x+2, y, R)));
  
  // img(x, y, B) = (h_sharp_k(0) * img_d(x, y, B) +
  // 		  h_sharp_k(1) * (img_d(x-1, y, B) +
  // 				  img_d(x+1, y, B)) +
  // 		      h_sharp_k(2) * (img_d(x-2, y, B) +
  // 				      img_d(x+2, y, B)) +
  // 		  h_sharp_k(3) * img_d(x-2, y, B) +
  // 		  img_d(x+2, y, B));
  
  // //img(x, y, G) = f32(in_bounded(x, y, G));
  // img(x, y, G) = f32(0);

  // Func img_u8;
  // img_u8(x, y, c) = u8(img(x, y , c));
  
  // double t1 = current_time(); 
  // Buffer<uint8_t> output = img_u8.realize(in.width(), in.height(), in.channels());
  // double t2 = current_time();
  // cout<<"Time: "<< t2 - t1 <<endl;
  // save_image(output, "../images/outputs/conv_flower.jpg"); 
  
} 



