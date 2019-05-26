#include <stdio.h>
#include <math.h>

#include "Halide.h"
#include "halide_image_io.h"
//#include "denoise.h"
#include "clock.h"
#define R 0
#define G 1
#define B 2


using namespace Halide;
using namespace Halide::ConciseCasts;
using namespace std;


int main(int argc, char **argv) {

	Buffer<uint8_t> in_b = Tools::load_image("../images/noise_test/bike.jpg");
	printf("channels: %d\n", in_b.channels());
	Buffer<uint8_t> tgt = Tools::load_image("../images/noise_test/bike.jpg");
	Func in = BoundaryConditions::repeat_edge(in_b);
	
	float sig_s;
	float sig_r;
	sig_s = 2.22f; //tune this
	sig_r = 2.22f; //tune thisb
	
        Expr L = i32(sig_r * 3.0f);
	Expr W = i32(sig_s * 3.0f);
        
	int W_ = 2.22 * 3.0f;
	
	Var x, y, c, t;

	// declare range and spatial filters
	
	RDom omega(0, 2 * W, 0, 2 * W), l_r(-L, L + 1, -L, L + 1);

	ImageParam g_sig_s(Float(64), 3);
        ImageParam g_sig_r(Float(64), 1);
		 
	Func box_filtered;
	box_filtered(x, y, c) = u8((float)(1) / ((2 * L + 1) * (2 * L + 1)) * sum(in(x + l_r.x, y + l_r.y, c)));

	Func denoised_num;
	denoised_num(x, y, c) += f32(g_sig_s(omega.x, omega.y, c) * g_sig_r(abs(box_filtered(x + omega.x - W, y + omega.y - W, c) - box_filtered(x, y, c))) * in(x + omega.x - W, y + omega.y - W, c));

	//denoised_num.trace_stores();
	Func denoised_den;
	denoised_den(x, y, c) += f32((g_sig_s(omega.x, omega.y, c)) * g_sig_r(abs(box_filtered(x + omega.x - W, y + omega.y - W, c) - box_filtered(x, y, c))));
	
	Func denoised;
	denoised(x, y, c) = f32(denoised_num(x, y, c) / denoised_den(x, y, c));
	denoised.trace_stores();
	
	//Buffer<uint8_t> output = denoised.realize(in_b.width(), in_b.height(), in_b.channels());
	
	RDom r(tgt);
	Func loss;
	loss() = f64(1.f);
	Expr diff_g = denoised(r.x, r.y, r.z) - tgt(r.x, r.y, r.z);
	loss() = f64(sum(diff_g * diff_g));
       
	loss.print_loop_nest();

	auto d_loss_d = propagate_adjoints(loss);
	
	auto e = Buffer<double>::make_scalar();

	Param<double> learning_rate;
	learning_rate.set(0.00001);
	
	Func new_g_sig_s;
	new_g_sig_s(x, y, c) = g_sig_s(x, y, c) - learning_rate * d_loss_d(g_sig_s)(x, y, c);
	//new_g_sig_s(x, y) = 0;

	Func new_g_sig_r;
	new_g_sig_r(t) = g_sig_r(t) - learning_rate * d_loss_d(g_sig_r)(t);
	
	Buffer<double> g_sig_r_b(256);
	g_sig_r.set(g_sig_r_b);

	
	for(int i = 0; i < 256; i++) {
	  g_sig_r_b(i) = exp(- i * i / (2 * sig_r * sig_r));
	  //printf("These values: %g\n", g_sig_r_b(i));
	}
	
	const int steps = 10;
	double initial_error = 0.0;
		
	Pipeline p({loss, new_g_sig_s, new_g_sig_r});
	
	//scheduling
	Var v;
	Func fs[] = {g_sig_s, loss};
	for (Func f : fs) {
	  
	  bool first = true;
	  for (Func df: d_loss_d.funcs(f)) {
	    if (first) {
	      first = false;
	      continue;
	    }
	    //df.compute_root().vectorize(df.args()[0], 4);
	    for (int i = 0; i < df.num_update_definitions(); i++) {
	      for (auto d: df.update(i).get_schedule().dims()) {
		if (d.is_pure()) {
		  df.update(i).vectorize(Var(d.var), 4);
		  break;
		}
	      }
	    }
	  }
	}
	
	// Gradient Descent Loop
	Buffer<double> g_sig_s_b(2 * W_, 2 * W_, 3);
	
	g_sig_s.set(g_sig_s_b);
        
	g_sig_s_b.fill(0);
	for (int i = 0; i <= 2 * W_; i++) {
	  for(int j = 0; j <= 2 * W_; j++) {
	    int i_ = i - W_;
	    int j_ = j - W_;
	    g_sig_s_b(i, j) = exp(-(i_ * i_ + j_ * j_)/(2 * sig_s * sig_s)); 
	    
	  }
	}
	
	for (int i = 0; i <= steps; i++) {
	  bool should_print = (i == 0 || i == steps/2 || i == steps);
	  if (should_print) {
	    printf("Iteration %d\n"
		   "Coefficients g_sig_s: ", i);
	    for (int j = 0; j <= W_; j++) {
	      for (int k = 0; k <= W_; k++) {
		printf("%g ", g_sig_s_b(j, k));
	      }
	    }
	    
	    printf("Iteration %d\n"
		   "Coefficients g_sig_r: ", i);
	    
	    for (int j = 0; j < 256; j++) {
	      printf("%g ", g_sig_r_b(j));
	    }
	    
	    printf("\n");
	  }
	  
	  double t1 = current_time();
	  
	  p.realize({e, g_sig_s_b, g_sig_r_b});
	  //loss_g.realize(e);
	  double t2 = current_time();
	  
	  cout<<"Time: "<< t2 - t1<<endl;
	  //Tools::save_image(output, "../images/noise_test/out300.jpg");
	  
	  if (should_print) {
	    printf("Err: %g\n", e());
	  }
	  
	  if (i == 0) {
	    initial_error = e();
	  }
	}
	
	double final_error = e();
	if (final_error<= 1e-10 && final_error < initial_error) {
	  printf("[fit function] Success!\n");
	  //return 0;
	} else {
	  printf("Did not converge!\n");
	  //return -1;
	}

	Func derv;
	derv(x, y, c) = (d_loss_d(g_sig_s)(x, y, c));
	derv.trace_stores();
	Buffer<double> b = derv.realize(in_b.width(), in_b.height(), in_b.channels());
	
}
