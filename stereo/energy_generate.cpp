//
// Created by Prikshet on 8/14/2019.
//

#include <stdio.h>
#include <string>
#include "Halide.h"
#include <math.h>
#include "clock.h"
#include <vector>


using namespace std;
using namespace Halide::ConciseCasts;
using namespace Halide;
//using namespace Halide::Tools;

//
//
//// A helper function to check if OpenCL, Metal or D3D12 is present on the host machine.
//
//#ifdef _WIN32
//#include <windows.h>
//#else
//#include <dlfcn.h>
//#endif
//
//Target find_gpu_target() {
//    // Start with a target suitable for the machine you're running this on.
//    Target target = get_host_target();
//
//    // Uncomment the following lines to try CUDA instead:
//    // target.set_feature(Target::CUDA);
//    // return target;
//
//#ifdef _WIN32
//    if (LoadLibraryA("d3d12.dll") != nullptr) {
//        target.set_feature(Target::D3D12Compute);
//    } else if (LoadLibraryA("OpenCL.dll") != nullptr) {
//        target.set_feature(Target::OpenCL);
//    }
//#elif __APPLE__
//    // OS X doesn't update its OpenCL drivers, so they tend to be broken.
//    // CUDA would also be a fine choice on machines with NVidia GPUs.
//    if (dlopen("/System/Library/Frameworks/Metal.framework/Versions/Current/Metal", RTLD_LAZY) != NULL) {
//        target.set_feature(Target::Metal);
//    }
//#else
//    if (dlopen("libOpenCL.so", RTLD_LAZY) != NULL) {
//        target.set_feature(Target::OpenCL);
//    }
//#endif
//
//    return target;
//}
//
//void tune_hyperparams() {
//
//    Param<int> order, samples;
//    Param<double> learning_rate;
//    Var x("x"), y("y");
//    ImageParam coeffs(Float(64), 1);
//
//
//    const int terms = 8;
//
//    Param<float> sigma;
//    Param<float> gamma;
//    Param<float> tau;
//
//
//    ImageParam ground_truth(Float(64), 2);
//    ImageParam result(Float(64), 2);
//
//    Func err;
//
//    err(x, y) = pow((result(x, y) - ground_truth(x, y)) / ground_truth(x, y), 2);
//
//    RDom d(1, samples - 1);
//    Func average_err;
//    average_err() = sum(err(d)) / samples;
//
//    auto d_err_d = propagate_adjoints(average_err);
//
//    // Compute the new coefficients in terms of the old.
//
//    Func new_coeffs;
//    new_coeffs(x) = coeffs(x) - learning_rate * d_err_d(coeffs)(x);
//
//
//    Buffer<double> c(terms);
//    order.set(terms);
//    samples.set(1000);
//    auto e = Buffer<double>::make_scalar();
//    coeffs.set(c);
//    Pipeline p({average_err, new_coeffs});
//    c.fill(0);
//    // Initialize to the Taylor series for sin about zero
//    c(0) = 1;
//    for (int i = 1; i < terms; i++) {
//        c(i) = -c(i-1)/(i*2*(i*2 + 1));
//    }
//    learning_rate.set(0.00001);
//    const int steps = 10000;
//    double initial_error = 0.0;
//    for (int i = 0; i <= steps; i++) {
//        bool should_print = (i == 0 || i == steps/2 || i == steps);
//        if (should_print) {
//            printf("Iteration %d\n"
//                   "Coefficients: ", i);
//            for (int j = 0; j < terms; j++) {
//                printf("%g ", c(j));
//            }
//            printf("\n");
//        }
//
//        p.compile_to_static_library("hyperparams", {sigma, gamma, tau, ground_truth, result}, "hyperparams");
//
//        if (should_print) {
//            printf("Err: %g\n", e());
//        }
//
//        if (i == 0) {
//            initial_error = e();
//        }
//    }
//}
//
//
//
//bool schedule_for_gpu(Func &f1, Func &f2, Func &f3) {
//    Var x, y, c, xo, yo, xi, yi;
//    Target target = fifitnd_gpu_target();
//    if (!target.has_gpu_feature()) {
//        return false;
//    }
//
//    // If you want to see all of the OpenCL, Metal, CUDA or D3D 12 API
//    // calls done by the pipeline, you can also enable the Debug flag.
//    // This is helpful for figuring out which stages are slow, or when
//    // CPU -> GPU copies happen. It hurts performance though, so we'll
//    // leave it commented out.
//    //target.set_feature(Target::Debug);
//
//    // We make the decision about whether to use the GPU for each
//    // Func independently. If you have one Func computed on the
//    // CPU, and the next computed on the GPU, Halide will do the
//    // copy-to-gpu under the hood. For this pipeline, there's no
//    // reason to use the CPU for any of the stages. Halide will
//    // copy the input image to the GPU the first time we run the
//    // pipeline, and leave it there to reuse on subsequent runs.
//
//    // As before, we'll compute the LUT once at the start of the
//    // pipeline.
//    f1.compute_root();
//
//
//
//    // This is a very common scheduling pattern on the GPU, so
//    // there's a shorthand for it:
//
//    // lut.gpu_tile(i, block, thread, 16);
//
//    // Func::gpu_tile behaves the same as Func::tile, except that
//    // it also specifies that the tile coordinates correspond to
//    // GPU blocks, and the coordinates within each tile correspond
//    // to GPU threads.
//
//    // Compute color channels innermost. Promise that there will
//    // be three of them and unroll across them.
//    f2.reorder(c, x, y)
//            .bound(c, 0, 3)
//            .unroll(c);
//
//    // Compute curved in 2D 8x8 tiles using the GPU.
//    f2.gpu_tile(x, y, xo, yo, xi, yi, 8, 8);
//
//    // This is equivalent to:
//    // curved.tile(x, y, xo, yo, xi, yi, 8, 8)
//    //       .gpu_blocks(xo, yo)
//    //       .gpu_threads(xi, yi);
//
//    // We'll leave sharpen as inlined into curved.
//
//    // Compute the padded input as needed per GPU block, storing
//    // the intermediate result in shared memory. In the schedule
//    // above xo corresponds to GPU blocks.
//    f3.compute_at(f2, xo);
//
//    // Use the GPU threads for the x and y coordinates of the
//    // padded input.
//    f3.gpu_threads(x, y);
//
//    // JIT-compile the pipeline for the GPU. CUDA, OpenCL, or
//    // Metal are not enabled by default. We have to construct a
//    // Target object, enable one of them, and then pass that
//    // target object to compile_jit. Otherwise your CPU will very
//    // slowly pretend it's a GPU, and use one thread per output
//    // pixel.
//    printf("Target: %s\n", target.to_string().c_str());
//    f2.compile_jit(target);
//
//    return true;
//}

// This functions finds the determinant of Matrix

Expr determinantOfMatrix(Expr mat[3][3])
{
    Expr ans;
    ans = mat[0][0] * (mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2])
          - mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0])
          + mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
    return ans;
}

/* A utility function to calculate area of triangle formed by (x1, y1),
   (x2, y2) and (x3, y3) */
Expr area(Expr x1, Expr y1, Expr x2, Expr y2, Expr x3, Expr y3) {
    return abs((x1*(y2-y3) + x2*(y3-y1)+ x3*(y1-y2))/2.f);
}

/* A function to check whether point P(x, y) lies inside the triangle formed
   by A(x1, y1), B(x2, y2) and C(x3, y3) */
Expr is_inside(vector<Expr> point, vector<Expr> triangle) {
    Expr x1 = triangle[0];
    Expr y1 = triangle[1];
    Expr x2 = triangle[2];
    Expr y2 = triangle[3];
    Expr x3 = triangle[4];
    Expr y3 = triangle[5];
    Expr x = point[0];
    Expr y = point[1];
    /* Calculate area of triangle ABC */
    Expr A = area (x1, y1, x2, y2, x3, y3);

    /* Calculate area of triangle PBC */
    Expr A1 = area (x, y, x2, y2, x3, y3);

    /* Calculate area of triangle PAC */
    Expr A2 = area (x1, y1, x, y, x3, y3);

    /* Calculate area of triangle PAB */
    Expr A3 = area (x1, y1, x2, y2, x, y);

    Expr inside = (A == A1 + A2 + A3);

    /* Check if sum of A1, A2 and A3 is same as A */
    return inside;
}




//
//vector<Expr> convert_to_vector(ImageParam points) {
//    vector<Expr> p;
//    for(int i = 0; i < points.get().dim(0).extent(); i++) {
//        p.push_back(points(i));
//    }
//    return p;
//}



// g++ energy*generate.cpp -g -std=c++11 -I ../include -L ../bin -lHalide -lpthread -ldl -o energy_generate && ./energy_generate

int main(int argc, char **argv) {
    const int terms = 8;
    Var x("x"), y("y"), i("i");

    Param<int> width;
    Param<int> height;
    Param<double> mu;
    Param<double> beta;
    Param<double> sigma;
    Param<double> gamma;
    Param<double> tau;
    Param<int> order, samples;
    Param<double> learning_rate;
    learning_rate.set(0.0001);
    RDom s(0, width, 0, height);

    s.where(abs(sqrt(s.x * s.x + s.y * s.y)  - mu) < 3 * sigma);

//    int n = 5; //number of features
    ImageParam disparities(Float(64), 2);
    ImageParam left_features(Float(64), 3);
    ImageParam right_features(Float(64), 3);
    ImageParam a(Float(64), 1);
    ImageParam b(Float(64), 1);
    ImageParam c(Float(64), 1);
    ImageParam point_tri_map(Int(32), 2);

    Expr map_lo_bnd = -10;
    Expr map_up_bnd = 10;
    Func mean;
    mean(x, y) = a(clamp(point_tri_map(x, y), map_lo_bnd, map_up_bnd)) * x + b(clamp(point_tri_map(x, y), map_lo_bnd, map_up_bnd)) * y + c(clamp(point_tri_map(x, y), map_lo_bnd, map_up_bnd));

    Expr diff_lo_bnd = -10;
    Expr diff_up_bnd = 10;
    Func diff;
    diff(x, y, i) = clamp(left_features(x, y, i) - right_features(x, y, i), diff_lo_bnd, diff_up_bnd); // for i = 0 to n

    Expr norm_lo_bnd = -10;
    Expr norm_up_bnd = 10;
    Func norm;
    RDom n_feat(0, 5);
    norm(x, y) = clamp(sum(abs(diff(x, y, n_feat))), norm_lo_bnd, norm_up_bnd);

    ImageParam d(Int(32), 2);
    int w = 10, h = 10;
    Func test;
    test(x, y) = 0;

    Expr mean_lo_bnd = -10;
    Expr mean_up_bnd = 10;
    Func energy;
    energy(x, y) = beta * norm(x - clamp(d(x, y), -w/2, w/2), y) - log (gamma + exp(- (pow((clamp(d(x, y), -w/2, w/2) - mean(x, y)), 2)) / (2 * sigma * sigma)));
//
//    energy(x, y) = beta * norm(x - d(x, y) , y);// - log (gamma + exp(-(pow((d(x, y) - mean(x, y)), 2)) / (2 * sigma * sigma)));
    RDom img_x(0, width), img_y(0, height);
    Func err;
//    err() = f64(0);
    Expr err_lo_bnd = -10;
    Expr err_up_bnd = 10;
    err() = clamp(sum(img_y, sum(img_x, energy(img_x, img_y))), err_lo_bnd, err_up_bnd);

    auto d_err_d = propagate_adjoints(err);

    Func new_d;
    new_d(x, y) = d(x, y) - learning_rate * d_err_d(d)(x, y);

    Func new_disparities;
    new_disparities(x, y) = f32(0);

    Pipeline p({err, new_d});
    p.compile_to_static_library("energy", {disparities, a, b, c, learning_rate, beta, gamma , width, height, sigma, point_tri_map, left_features, right_features, d}, "energy");
}