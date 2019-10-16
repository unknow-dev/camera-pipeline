//
// Created by Prikshet on 8/21/2019.
//

//
// Created by Prikshet on 8/14/2019.
//


// g++ energy*generate.cpp -g -std=c++11 -I ../include -L ../bin -lHalide -lpthread -ldl -o energy_generate && ./energy_generate
// g++ energy*gen.cpp -g -std=c++11 -I ../include -L ../bin -lHalide -lpthread -ldl -o energy_grad_des_gen && ./energy_grad_des_gen
#include <stdio.h>
#include <string>
#include "Halide.h"
#include <math.h>
#include "clock.h"
#include <vector>
#include "energy.h"


using namespace std;
using namespace Halide::ConciseCasts;
using namespace Halide;
//using namespace Halide::Tools;



// A helper function to check if OpenCL, Metal or D3D12 is present on the host machine.

#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif
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
//    Expr beta = 5.f, gamma = 10.f, tau = 10.f, learning_rate = 10.f ;
//
//    Param<int> order, samples;
//
//    Var x("x"), y("y");
//    ImageParam coeffs(Float(64), 1);
//
//    const int terms = 8;
//
//    ImageParam ground_truth(Float(64), 2);
//    ImageParam result(Float(64), 2);
//
//    Func err;
//
//    err(x, y) = pow((result(x, y) - ground_truth(x, y)) / ground_truth(x, y), 2);
//
//    RDom d_r(1, samples - 1);
//    Func average_err;
//    average_err() = sum(err(d_r)) / samples;
//
//    auto d_err_d = propagate_adjoints(average_err);
//
//    // Compute the new coefficients in terms of the old.
//
//    Func new_coeffs;
//    new_coeffs(x) = coeffs(x) - learning_rate * d_err_d(coeffs)(x);
//
//    Buffer<double> c_buffer(terms);
//    order.set(terms);
//    samples.set(1000);
//    auto e = Buffer<double>::make_scalar();
//    coeffs.set(c_buffer);
//
//    c_buffer.fill(0);
//    // Initialize to the Taylor series for sin about zero
//    c_buffer(0) = 1;
//    for (int i = 1; i < terms; i++) {
//        c_buffer(i) = -c_buffer(i-1)/(i*2*(i*2 + 1));
//    }
//    learning_rate.set(0.00001);
//    const int steps = 10000;
//    double initial_error = 0.0;
//    int err_aot;
//    Halide::Runtime::Buffer<double> disparities(100, 100);
//    for (int i = 0; i < 100; i++) {
//        for (int j = 0; j < 100; j++) {
//            disparities(i, j) = 0;
//        }
//    }
//    Halide::Runtime::Buffer<double> a(10);
//    Halide::Runtime::Buffer<double> b(10);
//    Halide::Runtime::Buffer<double> c(10);
//    Halide::Runtime::Buffer<double> point_tri_map(disparities.width(), disparities.height());
//    Halide::Runtime::Buffer<double> left_features(10);
//    Halide::Runtime::Buffer<double> right_features(10);
//    Halide::Runtime::Buffer<double> d(disparities.width(), disparities.height());
//
//
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
////        p.compile_to_static_library("energy", {disparities, a, b, c, beta, gamma , height, width,  \
////    sigma, point_tri_map, left_features, right_features, d}, "energy");
//        err_aot = p(disparities, a, b, c, learning_rate, disparities.width(), disparities.height(), point_tri_map, left_features, right_features, d);
////
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
//
//// This functions finds the determinant of Matrix
//
//Expr determinantOfMatrix(Expr mat[3][3])
//{
//    Expr ans;
//    ans = mat[0][0] * (mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2])
//          - mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0])
//          + mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
//    return ans;
//}
//
///* A utility function to calculate area of triangle formed by (x1, y1),
//   (x2, y2) and (x3, y3) */
//Expr area(Expr x1, Expr y1, Expr x2, Expr y2, Expr x3, Expr y3) {
//    return abs((x1*(y2-y3) + x2*(y3-y1)+ x3*(y1-y2))/2.f);
//}
//
///* A function to check whether point P(x, y) lies inside the triangle formed
//   by A(x1, y1), B(x2, y2) and C(x3, y3) */
//Expr is_inside(vector<Expr> point, vector<Expr> triangle) {
//    Expr x1 = triangle[0];
//    Expr y1 = triangle[1];
//    Expr x2 = triangle[2];
//    Expr y2 = triangle[3];
//    Expr x3 = triangle[4];
//    Expr y3 = triangle[5];
//    Expr x = point[0];
//    Expr y = point[1];
//    /* Calculate area of triangle ABC */
//    Expr A = area (x1, y1, x2, y2, x3, y3);
//
//    /* Calculate area of triangle PBC */
//    Expr A1 = area (x, y, x2, y2, x3, y3);
//
//    /* Calculate area of triangle PAC */
//    Expr A2 = area (x1, y1, x, y, x3, y3);
//
//    /* Calculate area of triangle PAB */
//    Expr A3 = area (x1, y1, x2, y2, x, y);
//
//    Expr inside = (A == A1 + A2 + A3);
//
//    /* Check if sum of A1, A2 and A3 is same as A */
//    return inside;
//}
//



//
//vector<Expr> convert_to_vector(ImageParam points) {
//    vector<Expr> p;
//    for(int i = 0; i < points.get().dim(0).extent(); i++) {
//        p.push_back(points(i));
//    }
//    return p;
//}

int main(int argc, char **argv) {
    int terms = 10; // what is terms?
    Buffer<double> co(terms);
    Param<double> order, samples, coeffs, learning_rate;
    order.set(terms);
    samples.set(1000);
    coeffs.set(co);

    co.fill(0);
    // Initialize to the Taylor series for sin about zero
    co(0) = 1;
    for (int i = 1; i < terms; i++) {
        co(i) = -co(i-1)/(i*2*(i*2 + 1));
    }

    // This gradient descent. We'll use a very slow learning rate and lots of
    // steps.

    samples.set(1000);
    const int steps = 10000;
    auto e = Buffer<double>::make_scalar();
    double initial_error = 0.0;
    float temp = 1000;
    learning_rate.set(0.0001);


    for (int i = 0; i <= steps; i++) {
//        learning_rate.set(exp(-temp)); //simulated annealing
        bool should_print = (i == 0 || i == steps/2 || i == steps);
        if (should_print) {
            printf("Iteration %d\n"
                   "Coefficients: ", i);
            for (int j = 0; j < terms; j++) {
                printf("%g ", co(j));
            }
            printf("\n");
        }



/* GPU scheduling */
//        Target target = find_gpu_target();
//        if (!target.has_gpu_feature()) {
//            return false;
//        }

// Let's compute the look-up-table using the GPU in 16-wide
// one-dimensional thread blocks. First we split the index
// into blocks of size 16:
        Var block, thread, j;
        new_disparities.split(j, block, thread, 16);
        err.split(j, block, thread, 16);
// Then we tell cuda that our Vars 'block' and 'thread'
// correspond to CUDA's notions of blocks and threads, or
// OpenCL's notions of thread groups and threads.
        new_disparities.gpu_blocks(block)
                .gpu_threads(thread);
        err.gpu_blocks(block)
                .gpu_threads(thread);

/* ------------ */



        if(err) {
            printf("Error\n");
        }

        if (should_print) {
            printf("Err: %g\n", e());
        }

        if (i == 0) {
            initial_error = e();
        }
    }

//    Target target;
//
//// And finally an iOS mach-o object file for one of Apple's 32-bit
//// ARM processors - the A6. It's used in the iPhone 5. The A6 uses
//// a slightly modified ARM architecture called ARMv7s. We specify
//// this using the target features field.  Support for Apple's
//// 64-bit ARM processors is very new in llvm, and still somewhat
//// flaky.
//    target.os = Target::IOS;
//    target.arch = Target::ARM;
//    target.bits = 32;
//    std::vector<Target::Feature> armv7s_features;
//    armv7s_features.push_back(Target::ARMv7s);
//    target.set_features(armv7s_features);
//    std::vector<Argument> args(5);
//    args[0] = disparities;
//    args[1] = left_features;
//    args[2] = right_features;
//    args[3] = support_points_l;
//    args[4] = support_points_r;
//    p.compile_to_file("support_host_ios", args, "support_host", target);
//
//// 32-bit arm iOS mach-o files start with the following magic bytes:
//    uint32_t arm_32_ios_magic[] = {0xfeedface, // Mach-o magic bytes
//                                   12,  // CPU type is ARM
//                                   11,  // CPU subtype is ARMv7s
//                                   1};  // It's a relocatable object file.
//    FILE *f = fopen("energy_generate_host_ios.o", "rb");
//    uint8_t header[32];
//    if (!f || fread(header, 32, 1, f) != 1) {
//        printf("Object file not generated\n");
//        return -1;
//    }
//    fclose(f);
//
//    if (memcmp(header, arm_32_ios_magic, sizeof(arm_32_ios_magic))) {
//        printf("Unexpected header bytes in 32-bit arm ios object file.\n");
//        return -1;
//    }
//
//// It looks like the object files we produced are plausible for
//// those targets. We'll count that as a success for the purposes
//// of this tutorial. For a real application you'd then need to
//// figure out how to integrate Halide into your cross-compilation
//// toolchain. There are several small examples of this in the
//// Halide repository under the apps folder. See HelloAndroid and
//// HelloiOS here:
//// https://github.com/halide/Halide/tree/master/apps/
//    printf("Success!\n");
//    return 0;
}


