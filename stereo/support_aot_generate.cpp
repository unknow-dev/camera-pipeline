#include<stdio.h>
#include<string>
#include "Halide.h"
#include <math.h>
#include "clock.h"
#include <vector>
// to compile: g++ support_aot_generate.cpp -g -std=c++11 -I ../include -L ../bin -lHalide -lpthread -ldl -o support_aot_generate && ./support_aot_generate
using namespace std;
using namespace Halide::ConciseCasts;
using namespace Halide;

int main(int argc, char **argv) {

    Var x("x"), y("y");


    Param<int> buffer_width;
    Param<int> buffer_height;

    ImageParam output_x(type_of<float>(), 2);
    ImageParam output_y(type_of<float>(), 2);


    std::vector<Argument> args(2);
    args[0] = output_x;
    args[1] = output_y;
    args[2] = buffer_width;
    args[3] = buffer_height;

    int k = 4; // # of support points

    vector <pair<Expr, Expr>> scan_range;
    pair <Expr, Expr> scan_height(0, (Expr) buffer_height);
    pair <Expr, Expr> scan_width((Expr) /*-buffer_width / 2*/ 0, (Expr) buffer_width / 2);

    scan_range = {scan_width, scan_height};
    RDom scanner(scan_range);
    Expr left_cond = (x < buffer_width / 4);
    Expr left_val = (scanner.x != x && scanner.y != y && scanner.x >= 0);
    Expr center_cond = (x >= buffer_width / 4 && x < buffer_width * 3 / 4);
    Expr center_val = (scanner.x != 0 && scanner.y != 0);
    Expr right_cond = (x >= buffer_width * 3 / 4);
    Expr right_val = (scanner.x != x && scanner.y != y && scanner.x < buffer_width / 2);

    scanner.where(select(left_cond, left_val, center_cond, center_val, select(right_cond, right_val, false)));

    std::vector<Expr> top_k(k * 3);
    for (int i = 0; i < k; i++) { // say we want top 4 support points.
        top_k[3*i] = 10000.0f;
        top_k[3*i+1] = 0.0f;
        top_k[3*i+2] = 0.0f;
    }

    Func argmin("argmin");
    argmin(x, y) = Tuple(top_k);
    Expr next_val = abs(output_x(x, y) - output_x(scanner.x, scanner.y)) + abs(output_y(x, y) - output_y(scanner.x, scanner.y));
    Expr next_x = f32(x + scanner.x);
    Expr next_y = f32(y + scanner.y);

    top_k = Tuple(argmin(x, y)).as_vector();
    // Insert a single element into a sorted list without actually branching
    top_k.push_back(next_val);
    top_k.push_back(next_x);
    top_k.push_back(next_y);

    for (int i = k; i > 0; i--) {
        Expr prev_val = top_k[(i-1)*3];
        Expr prev_x = top_k[(i-1)*3 + 1];
        Expr prev_y = top_k[(i-1)*3 + 2];
        Expr should_swap = top_k[i*3] < prev_val;

        top_k[(i-1)*3] = select(should_swap, top_k[i*3], prev_val);
        top_k[(i-1)*3 + 1] = select(should_swap, top_k[i*3 + 1], prev_x);
        top_k[(i-1)*3 + 2] = select(should_swap, top_k[i*3 + 2], prev_y);
        top_k[i*3] = select(should_swap, prev_val, top_k[i*3]);
        top_k[i*3 + 1] = select(should_swap, prev_x, top_k[i*3 + 1]);
        top_k[i*3 + 2] = select(should_swap, prev_y, top_k[i*3 + 2]);
    }

    // Discard the k+1th element
    top_k.pop_back(); top_k.pop_back(); top_k.pop_back();

    argmin(x, y) = Tuple(top_k);
//    argmin.trace_stores();

    argmin.compile_to_static_library("support", {output_x, output_y, buffer_width, buffer_height}, "argmin");
    argmin.compile_to_file("support_host", args, "support_host");

    // Target target;
    // // And finally an iOS mach-o object file for one of Apple's 32-bit
    // // ARM processors - the A6. It's used in the iPhone 5. The A6 uses
    // // a slightly modified ARM architecture called ARMv7s. We specify
    // // this using the target features field.  Support for Apple's
    // // 64-bit ARM processors is very new in llvm, and still somewhat
    // // flaky.
    // target.os = Target::IOS;
    // target.arch = Target::ARM;
    // target.bits = 32;
    // std::vector<Target::Feature> armv7s_features;
    // armv7s_features.push_back(Target::ARMv7s);
    // argmin.set_features(armv7s_features);
    // argmin.compile_to_file("support_host_ios", args, "support_host", target);

    // // 32-bit arm iOS mach-o files start with the following magic bytes:
    // uint32_t arm_32_ios_magic[] = {0xfeedface, // Mach-o magic bytes
    //                                12,  // CPU type is ARM
    //                                11,  // CPU subtype is ARMv7s
    //                                1};  // It's a relocatable object file.
    // f = fopen("support_host_ios.o", "rb");
    // if (!f || fread(header, 32, 1, f) != 1) {
    //     printf("Object file not generated\n");
    //     return -1;
    // }
    // fclose(f);

    // if (memcmp(header, arm_32_ios_magic, sizeof(arm_32_ios_magic))) {
    //     printf("Unexpected header bytes in 32-bit arm ios object file.\n");
    //     return -1;
    // }

    // // It looks like the object files we produced are plausible for
    // // those targets. We'll count that as a success for the purposes
    // // of this tutorial. For a real application you'd then need to
    // // figure out how to integrate Halide into your cross-compilation
    // // toolchain. There are several small examples of this in the
    // // Halide repository under the apps folder. See HelloAndroid and
    // // HelloiOS here:
    // // https://github.com/halide/Halide/tree/master/apps/
    // printf("Success!\n");
    // return 0;



}
