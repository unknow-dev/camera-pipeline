//
// Created by Prikshet on 8/4/2019.
//

#include<stdio.h>
#include<string>
#include "Halide.h"
#include <math.h>
#include "clock.h"
#include <vector>

#define R 0
#define G 1
#define B 2
using namespace std;
using namespace Halide::ConciseCasts;
using namespace Halide;

int main(int argc, char **argv) {
    Var x("x"), y("y");

    ImageParam left_buffer(type_of<uint8_t>(), 2);
    Expr x_clamped = clamp(x, 1, left_buffer.width() - 1);
    Expr y_clamped = clamp(y, 1, left_buffer.height() - 1);

    Func left_original("left_original");
    left_original(x, y) = f32(left_buffer(x_clamped, y_clamped));

    left_original.compute_root();

    // 3x3 sobel filter
    Buffer <float_t> sobel_1(3);
    sobel_1(0) = -1;
    sobel_1(1) = 0;
    sobel_1(2) = 1;

    Buffer <float_t> sobel_2(3);
    sobel_2(0) = 1;
    sobel_2(1) = 2;
    sobel_2(2) = 1;

    Func output_x_inter("output_x_inter");
    output_x_inter(x, y) = f32(left_original(x + 1, y) * sobel_1(0) + left_original(x, y) * sobel_1(1) + \
                               left_original(x - 1, y) * sobel_1(2));
    output_x_inter.compute_root();

    Func sobel_x("sobel_x");
    sobel_x(x, y) = f32(output_x_inter(x, y + 1) * sobel_2(0) + output_x_inter(x, y) * sobel_2(1) + \
                        output_x_inter(x, y - 1) * sobel_2(2));


    Func sobel_y_inter("sobel_y_inter");
    sobel_y_inter(x, y) = f32(left_original(x + 1, y) * sobel_2(0) + left_original(x, y) * sobel_2(1) + \
                              left_original(x - 1, y) * sobel_2(2));

    Func sobel_y("sobel_y");
    sobel_y(x, y) = f32(sobel_y_inter(x, y + 1) * sobel_1(0) + sobel_y_inter(x, y) * sobel_1(1) + \
                        sobel_y_inter(x, y - 1) * sobel_1(2));
    sobel_y.compute_root();

    Func sobel("sobel");
    sobel(x, y) = f32(sqrt(pow(sobel_x(x, y), 2) + pow(sobel_y(x, y), 2)));
    sobel.compute_root();

    RDom img(0, left_buffer.width(), 0, left_buffer.height());

//    Func sobel("sobel");
//    sobel(x, y) = f32(sqrt(sobel_x(x, y) * sobel_x(x, y) + sobel_y(x, y) * sobel_y(x, y)));

    Func max_x("max");
    max_x() = maximum(sobel_x(img.x, img.y));
    max_x.compute_root();

    Func min_x("min");
    min_x() = minimum(sobel_x(img.x, img.y));
    min_x.compute_root();

    Func max_y("max");
    max_y() = maximum(sobel_y(img.x, img.y));
    max_y.compute_root();

    Func min_y("min");
    min_y() = minimum(sobel_y(img.x, img.y));
    min_y.compute_root();

    Func max("max");
    max() = maximum(sobel(img.x, img.y));
    max.compute_root();

    Func min("min");
    min() = minimum(sobel(img.x, img.y));
    min.compute_root();

    Func sobel_x_norm("sobel_x_norm");
    sobel_x_norm(x, y) = f32(sobel_x(x, y) * 255 / (max_x() - min_x()));

    Func sobel_y_norm("sobel_y_norm");
    sobel_y_norm(x, y) = f32(sobel_y(x, y) * 255 / (max_y() - min_y()));

    Func sobel_norm("sobel_norm");
    sobel_norm(x, y) = f32(sobel(x, y) * 255 / (max() - min()));


    Func sobel_x_out("sobel_x_out");
    sobel_x_out(x, y) = f32(select(sobel_x_norm(x, y) <= 70, 0, sobel_x_norm(x, y)));
    sobel_x_out.compute_root();

    Func sobel_y_out("sobel_y_out");
    sobel_y_out(x, y) = f32(select(sobel_y_norm(x, y) <= 70, 0, sobel_y_norm(x, y)));
    sobel_y_out.compute_root();

    Func sobel_out("sobel_out");
    sobel_out(x, y) = f32(select(sobel_norm(x, y) <= 70, 0, sobel_norm(x, y)));
    sobel_out.compute_root();


    sobel_x_out.compile_to_static_library("sobel_x_out", {left_buffer}, "sobel_x_out");

    sobel_y_out.compile_to_static_library("sobel_y_out", {left_buffer}, "sobel_y_out");

    sobel_out.compile_to_static_library("sobel_out", {left_buffer}, "sobel_out");



    Target target;
    // And finally an iOS mach-o object file for one of Apple's 32-bit
    // ARM processors - the A6. It's used in the iPhone 5. The A6 uses
    // a slightly modified ARM architecture called ARMv7s. We specify
    // this using the target features field.  Support for Apple's
    // 64-bit ARM processors is very new in llvm, and still somewhat
    // flaky.
    target.os = Target::IOS;
    target.arch = Target::ARM;
    target.bits = 32;
    std::vector<Target::Feature> armv7s_features;
    armv7s_features.push_back(Target::ARMv7s);
    target.set_features(armv7s_features);
    sobel_x_out.compile_to_file("sobel_x_host_ios", args, "sobel_x_host_ios", target);
    sobel_y_out.compile_to_file("sobel_y_host_ios", args, "sobel_y_host_ios", target);
    sobel_out.compile_to_file("sobel_host_ios", args, "sobel_host_ios", target);

    // 32-bit arm iOS mach-o files start with the following magic bytes:
    uint32_t arm_32_ios_magic[] = {0xfeedface, // Mach-o magic bytes
                                   12,  // CPU type is ARM
                                   11,  // CPU subtype is ARMv7s
                                   1};  // It's a relocatable object file.
    f = fopen("sobel_x_host_ios.o", "rb");
    if (!f || fread(header, 32, 1, f) != 1) {
        printf("Object file not generated\n");
        return -1;
    }
    fclose(f);


    if (memcmp(header, arm_32_ios_magic, sizeof(arm_32_ios_magic))) {
        printf("Unexpected header bytes in 32-bit arm ios object file.\n");
        return -1;
    }

    f = fopen("sobel_y_host_ios.o", "rb");
    if (!f || fread(header, 32, 1, f) != 1) {
        printf("Object file not generated\n");
        return -1;
    }
    fclose(f);


    if (memcmp(header, arm_32_ios_magic, sizeof(arm_32_ios_magic))) {
        printf("Unexpected header bytes in 32-bit arm ios object file.\n");
        return -1;
    }

    f = fopen("sobel_host_ios.o", "rb");
    if (!f || fread(header, 32, 1, f) != 1) {
        printf("Object file not generated\n");
        return -1;
    }
    fclose(f);

    if (memcmp(header, arm_32_ios_magic, sizeof(arm_32_ios_magic))) {
        printf("Unexpected header bytes in 32-bit arm ios object file.\n");
        return -1;
    }

    // It looks like the object files we produced are plausible for
    // those targets. We'll count that as a success for the purposes
    // of this tutorial. For a real application you'd then need to
    // figure out how to integrate Halide into your cross-compilation
    // toolchain. There are several small examples of this in the
    // Halide repository under the apps folder. See HelloAndroid and
    // HelloiOS here:
    // https://github.com/halide/Halide/tree/master/apps/
    printf("Success!\n");
    return 0;

}
