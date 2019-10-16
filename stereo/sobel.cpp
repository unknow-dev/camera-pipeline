//
// Created by Prikshet on 7/14/2019.
//

#include<stdio.h>
#include<string>
#include "Halide.h"
#include "halide_image_io.h"
#include <math.h>
#include "clock.h"
#include <vector>
#define R 0
#define G 1
#define B 2
using namespace std;
using namespace Halide::Tools;
using namespace Halide::ConciseCasts;
using namespace Halide;


int main(int argc, char **argv) {
    Var x, y, k, c;
    Buffer<uint8_t> left_buffer = load_image("images/stereo/bike.jpg");
    Expr clamped_x = clamp(x, 0, left_buffer.width() - 1);
    Expr clamped_y = clamp(y, 0, left_buffer.height() - 1);
    Func left_original("left_original");
    left_original(x, y) = f32(left_buffer(clamped_x, clamped_y));
    left_original.compute_root();

    // 3x3 sobel filter
    Buffer<float_t> sobel_1(3);
    sobel_1(0) = -1;
    sobel_1(1) = 0;
    sobel_1(2) = 1;

    Buffer<float_t> sobel_2(3);
    sobel_2(0) = 1;
    sobel_2(1) = 2;
    sobel_2(2) = 1;

    Func output_x_inter("output_x_inter");
    output_x_inter(x, y) = f32(left_original(x + 1, y) * sobel_1(0) + left_original(x, y) * sobel_1(1) + left_original(x - 1, y) * sobel_1(2));
    output_x_inter.compute_root();
//    output_x_inter.trace_stores();
    Func output_x("output_x");
    output_x(x, y) = f32(output_x_inter(x, y + 1) * sobel_2(0) + output_x_inter(x,  y) * sobel_2(1) + output_x_inter(x, y - 1) * sobel_2(2));
    output_x.compute_root();
//    output_x.trace_stores();
    RDom img(0, left_buffer.width(), 0, left_buffer.height());

    Func output_y_inter("output_y_inter");
    output_y_inter(x, y) = f32(left_original(x + 1, y) * sobel_2(0) + left_original(x, y) * sobel_2(1) + left_original(x - 1, y) * sobel_2(2));
    Func output_y("output_y");
    output_y(x, y) = f32(output_y_inter(x, y + 1) * sobel_1(0) + output_y_inter(x, y) * sobel_1(1) + output_y_inter(x, y - 1) * sobel_1(2));
    output_y.compute_root();

    Func output("output");
    output(x, y) = f32(sqrt(output_x(x, y) * output_x(x, y) + output_y(x, y) * output_y(x, y)));
    output.compute_root();
//    output.trace_stores();
    Func max("max");
    max() = f32(maximum(output(img.x, img.y)));
    max.compute_root();
//    max.trace_stores();
    Func min("min");
    min() = f32(minimum(output(img.x, img.y)));
    min.compute_root();
//    min.trace_stores();
    // output_inter for scaling
    Func output_inter("output_inter");
    output_inter(x, y) = f32((output(x, y) - min()) * 255 / (max() - min()));
    output_inter.compute_root();
//    output_u8_inter.trace_stores();

    Func output_u8("output_u8");
    output_u8(x, y) = u8(select(output_inter(x, y) <= 20, 0, output_inter(x, y)));
    output_u8.compute_root();
//    output_u8.trace_stores();

    Func gradient("gradient");
    gradient(x, y) = f32(atan(output_y(x, y) / output_x(x, y)));
//    gradient.trace_stores();


    Buffer<uint8_t> output_buff = output_u8.realize(left_buffer.width(), left_buffer.height());
    save_image(output_buff, "images/stereo/sobel/out.png");

//    Buffer<float_t> gradient_buff = gradient.realize(left_buffer.width(), left_buffer.height());

}