//
// Created by Prikshet on 6/19/2019.
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

    Param<float> scan_range;
    Param<double> learning_rate;
    Var x, y, c;

    Buffer<uint8_t> left_buffer = load_image("images/stereo/left.png");
    Buffer<uint8_t> right_buffer = load_image("images/stereo/right.png");
    Buffer<uint8_t> ground_truth_buffer = load_image("images/stereo/ground_truth.png");

    RDom block(-9, 10, -9, 10);
    RDom scan(-25, 25);

    Func left_original = BoundaryConditions::repeat_edge(left_buffer);
    Func right_original = BoundaryConditions::repeat_edge(right_buffer);
    Func patch;
    patch(x, y, c) =

    auto d = propagate_adjoints(patch);
    Func left = d(left_original);
    Func right = d(right_original);


    left.trace_stores();
    Func color_average_left;
    RDom channel(0, 3);
    color_average_left(x, y) = sum(left(x, y, channel)) / 3;
    Func color_average_right;
    color_average_right(x, y) = sum(right(x, y, channel)) / 3;

    Func disparity_i16("disparity_i16");
    disparity_i16(x, y, c) = argmin(scan, sum(block, abs(i16(left(x + block.x, y + block.y , c)) - i16(right(x + block.x + scan, y + block.y, c)))))[0];
    disparity_i16(x, y, 3) = argmin(scan, sum(block, abs(i16(color_average_left(x + block.x, y + block.y)) - i16(color_average_right(x + block.x + scan, y + block.y)))))[0];
    disparity_i16.compute_root();
    disparity_i16.trace_stores();

    Func disparity;

    disparity(x, y, c) = u8(disparity_i16(x, y, c));
    disparity(x, y, 3) = u8(disparity_i16(x, y, 3));
    disparity.compute_root();
    Func err;
    err(x, y) = disparity(x, y, 3) - ground_truth_buffer(x, y);

    err.trace_stores();

    Buffer<uint8_t> err_buffer = err.realize(ground_truth_buffer.width(), ground_truth_buffer.height());
    Buffer<uint8_t> disparity_buffer = disparity.realize(ground_truth_buffer.width(), ground_truth_buffer.height());

    save_image(err_buffer, "images/stereo/err.png");
    save_image(disparity_buffer, "images/stereo/disparity.png");

    auto d_err_d = propagate_adjoints(err);

    //gradient descent

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
//        p.realize({e, c});
//
//        if (should_print) {
//            printf("Err: %g\n", e());
//        }
//
//        if (i == 0) {
//            initial_error = e();
//        }
//    }
//
//    double final_error = e();
//    if (final_error <= 1e-10 && final_error < initial_error) {
//        printf("[fit_function] Success!\n");
//        return 0;
//    } else {
//        printf("Did not converge\n");
//        return -1;
//    }
//}
//


    //disparity.trace_stores();

    int disparity_depth_mapping[116][3] ={{0, 0, 128}, {0, 0, 144}, {0, 0, 160}, {0, 0, 176}, {0, 0, 184}, {0, 0, 192}, {0, 0, 203}, {0, 0, 208}, {0, 0, 225}, {0, 0, 241},
                                   {0, 2, 255}, {0, 4, 255}, {0, 18, 255}, {0, 29, 255}, {0, 34, 255}, {0, 49, 255}, {0, 51, 255}, {0, 67, 255}, {0, 83, 255}, {0, 85, 255},
                                   {0, 99, 255}, {0, 111, 255}, {0, 115, 255}, {0, 131, 255}, {0, 132, 255}, {0, 148, 255}, {0, 164, 255}, {0, 167, 255}, {0, 180, 255}, {0, 196, 255},
                                   {0, 211, 255}, {0, 212, 255}, {0, 229, 255}, {0, 245, 255}, {1, 247, 254}, {6, 255, 249}, {12, 255, 243}, {22, 255, 233}, {31, 255, 224}, {38, 255, 217},
                                   {55, 255, 200}, {71, 255, 184}, {74, 255, 181}, {87, 255, 168}, {93, 255, 161}, {103, 255, 152}, {113, 255, 142}, {119, 255, 136}, {121, 255, 134}, {133, 255, 122},
                                   {136, 255, 119}, {152, 255, 103}, {156, 255, 99}, {168, 255, 87}, {175, 255, 80}, {184, 255, 71}, {194, 255, 61}, {200, 255, 55}, {202, 255, 53}, {214, 255,41},
                                   {217, 255, 38}, {218, 255, 37}, {233, 255, 22}, {237, 255, 18}, {249, 255, 6}, {255, 245, 0}, {255, 235, 0}, {255, 229, 0}, {255, 226, 0}, {255, 214,},
                                   {255, 212, 0}, {255, 206, 0}, {255, 196, 0}, {255, 192, 0}, {255, 180, 0}, {255, 164, 0}, {255, 153, 0}, {255, 148, 0}, {255, 145, 0}, {255, 134, 0},
                                   {255, 132, 0}, {255, 125, 0}, {255, 115, 0}, {255, 105, 0}, {255, 99, 0}, {255, 91, 0}, {255, 83, 0}, {255, 72, 0}, {255, 67, 0}, {255, 64,},
                                   {255, 53, 0}, {255, 51, 0}, {255, 44, 0}, {255, 34, 0}, {255, 24, 0}, {255, 18, 0}, {255, 10, 0}, {255, 2, 0}, {245, 1, 0}, {241, 0, 0},
                                   {237, 0, 0}, {225, 0, 0}, {218, 0, 0}, {208, 0, 0}, {198, 0, 0}, {192, 0, 0}, {179, 0, 0}, {176, 0, 0}, {175, 0, 0}, {164, 0, 0},
                                   {160, 0, 0}, {156, 0, 0}, {145, 0, 0}, {144, 0, 0}, {137, 0, 0}, {128, 0, 0}}; //116


//    Buffer<uint8_t> disparity_buffer = disparity.realize(left_buffer.width(), left_buffer.height());
//    save_image(disparity_buffer, "../images/stereo/disparity.png");
    Var d;
    Func disparity_depth_mapping_func;
    disparity_depth_mapping_func(d, c) = 0;

    for (int displacement = 0; displacement < 116; displacement++) {
        for (int color = 0; color <= 3; color++) { // last one is average
            disparity_depth_mapping_func(displacement, color) = disparity_depth_mapping[displacement][color];
        }
    }

    Var depth;
    Func depth_map;
    depth_map(x, y, c, depth) = u8(disparity_depth_mapping_func(u8(2.32f*disparity(x, y, c)), depth));


    // Saving to buffer
    Buffer<uint8_t> depth_map_buffer = depth_map.realize(left_buffer.width(), left_buffer.height(), 4, 3);

    string color_channel[4] = {"red", "blue", "green", "average"};

    for (int i = 0; i < 4; i++) {
        auto b = depth_map_buffer.sliced(2, i);
        save_image(b, "images/stereo/depth_map_"+color_channel[i]+"1.png");
    }

}
