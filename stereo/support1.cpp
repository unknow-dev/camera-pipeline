

//////////////
//
// Created by Prikshet on 7/16/2019.
//
// TODO Delaunay triangulation
#include "support.h"
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

    cout<<"executed"<<endl;
    Var x("x"), y("y");

    Buffer<uint8_t> left_buffer = load_image("images/stereo/bike_small.jpg");

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

    Func output_x("output_x");
    output_x(x, y) = f32(output_x_inter(x, y + 1) * sobel_2(0) + output_x_inter(x,  y) * sobel_2(1) + output_x_inter(x, y - 1) * sobel_2(2));
    output_x.compute_root();

    Func output_y_inter("output_y_inter");
    output_y_inter(x, y) = f32(left_original(x + 1, y) * sobel_2(0) + left_original(x, y) * sobel_2(1) + left_original(x - 1, y) * sobel_2(2));

    Func output_y("output_y");
    output_y(x, y) = f32(output_y_inter(x, y + 1) * sobel_1(0) + output_y_inter(x, y) * sobel_1(1) + output_y_inter(x, y - 1) * sobel_1(2));
    output_y.compute_root();
//    Var n("n");



// Calculate support pixel for each
    Func support("support");
    support(x, y) = Tuple(i32(0), i32(0), f32(0));

    RDom scan_left(0, left_buffer.width() / 4, 0, left_buffer.height());
    scan_left.where(scan_left.x != x && scan_left.y != y);
//    Expr first = argmin(abs(output_x(x, y) - output_x(scan_left.x, scan_left.y)) + abs(output_y(x, y) - output_y(scan_left.x, scan_left.y)));

    RDom scan_center(-left_buffer.width() / 4, left_buffer.width() / 2, 0, left_buffer.height());
    scan_center.where(scan_center.x != 0 && scan_center.y != 0);
//    Expr second = argmin(abs(output_x(x, y) - output_x(x + scan_center.x, scan_center.y)) + abs(output_y(x, scan_center.y) - output_y(x + scan_center.x, scan_center.y)));

    RDom scan_right(left_buffer.width() * 3/4, left_buffer.width() / 4, 0, left_buffer.height());
    scan_right.where(scan_right.x != x && scan_right.y != y);
//    Expr third = argmin(abs(output_x(x, y) - output_x(scan_right.x, scan_right.y)) + abs(output_y(x, y) - output_y(scan_right.x, scan_right.y)));

    // select based on x value
//    support(x, y) = select(x < width / 4, first, x < width * 3 / 4, second, third);

    int k = 4;
    Func argmin("argmin");

    std::vector<Expr> top_k[left_buffer.height()][left_buffer.width()];


    for (int xi = 0; xi < left_buffer.width(); xi++) {
        for (int yi = 0; yi < left_buffer.height(); yi++) {
            vector<Expr> temp(k * 3);
            for (int i = 0; i < k; i++) { // say we want top 4 support points.
                temp[3*i] = 10000.0f;
                temp[3*i+1] = 0;
                temp[3*i+2] = 0;
            }
            top_k[yi][xi] = temp;
        }
    }
    for (int i = 0; i < left_buffer.height(); i++) {
        for (int j = 0; j < left_buffer.width(); j++) {
            cout << "Val: "<<top_k[i][j][0]<<endl;
            cout << "x: "<<top_k[i][j][1]<<endl;
            cout << "y: "<<top_k[i][j][2]<<endl;

        }
    }


    argmin(x, y) = Tuple(top_k[0][0]);

    RDom H(0, left_buffer.height());
    for (int xi = 0; xi < left_buffer.width(); xi++) {
        for (int yi = 0; yi < left_buffer.height(); yi++) {
            cout<<xi<<" "<<yi<<endl;

            Expr next_val = abs(output_x(xi, yi) - output_x(xi + scan_center.x, scan_center.y)) + abs(output_y(xi, scan_center.y) - output_y(xi + scan_center.x, scan_center.y));

            Expr next_x = scan_center.x;
            Expr next_y = scan_center.y;

            top_k[yi][xi] = Tuple(argmin(xi, yi)).as_vector();
            // Insert a single element into a sorted list without actually branching
            top_k[yi][xi].push_back(next_val);
            top_k[yi][xi].push_back(next_x);
            top_k[yi][xi].push_back(next_y);
            for (int i = k; i > 0; i--) {
                cout<<"i= "<<i<<endl;
                Expr prev_val = top_k[yi][xi][(i-1)*3];
                Expr prev_x = top_k[yi][xi][(i-1)*3 + 1];
                Expr prev_y = top_k[yi][xi][(i-1)*3 + 2];
                Expr should_swap = top_k[yi][xi][i*3] < prev_val;

                top_k[yi][xi][(i-1)*3] = select(should_swap, top_k[yi][xi][i*3], prev_val);
                top_k[yi][xi][(i-1)*3 + 1] = select(should_swap, top_k[yi][xi][i*3 + 1], prev_x);
                top_k[yi][xi][(i-1)*3 + 2] = select(should_swap, top_k[yi][xi][i*3 + 2], prev_y);
                top_k[yi][xi][i*3] = select(should_swap, prev_val, top_k[yi][xi][i*3]);
                top_k[yi][xi][i*3 + 1] = select(should_swap, prev_x, top_k[yi][xi][i*3 + 1]);
                top_k[yi][xi][i*3 + 2] = select(should_swap, prev_y, top_k[yi][xi][i*3 + 2]);
            }
            // Discard the k+1th element
            top_k[yi][xi].pop_back(); top_k[yi][xi].pop_back(); top_k[yi][xi].pop_back();
        }
    }

    for (int xi = 0; xi < left_buffer.width(); xi++) {
        for (int yi = 0; yi < left_buffer.height(); yi++) {
            argmin(xi, yi) = Tuple(top_k[yi][xi]);
        }
    }

    argmin.compute_root();
    argmin.trace_stores();
    argmin.compile_to_lowered_stmt("argmin.html", {}, HTML);

    Realization real = argmin.realize(left_buffer.width(), left_buffer.height());

}
