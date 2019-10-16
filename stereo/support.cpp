//
// Created by Prikshet on 7/16/2019.
//
// TODO Delaunay triangulation
//#include "support.h"
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


    Var x("x"), y("y");

    Buffer<uint8_t> left_buffer = load_image("images/stereo/bike_smallest.jpg");

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

    Func sobel_x("sobel_x");
    sobel_x(x, y) = f32(output_x_inter(x, y + 1) * sobel_2(0) + output_x_inter(x,  y) * sobel_2(1) + output_x_inter(x, y - 1) * sobel_2(2));
//    sobel_x.compute_root(); didn't work
//    sobel_x.trace_stores();
    Buffer<float> output_x = sobel_x.realize(left_buffer.width(), left_buffer.height());
    output_x.set_name("output_x");

    Func sobel_y_inter("sobel_y_inter");
    sobel_y_inter(x, y) = f32(left_original(x + 1, y) * sobel_2(0) + left_original(x, y) * sobel_2(1) + left_original(x - 1, y) * sobel_2(2));

    Func sobel_y("sobel_y");
    sobel_y(x, y) = f32(sobel_y_inter(x, y + 1) * sobel_1(0) + sobel_y_inter(x, y) * sobel_1(1) + sobel_y_inter(x, y - 1) * sobel_1(2));
    sobel_y.compute_root();
//    sobel_y.trace_stores();
    Buffer<float> output_y = sobel_y.realize(left_buffer.width(), left_buffer.height());
    output_y.set_name("output_y");


    int k = 4; // # of support points
    vector<pair<Expr, Expr>> support_points(k * left_buffer.width() * left_buffer.height());
    // Calculate support pixel for each
    Func support("support");
    support(x, y) = Tuple(i32(0), i32(0), f32(0));


    vector <pair<Expr, Expr>> scan_range(2);
    pair <Expr, Expr> scan_height(0, (Expr) left_buffer.height());
    pair <Expr, Expr> scan_width;

    std::vector<Expr> top_k(k * 3);
    Expr predicate[3];
    int count = 0;
    RDom scanner;
    bool left, center, right;
    int which_pred = 0;
    Expr next_val, next_x, next_y, prev_val, prev_x, prev_y;
    Expr should_swap;
    pair<Expr, Expr> c;


    for (int yi = 0; yi < left_buffer.height(); yi++) {
        for (int xi = 0; xi < left_buffer.width() - 2; xi++) {
            double t1 = current_time();
            left = xi < left_buffer.width() / 4;
            center = (xi >= left_buffer.width() / 4 && xi < left_buffer.width() * 3 / 4);
            right = xi >= left_buffer.width() * 3 / 4;

            if (left) {

                    scan_width = make_pair((Expr) 0, (Expr) left_buffer.width() / 2);
                    which_pred = 0;
            }
            else if (center) {
                    scan_width = make_pair((Expr) xi - left_buffer.width() / 4, (Expr) left_buffer.width() / 2);
                    which_pred = 1;
            }
            else if (right) {
                    scan_width = make_pair((Expr) left_buffer.width() / 2, (Expr) left_buffer.width() / 2);
                    which_pred = 2;
            }
            else {
                cout<<"Error"<<endl;
            }

            scan_range = {scan_width, scan_height};

            // so far 1e-05 per iteration

            scanner = (scan_range); //slow 0.006

            // these three kind of slow 0.002

            predicate[0] = scanner.x != xi && scanner.y != yi;
            predicate[1] = scanner.x != 0 && scanner.y != 0;
            predicate[2] = scanner.x != xi && scanner.y != yi;
            scanner.where(predicate[which_pred]);


            // this loop is 1e-05
            for (int i = 0; i < k; i++) { // say we want top 4 support points.
                top_k[3*i] = 10000.0f;
                top_k[3*i+1] = 0;
                top_k[3*i+2] = 0;
            }

            Func argmin("argmin");

            // 0.003
            argmin() = Tuple(top_k);


            next_val = abs(output_x(xi, yi) - output_x(scanner.x, scanner.y)) + abs(output_y(xi, yi) - output_y(scanner.x, scanner.y));

            next_x = scanner.x;
            next_y = scanner.y;

            top_k = Tuple(argmin()).as_vector();

            // Insert a single element into a sorted list without actually branching
            top_k.push_back(next_val);
            top_k.push_back(next_x);
            top_k.push_back(next_y);

            for (int i = k; i > 0; i--) {
                prev_val = top_k[(i-1)*3];
                prev_x = top_k[(i-1)*3 + 1];
                prev_y = top_k[(i-1)*3 + 2];
                should_swap = top_k[i*3] < prev_val;

                top_k[(i-1)*3] = select(should_swap, top_k[i*3], prev_val);
                top_k[(i-1)*3 + 1] = select(should_swap, top_k[i*3 + 1], prev_x);
                top_k[(i-1)*3 + 2] = select(should_swap, top_k[i*3 + 2], prev_y);
                top_k[i*3] = select(should_swap, prev_val, top_k[i*3]);
                top_k[i*3 + 1] = select(should_swap, prev_x, top_k[i*3 + 1]);
                top_k[i*3 + 2] = select(should_swap, prev_y, top_k[i*3 + 2]);
            }

            // Discard the k+1th element
            top_k.pop_back(); top_k.pop_back(); top_k.pop_back();

            argmin() = Tuple(top_k);
            argmin.compute_root();
//            argmin.trace_stores();

            argmin.compile_to_lowered_stmt("argmin.html", {}, HTML);

            Realization real = argmin.realize(); //1.5 ms
            for (int i = 0; i < k; i++) {
                c.first = top_k[3*i+1];
                c.second = top_k[3*i+2];
                support_points.push_back(c);
            }

            count++;
            double t2 = current_time();
            cout<<count<<" "<<(t2-t1)/100<<" ms"<<endl;
        }
    }

    cout<<"executed"<<endl;
}