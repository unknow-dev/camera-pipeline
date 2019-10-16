//
// Created by Prikshet on 8/4/2019.
//

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
#include "support.h"

#define R 0
#define G 1
#define B 2
using namespace std;
using namespace Halide::Tools;
using namespace Halide::ConciseCasts;
using namespace Halide;

int main(int argc, char **argv) {

    int k = 4;
    vector<pair<Expr, Expr>> support_points(k * left_buffer.width() * left_buffer.height());
    Buffer<int> coord_x[4];
    Buffer<int> coord_y[4];
    Buffer<float> val[4];
    // argmin.compile_to_static_library("support", {output_x, output_y, buffer_width, buffer_height}, "argmin");

    int error = argmin(output_x, output_y, left_buffer.width(), left_buffer.height(), \
    coord_x[0], coord_y[0], val[0],\
    coord_x[1], coord_y[1], val[1],\
    coord_x[2], coord_y[2], val[2],\
    coord_x[3], coord_y[3], val[3]);

    if (error) {
        cout << "Halide returned error:" << error << endl;
        return -1;
    }

    for (int yi = 0; yi < left_buffer.height(); yi++) {
        for (int xi = 0; xi < left_buffer.width() - 2; xi++) {
//            cout << argmin(xi, yi) << endl;
        }
    }
}