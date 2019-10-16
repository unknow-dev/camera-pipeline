//
// Created by Prikshet on 8/16/2019.
//


#include <stdio.h>
#include <string>
#include <iostream>
#include <random>
#include <math.h>
#include <algorithm>
#include <vector>
#include <list>
#include "Halide.h"
#include "halide_image_io.h"
#include "temporal.h"
#include <unordered_map>

using namespace std;
using namespace Halide;
using namespace Halide::ConciseCasts;


int main() {

    ImageParam frames(float(64), 4);

    Expr s = 10;
    Func D[n_frames], I[n_frames], Dx_right_in[n_frames], Dx_left_in[n_frames],
            Dy_up_in[n_frames], Dy_down_in[n_frames], Dx_right_out[n_frames], Dx_left_out[n_frames],
            Dy_up_out[n_frames], Dy_down_out[n_frames], weighted_ssd[n_frames][n_frames][K];

    Var x, y, x_i, y_i, c, i, xo, yo, xi, yi;

    RDom u(-s, s, -s, s);

    RDom dx_left_out(-s - 1, -s, -s, s);
    RDom dx_left_in(-s, -s + 1, -s, s);

    RDom dx_right_out(s, s + 1, -s, s);
    RDom dx_right_in(s - 1, s, -s, s);

    RDom dy_up_out(-s, s, -s - 1, -s);
    RDom dy_up_in(-s, s, -s, -s + 1);

    RDom dy_down_out(-s, s, s, s + 1);
    RDom dy_down_in(-s, s, s - 1, s);
    Func I = Tuple(0, 0, 0, 0, 0, 0, 0, 0, 0);

    I[frame](x, y, c) = frames(frame, clamp(x, s, width - s), clamp(y, s, height - s), c);

    D[frame](x, y, x_i, y_i, c) = i16(sum(pow((I[frame](x + u.x, y + u.y, c)
                                               - I[frame](x_i + u.x, y_i + u.y, c)), 2)));

    Dx_left_in[frame](x, y, x_i, y_i, c) = i16(sum(pow((I[frame](x + dx_left_in.x, y + dx_left_in.y, c)
                                                        - I[frame](x_i + dx_left_in.x, y_i + dx_left_in.y, c)), 2)));

    Dx_right_in[frame](x, y, x_i, y_i, c) = i16(sum(pow((I[frame](x + dx_right_in.x, y + dx_right_in.y, c)
                                                         - I[frame](x_i + dx_right_in.x, y_i + dx_right_in.y, c)), 2)));

    Dy_up_in[frame](x, y, x_i, y_i, c) = i16(sum(pow((I[frame](x + dy_up_in.x, y + dy_up_in.y, c)
                                                      - I[frame](x_i + dy_up_in.x, y_i + dy_up_in.y, c)), 2)));

    Dy_down_in[frame](x, y, x_i, y_i, c) = i16(sum(pow((I[frame](x + dy_down_in.x, y + dy_down_in.y, c)
                                                        - I[frame](x_i + dy_down_in.x, y_i + dy_down_in.y, c)), 2)));

    Dx_left_out[frame](x, y, x_i, y_i, c) = i16(sum(pow((I[frame](x + dx_left_out.x, y + dx_left_out.y, c)
                                                         - I[frame](x_i + dx_left_out.x, y_i + dx_left_out.y, c)), 2)));

    Dx_right_out[frame](x, y, x_i, y_i, c) = i16(sum(pow((I[frame](x + dx_right_out.x, y + dx_right_out.y, c)
                                                          - I[frame](x_i + dx_right_out.x, y_i + dx_right_out.y, c)), 2)));

    Dy_up_out[frame](x, y, x_i, y_i, c) = i16(sum(pow((I[frame](x + dy_up_out.x, y + dy_up_out.y, c)
                                                       - I[frame](x_i + dy_up_out.x, y_i + dy_up_out.y, c)), 2)));

    Dy_down_out[frame](x, y, x_i, y_i, c) = i16(sum(pow((I[frame](x + dy_down_out.x, y + dy_down_out.y, c)
                                                         - I[frame](x_i + dy_down_out.x, y_i + dy_down_out.y, c)), 2)));

    Dx_left_in.compile_to_static_library("x_l_in", {frames}, "x_l_in")
    Dx_right_in.compile_to_static_library("x_r_in",  {frames}, "x_r_in");
    Dy_up_in.compile_to_static_library("y_u_in", {frames}, "y_u_in");
    Dy_down_in.compile_to_static_library("y_d_in", {frames}, "y_d_in")
    Dx_left_out.compile_to_static_library("x_l_out", {frames}, "x_l_out");
    Dx_right_out.compile_to_static_library("x_r_out", {frames}, "x_r_out");
    Dy_up_out.compile_to_static_library("y_u_out", {frames}, "y_up_out");
    Dy_down_out.compile_to_static_library("y_d_out", {frames}, "y_d_out");
}
