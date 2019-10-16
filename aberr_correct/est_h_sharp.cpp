//
// Created by Prikshet on 8/16/2019.
//

// est_h_sharp is for optimizing the two parameters of the h_sharp_inv kernel.
#include<stdio.h>
#include<string>
#include "Halide.h"
#include "halide_image_io.h"
#include <math.h>
#include "clock.h"
#define R 0
#define G 1
#define B 2
using namespace std;
using namespace Halide::Tools;
using namespace Halide::ConciseCasts;
using namespace Halide;

int main(int argc, char **argv) {
    Halide::Buffer<uint8_t> in = load_image("../images/chroma_test/in256.jpg");
    Halide::Buffer<uint8_t> in_green(in.width(), in.height());

    for (int i = 0; i < in.width(); i++) {
        for(int j = 0; j < in.height(); j++) {
            in_green(i, j) = in(i, j, G);
        }
    }
    Expr xi = 9;
    Expr e = 2;

    Func in_b = BoundaryConditions::repeat_edge(in_green);

    ImageParam lambda(Float(64), 1);
    ImageParam sigma(Float(64), 1);

    Buffer<double> lambda_b(1);
    lambda_b(0) = 1;
    Buffer<double> sigma_b(1);
    sigma_b(0) = 1;
    lambda.set(lambda_b);
    sigma.set(sigma_b);



    RDom left(-(xi + e), -e);
    RDom right(e + 1, xi + e + 1);

    Var x, y;
    Func sum_left;
    sum_left(x, y) = sum(in_b(x + left, y));

    Func sum_right;
    sum_right(x, y) = sum(in_b(x + right, y));
    Func mean_left;
    mean_left(x, y) = sum_left(x, y) / 9.0f;
    Func mean_right;
    mean_right(x, y) = sum_right(x, y) / 9.0f;

    Func sigma_left;
    sigma_left(x, y) = sqrt(sum(pow(in_b(x + left, y) - mean_left(x, y), 2)) / 9.0f);

    Func sigma_right;
    sigma_right(x, y) = sqrt(sum(pow(in_b(x + right, y) - mean_right(x, y), 2)) / 9.0f);

    Func E_ring; //E_ring is the loss function that we are minimizing
    E_ring(x, y) = (sigma_left(x, y) + sigma_right(x, y)) / abs(mean_left(x, y) - mean_right(x, y));

    auto d_E_ring_d = Halide::propagate_adjoints(E_ring);
    Func d_E_ring_lambda = d_E_ring_d(lambda);
    Func d_E_ring_sigma = d_E_ring_d(sigma);

}