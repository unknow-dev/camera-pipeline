//You asked it to compute the Func "support" over the same size as "left_buffer" when you call
//"support.realize(left_buffer.width(), left_buffer.height())", and it did that, hence the stores.
//Your boundary condition prevents it from loading out of bounds on the input by clamping the center coordinate of your window,
//which is why all of those values are three - they're all looking at the same window.
//If you wanted to realize support over the smaller inset window you can, check this tutorial for how:
//https://halide-lang.org/tutorials/tutorial_lesson_06_realizing_over_shifted_domains.html
//
//Getting the kth largest value in a language like Halide is not easy because there's no data-dependent control flow.
// One way to sort without branching is to use a sorting network (e.g. bitonic).
// You maintain a tuple containing the top k elements, and for each new batch of n items from the input
// you concatenate it to get a tuple (or probably a std::vector<Expr>) of size (k+n),
// then apply a sorting network so that the top k are in the first part of the vector,
// then you store just that portion back to the Func you're reducing into. If n == 1 you effectively
// get insertion sort and your sorting network can be quite simple (sketch of the code below).
// If n is proportionate to k you get better computational complexity, but the code is a little more annoying to write.
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
    int k = 4;
    Var x, y;
    Func source;
    source(x, y) = f32(x + y);

    Func argmin;
    std::vector<Expr> top_k(k * 3);
    for (int i = 0; i < k; i++) {
        top_k[3*i] = 10000.0f;
        top_k[3*i+1] = 0;
        top_k[3*i+2] = 0;
    }
    argmin() = Tuple(top_k);

    RDom r(0, 10, 0, 10);
    Expr next_val = source(r.x, r.y);
    Expr next_x = r.x;
    Expr next_y = r.y;
    top_k = Tuple(argmin()).as_vector();
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

    argmin() = Tuple(top_k);
    argmin.trace_stores();
    Realization real = argmin.realize();
}

