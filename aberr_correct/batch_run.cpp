//
// Created by Prikshet on 8/16/2019.
//



#include<stdio.h>
#include<string>
#include "Halide.h"
#include "halide_image_io.h"
#include <math.h>
#include "../../../Desktop/halide/tutorial/clock.h"
#include <vector>
#include "batch.h"

#include "red.h"
#include "blue.h"

#define R 0
#define G 1
#define B 2
using namespace std;
using namespace Halide::Tools;
using namespace Halide::ConciseCasts;
using namespace Halide;

int main(int argc, char **argv) {

    int batch_size = 2, width = 1920, height = 1200;
    Halide::Runtime::Buffer<uint8_t> buffer[batch_size];
    Halide::Runtime::Buffer<uint8_t> red_output[batch_size];
    Halide::Runtime::Buffer<uint8_t> blue_output[batch_size];
    Halide::Runtime::Buffer<uint8_t> green[batch_size];
    Halide::Runtime::Buffer<uint8_t> final[batch_size];
    for (int i = 0; i < batch_size; i++) {
        buffer[i] = load_image("/home/zendevil/Desktop/Halide/images/batch_test/" + to_string(i) + ".jpg");
    }

    int err_red = red({batch, red_output});
    if(err_red) {
        printf("error with red\n");
    }

    int err_blue = blue({batch, blue_output]});
    if(err_blue) {
        printf("error with blue\n");
    }

    int err_combine = combine({red_output, blue_output, green, final)})

    save_image(final, "corrected.png");

}
