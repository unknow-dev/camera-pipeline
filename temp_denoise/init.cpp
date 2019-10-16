#include "Halide.h"
#include "halide_image_io.h"
#include <stdio.h>
#include "clock.h"
#include <vector>
using namespace Halide;
using namespace Halide::Tools;
using namespace Halide::ConciseCasts;

std::vector<std::vector<double>> random_offsets(int x, int y) {


}

std::vector<std::vector<double>> heapify(std::vector<std::vector<double>> offsets) {


}

std::vector<std::vector<double>>** init(int f) {
  Halide::Runtime::Buffer<uint8_t> frame[2H + 1];
  
  for (int i = -H; i < H; i++) {
    frame[i] = load_image("test1/small/" + std::to_string(f + i) + ".png");
  }
  
  int W = frame[0].width();
  int H = frame[0].height();

  std::vector<std::vector<double>> offsets[H][W];
  
  for (int i = 0; i < 2H + 1; i++) { // for all frames
    // vector of (x, y, d)'s
    offsets[y][x] = heapify(random_offsets(x, y));
  }
  
}
