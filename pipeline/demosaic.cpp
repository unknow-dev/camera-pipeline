/* 
   Author: Ziqi Zhang
   University of Rochester
   zzhang73 at u.rochester.edu
*/


#include "Halide.h"
#include <stdio.h>
#include "halide_image_io.h"
#include <float.h> 
#include <math.h>
#include <string>

#define R 0
#define G 1
#define B 2

using namespace Halide::Tools;
using namespace Halide;
using namespace std;
using namespace Halide::ConciseCasts;
/*
Compile with: g++ demosaic.cpp -g -I ../include -I ../tools -L ../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -o demosaic -std=c++11
DYLD_LIBRARY_PATH=../bin ./demosaic 
*/

/*
	The bayer filter pattern of test case from the reverse pipeline 
	G R G R G
	B G B G B
	G R G R G
	B G B G B

*/

int main(int argc, char **argv){
	// Buffer<uint8_t> input = Tools::load_image("images/rgb.jpg");
	Halide::Buffer<uint8_t> input = load_image("../images/demosaic_test/1out.jpg");
	// Halide::Buffer<uint8_t> output(input.width(), input.height());
	// input.set_min(1, 1);
	//First find the brightest pixel

	//float max = -1, mag, max_r, max_g, max_b;
	//	uint8_t offset_r, offset_g, offset_b;

	Halide::Func demosaic("demosaic");
	// demosaic.trace_stores();

	Var x("x"), y("y"), c("c");
	
	// Halide::Expr value;

	demosaic(x, y, c) = input(x, y, c);
	//Red

	RDom r(1, input.width() - 2, 1, input.height() - 2);

	
	r.where(r.x % 2 == 0);
	r.where(r.y % 2 == 0);

	demosaic(r.x, r.y, G) = (input(r.x, r.y + 1, G) / 4 + input(r.x, r.y - 1, G) / 4 + input(r.x + 1, r.y, G) / 4 + input(r.x - 1, r.y, G) / 4);
	demosaic(r.x, r.y, B) = (input(r.x + 1, r.y + 1, B) / 4 + input(r.x + 1, r.y - 1, B) / 4 + input(r.x - 1, r.y + 1, B) / 4 + input(r.x - 1, r.y - 1, B) / 4);


	//Blue

	RDom b(1, input.width() - 2, 1, input.height() - 2);

	b.where(b.x % 2 == 1);
	b.where(b.y % 2 == 1);

	demosaic(b.x, b.y, R) = (input(b.x + 1, b.y + 1, R) / 4 + input(b.x - 1, b.y + 1, R) / 4 + input(b.x + 1, b.y - 1, R) / 4 + input(b.x - 1, b.y - 1, R) / 4);
	demosaic(b.x, b.y, G) = (input(b.x + 1, b.y, G) / 4 + input(b.x, b.y + 1, G) / 4 + input(b.x - 1, b.y, G) / 4 + input(b.x, b.y - 1, G) / 4);

	//Green 01

	RDom g1(1, input.width() - 2, 1, input.height() - 2);
	g1.where((g1.x % 2 == 0) && (g1.y % 2 == 1));


	demosaic(g1.x, g1.y, R) = (input(g1.x, g1.y + 1, R) / 2 + input(g1.x, g1.y - 1, R) / 2);
	demosaic(g1.x, g1.y, B) = (input(g1.x + 1, g1.y, B) / 2 + input(g1.x - 1, g1.y, B) / 2);


	//Green 02

	RDom g2(1, input.width() - 2, 1, input.height() - 2);
	g2.where((g2.x % 2) == 1 && (g2.y % 2 == 0));

	demosaic(g2.x, g2.y, R) = (input(g2.x + 1, g2.y, R) / 2 + input(g2.x - 1, g2.y, R) / 2);
	demosaic(g2.x, g2.y, B) = (input(g2.x, g2.y + 1, B) / 2 + input(g2.x, g2.y - 1, B) / 2);
	demosaic.trace_stores();
	Halide::Buffer<uint8_t> output =
        demosaic.realize(input.width(), input.height(), input.channels());
	
	save_image(output, "balanced.png");


	return 0;
}
