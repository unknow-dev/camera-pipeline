optics: optics.cpp
	g++ optics.cpp -g -I ../include -I ../tools -L ../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -lpthread -ldl -o ./exec/optics -std=c++11 -Wall -Werror

supersample: supersample.cpp
	g++ supersample.cpp -g -I ../include -I ../tools -L ../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -lpthread -ldl -o ./exec/supersample -std=c++11 -Wall -Werror

chroma: chroma.cpp
	g++ chroma.cpp -g -I ../include -I ../tools -L ../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -lpthread -ldl -o ./exec/chroma -std=c++11 -Wall -Werror

chroma_: chroma_.cpp
	g++ chroma_.cpp -g -I ../include -I ../tools -L ../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -lpthread -ldl -o ./exec/chroma_ -std=c++11 -Wall -Werror

spatial_: spatial_denoise_.cpp
	g++ spatial_denoise_.cpp -g -I ../include -I ../tools -L ../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -lpthread -ldl -o ./exec/spatial_ -std=c++11 -Wall -Werror

white: white_balance.cpp
	g++ white_balance.cpp -g -I ../include -I ../tools -L ../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -o ./exec/white -std=c++11


gamma: gamma_correction.cpp
	g++ gamma_correction.cpp -g -I ../include -I ../tools -L ../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -o ./exec/gamma -std=c++11

test: test_reduction.cpp
	g++ test_reduction.cpp -g -I ../include -I ../tools -L ../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -lpthread -ldl -o ./exec/test -std=c++11 -Wall -Werror

fit: fit_function.cpp
	g++ fit_function.cpp -g -I ../include -I ../tools -L ../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -lpthread -ldl -o ./exec/fit -std=c++11 -Wall -Werror

est_h_sharp: est_h_sharp.cpp
	g++ est_h_sharp.cpp -g -I ../include -I ../tools -L ../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -lpthread -ldl -o ./exec/est_h_sharp -std=c++11 -Wall -Werror

diff_gauss: diff_gauss.cpp
	g++ diff_gauss.cpp -g -I ../include -I ../tools -L ../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -lpthread -ldl -o ./exec/diff_gauss -std=c++11 -Wall -Werror


diff_gauss_: diff_gauss_.cpp
	g++ diff_gauss_.cpp -g -I ../include -I ../tools -L ../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -lpthread -ldl -o ./exec/diff_gauss_ -std=c++11 -Wall -Werror

diff_gauss_bw: diff_gauss_bw.cpp
	g++ diff_gauss_bw.cpp -g -I ../include -I ../tools -L ../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -lpthread -ldl -o ./exec/diff_gauss_bw -std=c++11 -Wall -Werror

pre: pre_demosaic.cpp
	g++ pre_demosaic.cpp -g -I ../include -I ../tools -L ../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -lpthread -ldl -o ./exec/pre -std=c++11 -Wall -Werror
batch: batch.cpp
	g++ batch.cpp -g -I ../include -I ../tools -L ../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -lpthread -ldl -o ./exec/batch -std=c++11 -Wall -Werror
