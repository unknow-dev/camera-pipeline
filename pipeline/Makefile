temporal: temporal.cpp heap.cpp aknn.cpp optical_flow.cpp
	g++ temporal.cpp heap.cpp aknn.cpp optical_flow.cpp -g -I ../include -I ../tools -L ../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -lpthread -ldl -o executables/temporal -std=c++11 `pkg-config opencv --cflags --libs` -Wall -Werror


ssd: ssd.cpp
	g++ ssd.cpp -g -I ../include -I ../tools -L ../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -lpthread -ldl -o executables/ssd -std=c++11 -Wall -Werror


pipeline: pipeline.cpp main.cpp
	g++ pipeline.cpp main.cpp -g -I ../include -I ../tools -L ../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -lpthread -ldl -o executables/pipeline -std=c++11 -Wall -Werror


spatial: spatial_denoise.cpp
	g++ spatial_denoise.cpp -g -I ../include -I ../tools -L ../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -lpthread -ldl -o exec/spatial -std=c++11 -Wall -Werror

spatial_: spatial_denoise_.cpp
	g++ spatial_denoise_.cpp -g -I ../include -I ../tools -L ../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -lpthread -ldl -o exec/spatial_ -std=c++11 -Wall -Werror

spatial_bw: spatial_denoise_bw.cpp
	g++ spatial_denoise_bw.cpp -g -I ../include -I ../tools -L ../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -lpthread -ldl -o exec/spatial_bw -std=c++11 -Wall -Werror

demosaic: demosaic.cpp
	g++ demosaic.cpp -g -I ../include -I ../tools -L ../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -lpthread -ldl -o exec/demosaic -std=c++11 -Wall -Werror


demosaic2: demosaic_2.cpp
	g++ demosaic_2.cpp -g -I ../include -I ../tools -L ../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -lpthread -ldl -o exec/demosaic2 -std=c++11 -Wall -Werror
