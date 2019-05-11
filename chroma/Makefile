optics: optics.cpp
	g++ optics.cpp -g -I ../include -I ../tools -L ../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -lpthread -ldl -o ./exec/optics -std=c++11 -Wall -Werror

supersample: supersample.cpp
	g++ supersample.cpp -g -I ../include -I ../tools -L ../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -lpthread -ldl -o ./exec/supersample -std=c++11 -Wall -Werror

chroma: chroma.cpp
	g++ chroma.cpp -g -I ../include -I ../tools -L ../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -lpthread -ldl -o ./exec/chroma -std=c++11 -Wall -Werror

est_h_sharp: est_h_sharp.cpp
	g++ est_h_sharp.cpp -g -I ../include -I ../tools -L ../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -lpthread -ldl -o ./exec/est_h_sharp -std=c++11 -Wall -Werror

diff_gauss: diff_gauss.cpp
	g++ diff_gauss.cpp -g -I ../include -I ../tools -L ../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -lpthread -ldl -o ./exec/diff_gauss -std=c++11 -Wall -Werror


diff_gauss_: diff_gauss_.cpp
	g++ diff_gauss_.cpp -g -I ../include -I ../tools -L ../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -lpthread -ldl -o ./exec/diff_gauss_ -std=c++11 -Wall -Werror

diff_gauss_bw: diff_gauss_bw.cpp
	g++ diff_gauss_bw.cpp -g -I ../include -I ../tools -L ../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -lpthread -ldl -o ./exec/diff_gauss_bw -std=c++11 -Wall -Werror
