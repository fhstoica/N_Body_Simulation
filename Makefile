GCC = g++
GLINK = g++

LIBPTHREAD = -lpthread
INCLUDEPTHREAD = -I/usr/include
LIBDIRPTHREAD = -L/usr/lib

CCFLAGS = -Wall -O3 -DLINUX 
GCCFLAGS = -Wall -O3 -DLINUX
DELAY = 20

all:star_cluster.srl

%.srl:%.cpp
	$(GCC) $(GCCFLAGS)  -o $@ $< $(INCLUDEPTHREAD) $(LIBPTHREAD)

clean:
	rm -f *.srl *~ *.ps *.mng *.o *.srk *.png *.txt

veryclean:
	rm -f *.srl *~ *.ps *.txt *.png *.o *.dat *.srk 
