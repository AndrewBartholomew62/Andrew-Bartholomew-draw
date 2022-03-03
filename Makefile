
OBJFILES = util.o debug.o main.o drawfns.o class-control.o input.o \
           bigint-scalar.o triangulate.o circle_pack.o KS_circle_pack.o generic-code-util.o gauss-to-peer.o \
           force_direction.o magnify.o convex.o shrink_regions.o edge_distribution.o plot.o \
            gauss-orientation.o reidemeister.o laces.o
DEPS     = ./include/* 

COMPILE  = g++ -Wall -Wno-misleading-indentation -std=c++11 -I ./include -g -c $< -o $@ 

all: $(OBJFILES)
	g++ -o draw $(OBJFILES) 

util.o: ./src/util.cpp ./include/util.h $(DEPS)
	$(COMPILE)

debug.o: ./src/debug.cpp $(DEPS)
	$(COMPILE)

class-control.o: ./src/class-control.cpp $(DEPS)
	$(COMPILE)
	
bigint-scalar.o: ./src/bigint-scalar.cpp ./include/bigint-scalar.h	
	$(COMPILE)

input.o: ./src/input.cpp $(DEPS)
	$(COMPILE)

gauss-to-peer.o: ./src/gauss-to-peer.cpp $(DEPS)
	$(COMPILE)

generic-code-util.o: ./src/generic-code-util.cpp $(DEPS)
	$(COMPILE)

gauss-orientation.o: ./src/gauss-orientation.cpp $(DEPS)
	$(COMPILE)

reidemeister.o: ./src/reidemeister.cpp $(DEPS)
	$(COMPILE)

main.o: ./src/main.cpp $(DEPS)
	$(COMPILE)

drawfns.o: ./src/drawfns.cpp $(DEPS)
	$(COMPILE)

triangulate.o: ./src/triangulate.cpp $(DEPS)
	$(COMPILE)

#KS_circle_pack.o: ./src/KS_circle_pack.c
#	/usr/bin/gcc -g -O2 -s -std=gnu89 -fomit-frame-pointer -c $< -o $@

KS_circle_pack.o: ./src/KS_circle_pack.cpp
	$(COMPILE)

circle_pack.o: ./src/circle_pack.cpp $(DEPS)
	$(COMPILE)

force_direction.o: ./src/force_direction.cpp $(DEPS)
	$(COMPILE)

shrink_regions.o: ./src/shrink_regions.cpp $(DEPS)
	$(COMPILE)

edge_distribution.o: ./src/edge_distribution.cpp $(DEPS)
	$(COMPILE)

plot.o: ./src/plot.cpp $(DEPS)
	$(COMPILE)

magnify.o: ./src/magnify.cpp $(DEPS)
	$(COMPILE)

convex.o: ./src/convex.cpp $(DEPS)
	$(COMPILE)
	
laces.o: ./src/laces.cpp $(DEPS)
	$(COMPILE)
	

clean:
	rm -f $(OBJFILES)
