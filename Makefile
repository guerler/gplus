GXXFLAGS = -O3 -ffast-math -fprefetch-loop-arrays -mfpmath=sse -funswitch-loops -DNDEBUG -Wno-deprecated

all:
	g++ -c -o gplus.o gplus.cpp $(GXXFLAGS)
	g++ -o gplus $(GXXFLAGS) *.o

clean:	
	rm -f *.o
