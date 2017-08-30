PROGRAM=../sphdec
OBJS=main.o misc.o config.o algorithms.o sphdec.o
CPPFLAGS=-I./fplll/include -L./fplll/lib$
FLAGS=-std=c++11 -larmadillo -llapack -lblas -Wall -fopenmp -D_GLIBCXX_USE_CXX11_ABI=0 -lboost_iostreams -lboost_system -lboost_filesystem -lgmp -lmpfr -I./fplll/include -L./fplll/lib -lfplll
MODE=-O3
# debug: MODE=-g
# production: MODE=-O3
# use GPU: -lnvblas

$(PROGRAM): $(OBJS)
	g++ -o $(PROGRAM) $(OBJS) $(FLAGS)

main.o: main.cpp
	g++ -c $(MODE) -o main.o main.cpp $(FLAGS)

config.o: config.cpp
	g++ -c $(MODE) -o config.o config.cpp $(FLAGS)

misc.o: misc.cpp
	g++ -c $(MODE) -o misc.o misc.cpp $(FLAGS)

algorithms.o: algorithms.cpp
	g++ -c $(MODE) -o algorithms.o algorithms.cpp $(FLAGS)

sphdec.o: sphdec.cpp
	g++ -c $(MODE) -o sphdec.o sphdec.cpp $(FLAGS)

.PHONY: clean
clean:
	-rm $(OBJS) $(PROGRAM)
