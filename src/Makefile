PROGRAM=../sphdec
OBJS=main.o misc.o config.o algorithms.o sphdec.o
FLAGS=-std=c++14 -larmadillo -llapack -lblas -Wall -fopenmp -D_GLIBCXX_USE_CXX11_ABI=0

$(PROGRAM): $(OBJS)
	g++ -o $(PROGRAM) $(OBJS) $(FLAGS)

config.o: config.cpp
	g++ -c -O2 -g -o config.o config.cpp $(FLAGS)

misc.o: misc.cpp
	g++ -c -O2 -g -o misc.o -std=c++14 -Wall misc.cpp -D_GLIBCXX_USE_CXX11_ABI=0 -fopenmp

algorithms.o: algorithms.cpp
	g++ -c -O2 -g -o algorithms.o algorithms.cpp $(FLAGS)

sphdec.o: sphdec.cpp
	g++ -c -O2 -g -o sphdec.o sphdec.cpp $(FLAGS)

main.o: main.cpp
	g++ -c -O2 -g -o main.o main.cpp $(FLAGS)

.PHONY: clean
clean:
	-rm $(OBJS) $(PROGRAM)