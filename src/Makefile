PROGRAM=../sphdec
CXX=g++
OBJS=main.o misc.o config.o algorithms.o sphdec.o
FLAGS=-std=c++11 -larmadillo -llapack -Wall -fopenmp -D_GLIBCXX_USE_CXX11_ABI=0
MODE=-O3
# debug: MODE=-g
# production: MODE=-O3

# make with=plotting
ifeq ($(findstring plotting, $(with)), plotting)
FLAGS+= -lboost_iostreams -lboost_system -lboost_filesystem -D PLOTTING
endif

# make with=lll
ifeq ($(findstring lll, $(with)), lll)
FLAGS+= -lgmp -lmpfr -I./fplll/include -L./fplll/lib -lfplll -D USE_LLL
endif

# make with=gpu
ifeq ($(findstring gpu, $(with)), gpu)
FLAGS+= -lnvblas
endif

# to compile with everything use something like
# make with=plotting+lll+gpu

# needs to come after -lnvblas flag
FLAGS+= -lblas


$(PROGRAM): $(OBJS)
	$(CXX) -o $(PROGRAM) $(OBJS) $(FLAGS)

main.o: main.cpp
	$(CXX) -c $(MODE) -o main.o main.cpp $(FLAGS)

config.o: config.cpp
	$(CXX) -c $(MODE) -o config.o config.cpp $(FLAGS)

misc.o: misc.cpp
	$(CXX) -c $(MODE) -o misc.o misc.cpp $(FLAGS)

algorithms.o: algorithms.cpp
	$(CXX) -c $(MODE) -o algorithms.o algorithms.cpp $(FLAGS)

sphdec.o: sphdec.cpp
	$(CXX) -c $(MODE) -o sphdec.o sphdec.cpp $(FLAGS)

.PHONY: clean
clean:
	-rm $(OBJS) $(PROGRAM)
