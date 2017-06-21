PROGRAM=sphdec
OBJS=main.o misc.o config.o algorithms.o

$(PROGRAM): $(OBJS)
	g++ -o $(PROGRAM) $(OBJS) -std=c++14 -larmadillo -llapack -lblas -Wall -D_GLIBCXX_USE_CXX11_ABI=0

config.o: config.cpp
	g++ -c -O2 -g -o config.o -std=c++14 -larmadillo -llapack -lblas -Wall config.cpp -D_GLIBCXX_USE_CXX11_ABI=0

misc.o: misc.cpp
	g++ -c -O2 -g -o misc.o -std=c++14 -Wall misc.cpp -D_GLIBCXX_USE_CXX11_ABI=0

algorithms.o: algorithms.cpp
	g++ -c -O2 -g -o algorithms.o -std=c++14 -larmadillo -llapack -lblas -Wall algorithms.cpp -D_GLIBCXX_USE_CXX11_ABI=0

main.o: main.cpp
	g++ -c -O2 -g -o main.o -std=c++14 -larmadillo -llapack -lblas -Wall main.cpp -D_GLIBCXX_USE_CXX11_ABI=0

.PHONY: clean
clean:
	-rm $(OBJS) $(PROGRAM)