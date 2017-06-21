PROGRAM=sphdec
OBJS=main.o misc.o config.o

$(PROGRAM): $(OBJS)
	g++ -o $(PROGRAM) $(OBJS) -std=c++14 -larmadillo -llapack -lblas -Wall

config.o: config.cpp
	g++ -c -O2 -g -o config.o -std=c++14 -larmadillo -llapack -lblas -Wall config.cpp

misc.o: misc.cpp
	g++ -c -O2 -g -o misc.o -std=c++14 -Wall misc.cpp

main.o: main.cpp
	g++ -c -O2 -g -o main.o -std=c++14 -larmadillo -llapack -lblas -Wall main.cpp

.PHONY: clean
clean:
	-rm $(OBJS) $(PROGRAM)