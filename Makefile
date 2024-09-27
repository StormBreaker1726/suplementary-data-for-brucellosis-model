CXX=gcc
INC=-I ./include
CFLAGS=-fopenmp $(INC)


exec: main.c ./source/*
	$(CXX) $(CFLAGS) -O3 -o $@ $^ -lm
debug: main.c ./source/*
	$(CXX) $(CFLAGS) -g -o $@ $^ -lm
all: exec debug
clean:
	rm -f exec debug
.PHONY: clean