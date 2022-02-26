CXXFLAGS=-O3 -Wall # optimize and complain

FILENAME=nn
TEST=irisData
NODES=1

o: ${FILENAME}.cpp mat.o rand.o
	g++ -o nn ${FILENAME}.cpp mat.o rand.o

mat.o: mat.cpp rand.o
	g++ -c -o mat.o mat.cpp

rand.o: rand.cpp
	g++ -c -o rand.o rand.cpp


run: clean o
	./nn < ${TEST}.in > nn.out

debug: clean o
	./nn < ${TEST}.in

test: run
	diff nn.out ${TEST}${NODES}.out; [ $$? -eq 1 ]

clean:
	rm -f nn *.o nn.out
