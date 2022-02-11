CXXFLAGS=-O3 -Wall # optimize and complain

FILENAME=nn
TEST=irisData

o: ${FILENAME}.cpp mat.o rand.o
	g++ -o a.o ${FILENAME}.cpp mat.o rand.o

mat.o: mat.cpp rand.o
	g++ -c -o mat.o mat.cpp

rand.o: rand.cpp
	g++ -c -o rand.o rand.cpp

run: o
	clear
	cat ${TEST}.in | ./a.o

test: o
	rm -f test.out
	cat ${TEST}.in | ./a.o >> test.out
	clear
	diff test.out ${TEST}.out

clean:
	rm *.o test.out
