CXXFLAGS=-O3 -Wall # optimize and complain

FILENAME=nn
TEST=logic

o: ${FILENAME}.cpp mat.o rand.o
	g++ -o nn ${FILENAME}.cpp mat.o rand.o

mat.o: mat.cpp rand.o
	g++ -c -o mat.o mat.cpp

rand.o: rand.cpp
	g++ -c -o rand.o rand.cpp

run: o
	clear
	cat ${TEST}.in | ./nn

test: o
	rm -f test.out
	cat ${TEST}.in | ./a.o >> test.out
	clear
	diff test.out ${TEST}.out

clean:
	rm nn *.o test.out
