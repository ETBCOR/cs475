FLAGS=-O3 # -Wall optimize and complain

FILENAME=kdtree

LIST=smallList.txt

TEST=smallTest.txt

o: ${FILENAME}.cpp mat.o rand.o
	g++ ${FLAGS} -o ${FILENAME} ${FILENAME}.cpp mat.o rand.o

mat.o: mat.cpp rand.o
	g++ ${FLAGS} -c -o mat.o mat.cpp

rand.o: rand.cpp
	g++ ${FLAGS} -c -o rand.o rand.cpp

run: clean o
	cat ${LIST} ${TEST} | ./${FILENAME}

test: o
	rm -f my.out
	./localtest > my.out
	diff -y my.out localtest.out

clean:
	rm -f ${FILENAME} *.o my.out
