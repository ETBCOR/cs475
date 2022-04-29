FLAGS=-O3 # -Wall optimize and complain

FILENAME=id7

INPUT=coke3.in

o: ${FILENAME}.cpp tree.o mat.o rand.o
	g++ ${FLAGS} -o ${FILENAME} ${FILENAME}.cpp mat.o rand.o tree.o

t: t.cpp mat.o rand.o
	g++ ${FLAGS} -o t t.cpp mat.o rand.o

tree.o: tree.cpp
	g++ ${FLAGS} -c -o tree.o tree.cpp

mat.o: mat.cpp rand.o
	g++ ${FLAGS} -c -o mat.o mat.cpp

rand.o: rand.cpp
	g++ ${FLAGS} -c -o rand.o rand.cpp

run: clean o
	./${FILENAME} < ${INPUT}

test: o
	rm -f z
	./localtest > z
	diff -sy z localtest.out

clean:
	rm -f ${FILENAME} *.o z
