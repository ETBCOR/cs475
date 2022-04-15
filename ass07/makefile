FLAGS=-O3 # -Wall optimize and complain

FILENAME=kmeans

K=5
T=2
INPUT=mix5.in

o: ${FILENAME}.cpp mat.o rand.o
	g++ ${FLAGS} -o ${FILENAME} ${FILENAME}.cpp mat.o rand.o

#kdtree.o: mat.o
	#g++ ${FLAGS} -c -o kdtree.o kdtree.cpp

mat.o: mat.cpp rand.o
	g++ ${FLAGS} -c -o mat.o mat.cpp

rand.o: rand.cpp
	g++ ${FLAGS} -c -o rand.o rand.cpp

run: clean o
	./${FILENAME} ${K} ${T} < $(INPUT)

test: o
	rm -f my.out
	./localtest > my.out
	diff -y my.out localtest.out

clean:
	rm -f ${FILENAME} *.o my.out
