CXXFLAGS=-O3 -Wall # optimize and complain

# abbeyRoadAlbumCover.ppm
# bambooXuBeihong.pgm
# girlWithPearlEarring.ppm
# girlWithPearlEarringJohansson.ppm

FILENAME=pca
IN=girlWithPearlEarring.ppm
K=40

o: ${FILENAME}.cpp mat.o rand.o
	g++ -o ${FILENAME} ${FILENAME}.cpp mat.o rand.o

mat.o: mat.cpp rand.o
	g++ -c -o mat.o mat.cpp

rand.o: rand.cpp
	g++ -c -o rand.o rand.cpp

run: clean o
	./${FILENAME} ${K} < ${IN}

clean:
	rm -f ${FILENAME} *.o z.ppm z.pgm
