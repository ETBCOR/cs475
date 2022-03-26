 FLAGS=-O3 # -Wall optimize and complain

FILENAME=kdtree

#colorsUniqList.txt
#crayola16dupList.txt
#crayolaCrayonsList.txt

LIST=crayola16upList.txt

TEST=colorsUniqTest.txt

o: ${FILENAME}.cpp mat.o rand.o
	g++ ${FLAGS} -o ${FILENAME} ${FILENAME}.cpp mat.o rand.o

mat.o: mat.cpp rand.o
	g++ ${FLAGS} -c -o mat.o mat.cpp

rand.o: rand.cpp
	g++ ${FLAGS} -c -o rand.o rand.cpp

run: clean o
	cat ${LIST} ${TEST} | ./${FILENAME}

clean:
	rm -f ${FILENAME}
