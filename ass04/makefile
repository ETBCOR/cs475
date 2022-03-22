 FLAGS=-O3 # -Wall optimize and complain

FILENAME=knn

K=1

# colorsUniqList.txt
# crayola16dupList.txt
# crayolaCrayonsList.txt
# smallList.txt

LIST=smallList.txt

# colorsUniqTest.txt
# crayolaCrayonsTest.txt
# randomColorsTest.txt
# smallTest.txt
# vandalTest.txt

TEST=smallTest.txt

o: ${FILENAME}.cpp mat.o rand.o
	g++ ${FLAGS} -o ${FILENAME} ${FILENAME}.cpp mat.o rand.o

mat.o: mat.cpp rand.o
	g++ ${FLAGS} -c -o mat.o mat.cpp

rand.o: rand.cpp
	g++ ${FLAGS} -c -o rand.o rand.cpp

run: clean o
	cat ${LIST} ${TEST} | ./${FILENAME} ${K}

clean:
	rm -f ${FILENAME}
