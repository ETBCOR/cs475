warmup.o: warmup.cpp mat.o rand.o
	g++ -o warmup warmup.cpp mat.o rand.o

mat.o: mat.cpp rand.o
	g++ -c -o mat.o mat.cpp

rand.o: rand.cpp
	g++ -c -o rand.o rand.cpp

run: warmup.o
	clear
	cat warmup1.in | ./warmup

test: warmup.o
	rm test.out
	cat warmup1.in | ./warmup >> test.out
	clear
	diff warmup.out test.out

clean:
	rm *.o test.out
