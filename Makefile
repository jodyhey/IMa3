CC1 = mpicxx
CC2 = g++
FLAGS = -O3 -fpermissive


all:
	$(CC1) *.cpp -D MPI_ENABLED -D NDEBUG -D IMA3RELEASE -o IMa3 $(FLAGS)

singlecpu:
	$(CC2) *.cpp -U MPI_ENABLED -D NDEBUG -D IMA3RELEASE -o IMa3_singlecpu $(FLAGS)

testbed:
	$(CC1) *.cpp -D MPI_ENABLED -D NDEBUG -D STDTEST -o IMa3_stdtest $(FLAGS)

debug:
	$(CC1) *.cpp -D MPI_ENABLED -D DEBUG -o IMa3 $(FLAGS)

clean:
	rm IMa3 IMa3_singlecpu IMa3_stdtest
	rm *.o