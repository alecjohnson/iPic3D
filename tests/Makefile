
INCLUDE = -I../include

OBJECTS = \
  test_arrays.o \
  ../utility/asserts.o \
  debug.o

# need -lrt to handle clock_gettime()
ARRAY_FLAGS = #-DCHAINED_ARRAYS #-DFLAT_ARRAYS #-DCHECK_BOUNDS
FLAGS = -O3 -DNO_MPI -fno-exceptions $(ARRAY_FLAGS) #-ggdb 

COMPILER = icpc #icpc # g++ # c++

testvbasic.mic: testvbasic.cpp
	icpc -mmic testvbasic.cpp -o testvbasic.mic

run.testvbasic.mic: testvbasic.mic
	ssh `hostname`-mic1 tests/testvbasic.mic

testvrinke.mic: testvrinke.cpp
	KMP_AFFINITY=balanced,granularity=thread icpc -mmic $(FLAGS) -fopenmp testvrinke.cpp $(INCLUDE) -o testvrinke.mic

run.testvrinke.mic: testvrinke.mic
	OMP_NUM_THREADS=50 mpirun -n 1 -host `hostname`-mic1 -env OMP_NUM_THREADS=50 ./testvrinke.mic

testv2: testv2.cpp Makefile
	$(COMPILER) $(FLAGS) -fopenmp testv2.cpp $(INCLUDE) -o testv2

#icpc -mmic -S -O3 testv2.cpp -o testv2.mic.s
testv2.mic: testv2.cpp Makefile
	KMP_AFFINITY=compact icpc -mmic $(FLAGS) -fopenmp testv2.cpp $(INCLUDE) -o testv2.mic

#icc -mmic -S -O3 testv2.cpp -std=c99
testv2.xeon: testv2.cpp Makefile
	KMP_AFFINITY=compact icpc $(FLAGS) -fopenmp testv2.cpp $(INCLUDE) -o testv2.xeon

run.testv2.mic: testv2.mic
	mpirun -n 1 -host `hostname`-mic1 -env OMP_NUM_THREADS=2 ./testv2.mic

run.testv2.xeon: testv2.xeon
	mpirun -n 1 -env OMP_NUM_THREADS=2 ./testv2.xeon

testv.mic: testv.cpp Makefile
	icpc testv.cpp -g -O3 -fopenmp -fp-model fast=2 -vec-report6 -mmic -o testv.mic

testv.xeon: testv.cpp Makefile
	icpc testv.cpp -g -O3 -fopenmp -vec-report6 -o testv.xeon

vtune.mic: testv.mic Makefile
	amplxe-cl -collect knc-general-exploration -knob enable-vpu-metrics=true -knob enable-tlb-metrics=true -- mpirun -host knc1-mic0 -np 1 -env OMP_NUM_THREADS 1 ./testv.mic > vtune.txt

test: clean_test_arrays test_arrays

test_arrays: clean_test_arrays $(OBJECTS)
	$(COMPILER) $(FLAGS) $(INCLUDE) $(OBJECTS) -o test_arrays

test_arrays.o:
	$(COMPILER) $(FLAGS) -c test_arrays.cpp $(INCLUDE) -o test_arrays.o


../utility/asserts.o:
	$(COMPILER) $(FLAGS) -c ../utility/asserts.cpp $(INCLUDE) -o ../utility/asserts.o

debug.o:
	$(COMPILER) $(FLAGS) -c ../utility/debug.cpp $(INCLUDE) -o debug.o

clean: clean_old_test_arrays clean_test_arrays clean_executables
	rm -f test_arrays $(OBJECTS)

clean_executables:
	rm -f *.mic *.xeon

clean_test_arrays:
	rm -f test_arrays test_arrays.o

old_test_arrays: clean_old_test_arrays
	$(COMPILER) $(FLAGS) old_test_arrays.cpp  -o old_test_arrays

clean_old_test_arrays:
	rm -f old_test_arrays old_test_arrays.o
