
INCLUDE = -I../include

OBJECTS = \
  test_arrays.o \
  ../utility/asserts.o \
  debug.o

ARRAY_FLAGS = #-DCHAINED_ARRAYS #-DFLAT_ARRAYS #-DCHECK_BOUNDS
FLAGS = -O3 -DNO_MPI -fno-exceptions $(ARRAY_FLAGS) #-ggdb 

COMPILER = c++ #icpc # g++

testv.mic: testv.cpp Makefile
	icpc testv.cpp -g -O3 -fopenmp -vec-report6 -mmic -o testv.mic

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

clean: clean_old_test_arrays clean_test_arrays
	rm -f test_arrays $(OBJECTS)

clean_test_arrays:
	rm -f test_arrays test_arrays.o

old_test_arrays: clean_old_test_arrays
	$(COMPILER) $(FLAGS) old_test_arrays.cpp  -o old_test_arrays

clean_old_test_arrays:
	rm -f old_test_arrays old_test_arrays.o
