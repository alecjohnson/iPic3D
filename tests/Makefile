
INCLUDE = -I../include

test_arrays_objects = \
  test_arrays.o \
  ../utility/asserts.o \
  debug.o

test_particles_objects = \
  debug.o

ARRAY_FLAGS = #-DCHAINED_ARRAYS #-DFLAT_ARRAYS #-DCHECK_BOUNDS
FLAGS = -O3 -DNO_MPI -fno-exceptions $(ARRAY_FLAGS) #-ggdb 

COMPILER = terminal-g++ #icpc # g++

test_particles: test_particles.cpp $(test_particles_objects)
	$(COMPILER) test_particles.cpp $(FLAGS) $(INCLUDE) $(test_particles_objects) -o test_particles

test: test.cpp ../include/Alloc.h
	$(COMPILER) test.cpp $(FLAGS) $(INCLUDE) -o test

test_arrays: clean_test_arrays $(test_arrays_objects)
	$(COMPILER) $(FLAGS) $(INCLUDE) $(test_arrays_objects) -o test_arrays

test_arrays.o:
	$(COMPILER) $(FLAGS) -c test_arrays.cpp $(INCLUDE) -o test_arrays.o


../utility/asserts.o:
	$(COMPILER) $(FLAGS) -c ../utility/asserts.cpp $(INCLUDE) -o ../utility/asserts.o

debug.o:
	$(COMPILER) $(FLAGS) -c ../utility/debug.cpp $(INCLUDE) -o debug.o

clean: clean_old_test_arrays clean_test_arrays
	rm -f test_arrays $(test_arrays_objects)

clean_test_arrays:
	rm -f test_arrays test_arrays.o

old_test_arrays: clean_old_test_arrays
	$(COMPILER) $(FLAGS) old_test_arrays.cpp  -o old_test_arrays

clean_old_test_arrays:
	rm -f old_test_arrays old_test_arrays.o
