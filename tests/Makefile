
INCLUDE = -I../include

OBJECTS = \
  test_arrays.o \
  ../utility/asserts.o \
  debug.o

FLAGS = -O3 -DNO_MPI #-DCHECK_BOUNDS
#FLAGS = -ggdb -O3 -DNO_MPI #-DCHECK_BOUNDS

test: clean_test_arrays test_arrays

test_arrays: clean_test_arrays $(OBJECTS)
	g++ $(FLAGS) $(INCLUDE) $(OBJECTS) -o test_arrays

test_arrays.o:
	g++ $(FLAGS) -c test_arrays.cpp $(INCLUDE) -o test_arrays.o


../utility/asserts.o:
	g++ $(FLAGS) -c ../utility/asserts.cpp $(INCLUDE) -o ../utility/asserts.o

debug.o:
	g++ $(FLAGS) -c ../utility/debug.cpp $(INCLUDE) -o debug.o
	
clean: clean_old_test_arrays clean_test_arrays
	rm -f test_arrays $(OBJECTS)

clean_test_arrays:
	rm -f test_arrays test_arrays.o

old_test_arrays: clean_old_test_arrays
	g++ $(FLAGS) old_test_arrays.cpp  -o old_test_arrays

clean_old_test_arrays:
	rm -f old_test_arrays old_test_arrays.o
