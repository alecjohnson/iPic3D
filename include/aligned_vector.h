
// this uses std::vector, which
// must include about 8000 lines
//
//#include "aligned_allocator.h"
//#include <vector> // needed for aligned_vector
//#define aligned_vector(type) std::vector<type, aligned_allocator<type, 64> >
//
// this approximate implementation of std::vector 
// includes only about 2800 lines
//
#include "Larray.h"
#define aligned_vector(type) Larray<type>
