
#include "aligned_allocator.h"
#include <vector> // needed for aligned_vector
#define aligned_vector(type) std::vector<type, aligned_allocator<type, 64> >
//#include "Larray.h"
//#define aligned_vector(type) Larray<type>
