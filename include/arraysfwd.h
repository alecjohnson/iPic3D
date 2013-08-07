/* forward declaration for array classes */
#ifndef arraysfwd_h
#define arraysfwd_h

template <class T>
class ConstArr3;
template <class T>
class ConstArr4;
template <class T>
class Arr1;
template <class T>
class Arr2;
template <class T>
class Arr3;
template <class T>
class Arr4;
template <class T>
class Array1;
template <class T>
class Array2;
template <class T>
class Array3;
template <class T>
class Array4;

// These aliases are defined for the following flexibilization purposes:
// - to avoid filling the code with template brackets
//   (i.e., to minimize explicitly template-dependent code).
// - so that they can be redefined according to the user's
//   preferred array implementation.
//
//typedef Arr1<int> intArr1;
//typedef Arr2<int> intArr2;
//typedef Arr3<int> intArr3;
//typedef Arr4<int> intArr4;
//typedef ConstArr1<double> doubleArr1;
//typedef ConstArr2<double> doubleArr2;
//
typedef ConstArr3<double> doubleCar3;
typedef ConstArr4<double> doubleCar4;
typedef Arr1<double> doubleArr1;
typedef Arr2<double> doubleArr2;
typedef Arr3<double> doubleArr3;
typedef Arr4<double> doubleArr4;
typedef Array1<double> doubleArray1;
typedef Array2<double> doubleArray2;
typedef Array3<double> doubleArray3;
typedef Array4<double> doubleArray4;
#endif
