# genvector
GenVector library source code
1) Clone the project first, please:
$git clone https://github.com/AndStorm/genvector.git
2) Then launch buildGenVector.sh file, please:
$./buildGenVector.sh
All the environment variables and compiler options are listed in this file.
It should build and install genvector library. It does if to use Intel compiler
(icc and icpc).
If to usePGI 19.4 (pgcc anf pgc++)
it gives a lot of such errors:
[  5%] Building CXX object CMakeFiles/genvector.dir/src/3DDistances.cxx.o
"/opt/pgi/linux86-64-llvm/19.4/include/edg/mmintrin.h", line 85: error:
          
argument of type "__v2si" is incompatible with parameter of type
 "__attribute((vector_size(8))) int"
return __builtin_ia32_vec_ext_v2si((__v2si)__m, 0);

^