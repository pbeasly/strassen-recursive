# strassen-recursive
Measuring performance of Strassen recursive function for matrix multiplication.

Project measures the computation performance in MFLOPs for matrix multiplication with increasing matrix sizes.  The plot 'strassen_chart.png' shows how the performance efficiency of the Strassen recursive method improves with larger matrix sizes and is related to powers of 2 scaling.

File descriptions:
- strassen_chart.png:  Plot of computation performance in MFLOPs with increasing matrix size.
- strassen_mm.hpp:  C++ header file
- strassen_mm.cpp:  C++ source file
- strass_test.cpp:  C++ code for testing the strassen function
- xstrass.exe:      Executable

**Compilation Notes:**
Compiled using GNU C++ compiler.

g++ -std=c++14 -o xstrass strassen_mm.cpp strass_test.cpp

