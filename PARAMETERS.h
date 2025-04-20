#pragma once

#define PI 3.1415926535897932384626433832795029

#define M 10000  // 10000

#define x1 0.
#define x2 10.//10.
#define L (x2 - x1)

#define t1 0.
#define t2 10.//10.

#define h 0.02 //0.02
#define t 0.01 //0.01
#define r (t/h)

#define N (L / h)
#define T ((t2 - t1) / t)

#define SIZE_X (int)(N)
#define SIZE_T (int)(T)

#define m 2.
#define g 100.

#define TMP 5.

#define al 15.
#define x_q 5.
#define sigma 0.085 // 10.*h/(2*sqrt(2*ln(2))) = 0.085 ; 10 = length in terms of coordinate steps
//#define sigma 0.17    // here 10 -> 20 - width