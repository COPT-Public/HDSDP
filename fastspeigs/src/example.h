#ifndef example_h
#define example_h

#include "speigs.h"

spint n = 5;

/* Zero matrix */
spint Ap0[] = {0, 0, 0, 0, 0, 0};
spint Ai0[0];
double Ax0[0];

/* Diagonal matrix */
spint Ap1[6] = {0, 1, 2, 3, 3, 4};
spint Ai1[5] = {0, 1, 2, 4};
double Ax1[5] = {1.0, 2.0, 3.0, 4.0};

/* Two-two matrix */
spint Ap2[6] = {0, 1, 1, 1, 1, 1};
spint Ai2[1] = {3};
double Ax2[5] = {10.0};

/* Rank-one matrix */
spint Ap3[6] = {0, 0, 0, 2, 3, 3};
spint Ai3[3] = {2, 3, 3};
double Ax3[3] = {1.0, -2.0, 4.0};

/* Sparse submatrix */
spint Ap4[6] = {0, 0, 0, 2, 3, 3};
spint Ai4[3] = {2, 3, 3};
double Ax4[3] = {-1.0, -2.0, 100.0};

/* General matrix*/
spint Ap5[6] = {0, 5, 9, 12, 14, 15};
spint Ai5[15] = {0, 1, 2, 3, 4,
                    1, 2, 3, 4,
                       2, 3, 4,
                          3, 4,
                             4};
double Ax5[15] = {-1.0, -2.0, 1.0, 3.5, 912,
                   12.5,  -1.0 , -1.3, 50.0, 10.1,
                    13.1, 0.01, 0.5, 9.5, 20.3 };

/* Eigen value and vector */

double e[5] = {0.0};
double V[25] = {0.0};

#endif /* example_h */
