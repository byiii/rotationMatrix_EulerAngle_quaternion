#ifndef HELPERS_H
#define HELPERS_H

#include <cfloat>
#include <cmath>

/// print a matrix
void printMatrix(const double *mat, unsigned rows, unsigned cols, bool rowMajor_p = true);


/// print a vector
void printVector(const double *vec, unsigned length, bool column_p=true);


/// get the sign of a double
inline double getSign(double x)
{
    if(fabs(x)>DBL_EPSILON)
    {
        if(x>0.0)
            return 1.0;
        else
            return -1.0;
    }
    else
        return 1.0;
}

#endif // HELPERS_H
