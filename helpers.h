#ifndef HELPERS_H
#define HELPERS_H

/// print a matrix
void printMatrix(const double *mat, unsigned rows, unsigned cols, bool rowMajor_p = true);


/// print a vector
void printVector(const double *vec, unsigned length, bool column_p=true);

#endif // HELPERS_H
