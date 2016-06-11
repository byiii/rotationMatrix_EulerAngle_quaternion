#include "helpers.h"

#include <iostream>
#include <cmath>

void printMatrix(const double *mat, unsigned rows, unsigned cols, bool rowMajor_p)
{
    // detect the maximum number of digits will be printed, so cout.width
    // can be adapted to it.
    double max_value = -1E30;
    double min_value = 1e30;
    for(size_t i=0; i<rows*cols; ++i)
    {
        if(mat[i]>max_value)
            max_value = mat[i];
        if(mat[i]<min_value)
            min_value = mat[i];
    }
    max_value = fabs(max_value); min_value = fabs(min_value);
    max_value = std::fmax(max_value, min_value);
    // determine the cout.width number
    int cc = 12;
    while(max_value>1.0)
    {
        max_value /= 10.0;
        cc += 1;
    }

    // do the print job
    std::cout.precision(8); // print most 8 decimal places
    std::cout.setf( std::ios::fixed, std:: ios::floatfield); // always print 8 decimal places
    std::cout << "\n";

    if(rowMajor_p)
    {
        for(size_t r=0; r<rows; ++r)
        {
            std::cout << "|";
            unsigned idx = r*cols;
            for(size_t c=0; c<cols; ++c)
            {
                std::cout.width(cc); // set width
                std::cout << std::right << mat[idx+c];
            }
            std::cout << " |" << std::endl;
        }
        std::cout << std::endl;
    }
    else
    {
        for(size_t r=0; r<rows; ++r)
        {
            std::cout << "|";
            for(size_t c=0; c<cols; ++c)
            {
                std::cout.width(cc); // set width
                std::cout << std::right << mat[c*cols+r];
            }
            std::cout << " |" << std::endl;
        }
        std::cout << std::endl;
    }
}


/// print a vector
void printVector(const double *vec, unsigned length, bool column_p)
{
    // do the print job
    std::cout.precision(8); // print most 8 decimal places
    std::cout.setf( std::ios::fixed, std:: ios::floatfield); // always print 8 decimal places
    std::cout << "\n[";
    for(size_t i=0; i<length; ++i)
    {
        std::cout << " " << vec[i];
    }
    if(column_p)
        std::cout << "]^T" << std::endl;
    else
        std::cout << "]" << std::endl;
    std::cout << std::endl;
}
