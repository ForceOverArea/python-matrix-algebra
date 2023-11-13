#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"

double getValAtIndex(size_t i, size_t j, void *m)
{
    return *getIndex(i, j, ((Matrix *)m));
}

void setValAtIndex(size_t i, size_t j, double val, Matrix *m)
{
    *getIndex(i, j, ((Matrix *)m)) = val;
}

void *getSubset(size_t istart, size_t jstart, size_t iend, size_t jend, void *m)
{
    return (void *)subset(istart, jstart, iend, jend, ((Matrix *)m));
}

void scaleMatrix(double scalar, void *m)
{
    Matrix *mat = ((Matrix *)m); 
    for (size_t i = 0; i < mat->rows; i++)
        scaleRow(scalar, i, mat);
}

bool matrixEquality(void *a, void *b)
{
    return matricesAreEqual((Matrix *)a, (Matrix *)b);
}

void *invertMatrix(void *m)
{
    return (void *)invert((Matrix *)m);
}

void killMatrix(void *m)
{
    free((Matrix *)m);
}

void *matrixProduct(void *a, void *b)
{
    return (void *)multiplyMatrices((Matrix *)a, (Matrix *)b);
}