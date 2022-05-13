//
// Created by Yimin on 5/8/22.
//

#ifndef LEVELSET_LINALG_LAP_H
#define LEVELSET_LINALG_LAP_H

#include "linalg.h"

using namespace bbfmm;

/*
 * A = P L U
 * 
 * P: pivoting, L U are lower/upper triangle matrices.
 */
void dgetrf(Matrix& A, int* ipiv); // factorization, in-place.
void dgetri(Matrix& A, int* ipiv); // inversion

void dgesv(Matrix& A, Vector& b, int* ipiv); // solve
void dgesv(Matrix& A, Matrix& B, int* ipiv); // solve, overloaded.


void dsyev(Matrix& A, Vector& w);// solve for both eigenvalue and eigevector.

#endif 