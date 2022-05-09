#include "lapacke_wrapper.h"

using namespace bbfmm;

void dgetrf(Matrix& A, int* ipiv) {
    int ret;
    ret = LAPACKE_dgetrf(LAPACK_COL_MAJOR, A.row(), A.col(), A.data(), A.row(), ipiv);
    assert(ret == 0); // success
}

void dgetri(Matrix& A, int* ipiv) {
    assert(A.row() == A.col());
    int ret;
    ret = LAPACKE_dgetri(LAPACK_COL_MAJOR, A.row(), A.data(), A.row(), ipiv);
    assert(ret == 0);
}


void dgesv(Matrix& A, Vector& b, int* ipiv) {
    assert(A.row() == A.col());
    assert(A.row() == b.row());
    int ret;
    ret = LAPACKE_dgesv(LAPACK_COL_MAJOR, A.row(), 1, A.data(), A.row(), ipiv, b.data(), b.row());
    // A, b are overwritten.
    assert(ret == 0);
}

void dgesv(Matrix& A, Matrix& B, int* ipiv) {
    assert(A.row() == A.col());
    assert(A.row() == B.row());
    int ret;
    ret = LAPACKE_dgesv(LAPACK_COL_MAJOR, A.row(), B.col(), A.data(), A.row(), ipiv, B.data(), B.row());
    // A, b are overwritten.
    assert(ret == 0);
}