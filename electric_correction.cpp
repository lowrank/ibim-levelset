#include "electric_correction.h"

using namespace bbfmm;

scalar_t biquintic_interpolation(Vector& X, Vector& Y, Matrix &wXY, Vector& xy) {

    index_t L = 6;

    int* ipiv = (int *)malloc(L * sizeof(int));

    Matrix Vx(L, L), Vy(L, L);
    Vector _x(L)   , _y(L)   ;

    for (int col = 0; col < L; col++) {
        for (int row = 0; row < L; row++) {
            Vx(row, col) = pow( X(col), row );
            Vy(row, col) = pow( Y(col), row );
        }
    }

    for (int row = 0; row < L; row++) {
        _x(row)=pow(xy(0), row);
        _y(row)=pow(xy(1), row);
    }

    dgesv(Vx, _x, ipiv);
    dgesv(Vy, _y, ipiv);

    Vector tmp(L);
    dgemv(1.0, wXY, _y, 0. , tmp) ;

    free(ipiv);

    return ddot(_x, tmp);
}

/// write file.
std::istream& operator>> (std::istream& in, std::vector<double>& v) {
    double d;
    while (in >> d) {
        v.push_back(d);
    }
    return in;
}

/// correction only for singular part within the vacant_radius
void electric_correction(Grid& g, levelset& ls, Surface& surf, Molecule& mol, scalar_t rescale, Config& cfg) {

    // load weights.
    std::ifstream inputFile{"../weights/w.txt"};
    std::vector<double> single_correction_weights;
    inputFile >> single_correction_weights;

}


