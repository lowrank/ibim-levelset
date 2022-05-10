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

/// write into vector.
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

    vector<point> source, target;
    vector<scalar_t > weight, normalX, normalY, normalZ;
    vector<scalar_t > curvatures;

    scalar_t dx = ls.dx / rescale;
    /*
     * map all points back to the actual protein surface.
     */
    for (index_t id = 0; id < surf.nodes.size(); ++id) {
        source.push_back(
                {
                    surf.nodes[id].data[0]/rescale,
                    surf.nodes[id].data[1]/rescale,
                    surf.nodes[id].data[2]/rescale
                }
        );

        target.push_back(
                {
                    surf.nodes[id].data[0]/rescale,
                    surf.nodes[id].data[1]/rescale,
                    surf.nodes[id].data[2]/rescale
                }
        );
        weight.push_back(surf.weights[id] * rescale * dx * SQR(dx));

        /*
         * normal vectors do not rescale.
         */
        normalX.push_back(surf.normals[id].data[0]);
        normalY.push_back(surf.normals[id].data[1]);
        normalZ.push_back(surf.normals[id].data[2]);

        curvatures.push_back(surf.curvatures[2 * id] * rescale);
        curvatures.push_back(surf.curvatures[2 * id + 1] * rescale);

    }

    std::cout << std::setw(15) <<"POINTS NUM"  << " " << std::setw(8) << source.size() << std::endl;



    


}


