#include "electric_correction.h"

using namespace bbfmm;

/*
 * biquintic interpolation
 */
scalar_t biquintic_interpolation(Vector& X, Vector& Y, Matrix &wXY, Vector& xy) {

    int* ipiv = (int *)malloc(BQSIZE * sizeof(int));

    Matrix Vx(BQSIZE, BQSIZE), Vy(BQSIZE, BQSIZE);
    Vector _x(BQSIZE)   , _y(BQSIZE)   ;

    for (int col = 0; col < BQSIZE; col++) {
        for (int row = 0; row < BQSIZE; row++) {
            Vx(row, col) = pow( X(col), row );
            Vy(row, col) = pow( Y(col), row );
        }
    }

    for (int row = 0; row < BQSIZE; row++) {
        _x(row)=pow(xy(0), row);
        _y(row)=pow(xy(1), row);
    }

    dgesv(Vx, _x, ipiv);
    dgesv(Vy, _y, ipiv);

    Vector tmp(BQSIZE);
    dgemv(1.0, wXY, _y, 0. , tmp) ;

    free(ipiv);

    return ddot(_x, tmp);
}

Vector weightsinterp_s2biquintic_ab_plus(Vector &alpha_beta, 
std::vector<scalar_t>& alpha_all, scalar_t h_alpha, 
std::vector<scalar_t>& beta_all,  scalar_t h_beta,
std::vector<scalar_t>& weight_all) 
{
    int i = int ( ceil( ( alpha_beta(0) - alpha_all[0]) / h_alpha )   ) ;
    int j = int ( ceil( ( alpha_beta(1)  - beta_all[0] ) / h_beta )   ) ;

    int na = alpha_all.size();
    int nb = beta_all.size();

    i = min(int(alpha_all.size() - 4), max(i, 2));
    j = min(int(beta_all.size()  - 4), max(j, 2));

    Vector X(BQSIZE);
    Vector Y(BQSIZE);

    for (int b_id = 0; b_id < BQSIZE; b_id++) {
        X(b_id) = alpha_all[i-2 + b_id]; // i-2 to i+3
        Y(b_id) = beta_all [j-2 + b_id];
    }

    Vector ret(NFOUR);
    Matrix wXY(BQSIZE, BQSIZE);

    for (int k_id; k_id < NFOUR; k_id++) {
        for (int row = 0; row < BQSIZE; row++) {
            for (int col = 0; col < BQSIZE; col++) {
                index_t ind = k_id + (i-2 + row) * NFOUR + (j-2 + col) * NFOUR * na;
                wXY(row, col) = weight_all[ind];
            }
        }

        // std::cout << wXY << std::endl;
        // std::cout << "alpha beta " << alpha_beta(0) << " " << alpha_beta(1) << std::endl;
        

        Vector v(2); v = alpha_beta; // copy
        ret(k_id) = biquintic_interpolation(X, Y, wXY, v);

        // std::cout << ret(k_id) << std::endl;
    }

    return ret;
}

scalar_t w_k0_ptilde1(Vector& alpha_beta, Matrix& AmatFI2, Vector& tmatFI2, 
std::vector<scalar_t>& alpha_all, scalar_t h_alpha, 
std::vector<scalar_t>& beta_all,  scalar_t h_beta,
std::vector<scalar_t>& weight_all,
std::function<scalar_t(scalar_t)> func) {
    Vector omall_w1 = weightsinterp_s2biquintic_ab_plus(alpha_beta, alpha_all, h_alpha, beta_all, h_beta, weight_all);
    
    // std::cout << omall_w1 << std::endl;

    Vector f_tmatFI2(tmatFI2.row());
    for (int i = 0; i < tmatFI2.row(); ++i) {
        f_tmatFI2(i) = func(tmatFI2(i));
    }

    int* ipiv = (int *)malloc(tmatFI2.row() * sizeof(int));
    Matrix X(AmatFI2.row(), AmatFI2.col());
    X = AmatFI2;
    dgesv(X, f_tmatFI2, ipiv);
    free(ipiv);

    // std::cout << "f"<<f_tmatFI2 << std::endl;

    return ddot(f_tmatFI2, omall_w1);
}

/// write into vector.
std::istream& operator>> (std::istream& in, std::vector<double>& v) {
    double d;
    while (in >> d) {
        v.push_back(d);
    }
    return in;
}

Matrix D0(scalar_t eta, scalar_t kappa_1, scalar_t kappa_2) {
    Matrix ret(CODIM, CODIM);
    ret(0, 0) = 1.0 / (1 + kappa_1 * eta);
    ret(1, 1) = 1.0 / (1 + kappa_2 * eta);
    return ret;
}

Matrix M0(scalar_t kappa_1, scalar_t kappa_2) {
    Matrix ret(CODIM, CODIM);
    ret(0, 0) = kappa_1;
    ret(1, 1) = kappa_2;
    return ret;
}

Matrix A0(Vector& v1, Vector& v2) {
    Matrix ret(CODIM, CODIM);
    ret(0, 0) = v1(0);
    ret(1, 0) = v1(1);
    ret(0, 1) = v2(0);
    ret(1, 1) = v2(1);
    return ret;
}

Vector auxiliary_func_0(scalar_t t_, Matrix& D_, Matrix& A_) {
    Vector v(2), u(2);
    v(0) = cos(t_); v(1) = sin(t_);
    dgemv(1.0, A_, v, 0., u); // u = A v
    dgemv(1.0, D_, u, 0., v); // v = D u
    return v;
}

scalar_t auxiliary_func_1(scalar_t t_, Matrix& D_, Matrix& A_, Matrix& M_) {
    Vector q(2);
    Vector v = auxiliary_func_0(t_, D_, A_);
    dgemv(1.0, M_, v, 0., q); // q = M v
    scalar_t num = ddot(q, v);
    scalar_t den = pow( nrm2(v), 3 );

    return 0.5 * num / den;

    //  0.5*(dot(D0*Amat*[cos(t);sin(t)], Mmat*D0*Amat*[cos(t);sin(t)]))/(norm(D0*Amat*[cos(t);sin(t)])^3)
}

scalar_t auxiliary_func_2(scalar_t t_, Matrix& D_, Matrix& A_) {
    Vector v = auxiliary_func_0(t_, D_, A_);
    return 1.0 / nrm2( v );
}

/// correction only for singular part within the vacant_radius
void electric_correction(Grid& g, levelset& ls, Surface& surf, Molecule& mol, scalar_t rescale, Config& cfg) {

    // load weights.
    std::ifstream alpha_inputFile{"../weights/alpha_all.txt"};
    std::vector<double> alpha_all;
    alpha_inputFile >> alpha_all;

    std::ifstream beta_inputFile{"../weights/beta_all.txt"};
    std::vector<double> beta_all;
    beta_inputFile >> beta_all;

    std::ifstream weight_inputFile{"../weights/weights_all.txt"};
    std::vector<double> weight_all;
    weight_inputFile >> weight_all;

    index_t Na = alpha_all.size();
    index_t Nb = beta_all.size(); 
    index_t Nw = weight_all.size();

    scalar_t h_alpha = alpha_all[1] - alpha_all[0];
    scalar_t h_beta  = beta_all[1]  - beta_all[0];

    // weights_all_shape = 45, 101, 101
    Matrix AmatFI2(NFOUR, NFOUR);
    setValue(AmatFI2, 1.0);

    Vector tmatFI2(NFOUR);
    for (int t_id = 0; t_id < NFOUR; t_id++) {
        tmatFI2(t_id) = t_id * M_PI/NFOUR;
    }

    scalar_t hmatFI2 = M_PI/NFOUR;

    for (int row = 0; row < NFOUR; row++) {
        for (int col = 1; col <= NFOUR / 2; col++) {
            AmatFI2(row, 2 * col - 1) = cos(2 * col * tmatFI2(row));
            AmatFI2(row, 2 * col    ) = sin(2 * col * tmatFI2(row));
        }
    }

    vector<point> source, target;
    vector<scalar_t > weight, normalX, normalY, normalZ;
    vector<scalar_t > curvatures;

    scalar_t dx = ls.dx / rescale;

    vector<vector<int>> _contrib_id(surf.nodes.size());

    vector<vector<scalar_t>> K11_contrib_v(surf.nodes.size());
    vector<vector<scalar_t>> K12_contrib_v(surf.nodes.size());
    vector<vector<scalar_t>> K21_contrib_v(surf.nodes.size());
    vector<vector<scalar_t>> K22_contrib_v(surf.nodes.size());
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
    
    std::cout << std::setw(15) <<"POINTS NUM"  << " " << std::setw(8) << source.size()  << std::endl;

    for (index_t id = 0; id < source.size(); id++) {
        scalar_t theta = acos(normalZ[id] / norm(surf.normals[id]));
        scalar_t phi   = atan2(normalY[id] , normalX[id]);

        scalar_t kappa_1 = curvatures[2 * id];
        scalar_t kappa_2 = curvatures[2 * id + 1];

        Matrix M = M0(kappa_1, kappa_2);
        Matrix A = A0(surf.eigenvectors[2 * id], surf.eigenvectors[2 * id + 1]);

        vector<scalar_t> ALPHA;
        vector<scalar_t> BETA;
        vector<scalar_t> DIST;
        vector<int>      ID;

        // get loc of projected node, round towards direction (-1,-1,-1).
        index_t node_x = (index_t)((surf.nodes[id].data[0] - ls.sx) / ls.dx);
        index_t node_y = (index_t)((surf.nodes[id].data[1] - ls.sy) / ls.dx);
        index_t node_z = (index_t)((surf.nodes[id].data[2] - ls.sz) / ls.dx);
        
        // decide rotation.
        if (fabs(tan(theta)) >= SQRT2 && fabs(tan(phi)) >= 1 ) {
        // 1. Tilted. y-planes
            for (int cur_id = -ls.thickness; cur_id <= ls.thickness + 1; cur_id++) {
                scalar_t cur_y = ls.sy + (node_y + cur_id) * ls.dx;
                scalar_t scale = ( cur_y - surf.nodes[id].data[1] ) / normalY[id];
                scalar_t cur_z = surf.nodes[id].data[2] + scale * normalZ[id];
                scalar_t cur_x = surf.nodes[id].data[0] + scale * normalX[id];

                // find closet node.
                index_t cur_node_y = node_y + cur_id;
                index_t cur_node_z = round((cur_z - ls.sz) / ls.dx) ;
                index_t cur_node_x = round((cur_x - ls.sx) / ls.dx) ;

                scalar_t dist  = g.get(cur_node_x, cur_node_y, cur_node_z);

                if (fabs( dist ) < ls.thickness * ls.dx) {
                    scalar_t alpha = (cur_z - ls.sz) / ls.dx - cur_node_z;
                    scalar_t beta  = (cur_x - ls.sx) / ls.dx - cur_node_x;

                    ALPHA.push_back(alpha);
                    BETA.push_back(beta);

                    DIST.push_back(dist / rescale);
                    ID.push_back(  cur_node_x * ls.Ny * ls.Nz + cur_node_y * ls.Nz + cur_node_z  );

                }                
            }
        }
        else if (fabs(tan(theta)) >= SQRT2 && fabs(tan(phi)) < 1) {
        // 2. Tilted. x-planes.
            for (int cur_id = -ls.thickness; cur_id <= ls.thickness + 1; cur_id++) {
                scalar_t cur_x = ls.sx + (node_x + cur_id) * ls.dx;
                scalar_t scale = ( cur_x - surf.nodes[id].data[0] ) / normalX[id];
                scalar_t cur_y = surf.nodes[id].data[1] + scale * normalY[id];
                scalar_t cur_z = surf.nodes[id].data[2] + scale * normalZ[id];


                index_t cur_node_x = node_x + cur_id;
                index_t cur_node_y = round((cur_y - ls.sy) / ls.dx) ;
                index_t cur_node_z = round((cur_z - ls.sz) / ls.dx) ;

                scalar_t dist  = g.get(cur_node_x, cur_node_y, cur_node_z);

                if (fabs( dist ) < ls.thickness * ls.dx) {
                    scalar_t alpha = (cur_y - ls.sy) / ls.dx - cur_node_y;
                    scalar_t beta  = (cur_z - ls.sz) / ls.dx - cur_node_z;

                    ALPHA.push_back(alpha);
                    BETA.push_back(beta);

                    DIST.push_back(dist / rescale);
                    ID.push_back(  cur_node_x * ls.Ny * ls.Nz + cur_node_y * ls.Nz + cur_node_z  );
                }

            }
        }
        else {
        // 3. If not too tilted, use default z-planes
            for (int cur_id = -ls.thickness; cur_id <= ls.thickness + 1; cur_id++) {
                scalar_t cur_z = ls.sz + (node_z + cur_id) * ls.dx;
                scalar_t scale = ( cur_z - surf.nodes[id].data[2] ) / normalZ[id];
                scalar_t cur_x = surf.nodes[id].data[0] + scale * normalX[id];
                scalar_t cur_y = surf.nodes[id].data[1] + scale * normalY[id];

                index_t cur_node_z = node_z + cur_id;
                index_t cur_node_x = round((cur_x - ls.sx) / ls.dx) ;
                index_t cur_node_y = round((cur_y - ls.sy) / ls.dx) ;

                scalar_t dist  = g.get(cur_node_x, cur_node_y, cur_node_z);

                if (fabs( dist ) < ls.thickness * ls.dx) {
                    scalar_t alpha = (cur_x - ls.sx) / ls.dx - cur_node_x;
                    scalar_t beta  = (cur_y - ls.sy) / ls.dx - cur_node_y;

                    ALPHA.push_back(alpha);
                    BETA.push_back(beta);

                    DIST.push_back(dist / rescale);
                    ID.push_back(  cur_node_x * ls.Ny * ls.Nz + cur_node_y * ls.Nz + cur_node_z  );
                }

            } // end for
        } // end else

        //get 4 values.
        
        scalar_t kappa =atof(cfg.options["solvent_kappa"].c_str());
        scalar_t dE = atof(cfg.options["solvent_dE"].c_str());
        scalar_t dI = atof(cfg.options["solvent_dI"].c_str());

        scalar_t secondt_const = 0.25/M_PI;
        scalar_t epsl_ratio = dE/dI;

        if (id % 5000 == 0) std::cout << id << std::endl;

        for (int _s = 0; _s < ALPHA.size(); ++ _s) {

            auto dist = DIST[_s];
            auto D0_now = D0(dist, kappa_1, kappa_2);

            auto lfun_1 = [&](scalar_t t) { return auxiliary_func_1(t, D0_now, A, M);};
            auto lfun_2 = [&](scalar_t t) { return auxiliary_func_2(t, D0_now, A);};

            Vector alpha_beta(2);
            alpha_beta(0) = ALPHA[_s];
            alpha_beta(1) = BETA[_s];

            auto wjj = w_k0_ptilde1(alpha_beta, 
            AmatFI2, tmatFI2, 
            alpha_all, h_alpha, 
            beta_all,  h_beta,
            weight_all,
            lfun_1);
            
            auto w21 = w_k0_ptilde1(alpha_beta, 
            AmatFI2, tmatFI2, 
            alpha_all, h_alpha, 
            beta_all,  h_beta,
            weight_all,
            lfun_2
            );


            _contrib_id[id].push_back(surf.mapping[ID[_s]]);
                
            K11_contrib_v[id].push_back( wjj * weight[surf.mapping[ID[_s]]] * secondt_const * (1 - epsl_ratio) );
            K22_contrib_v[id].push_back( wjj * weight[surf.mapping[ID[_s]]] * secondt_const * (1 - 1/epsl_ratio) );
            K21_contrib_v[id].push_back( w21*  weight[surf.mapping[ID[_s]]] * 0.5 * SQR(kappa) );
            K12_contrib_v[id].push_back( secondt_const * kappa * weight[surf.mapping[ID[_s]]] );
        }

// ij0_PB_a(t,D0,Amat,Mmat) = 0.5*(dot(D0*Amat*[cos(t);sin(t)], Mmat*D0*Amat*[cos(t);sin(t)]))/(norm(D0*Amat*[cos(t);sin(t)])^3)

// wjj = w_k0_ptilde1([α;β]; lfun=(t->ij0_PB_a(t,D0_now,Amat, Mmat)));
// w21 = w_k0_ptilde1([α;β]; lfun=(t->1/norm(D0_now*Amat*[cos(t);sin(t)])))
// w_K11_single[ m ][ itmp ] = wjj * secondt_const * (1-epsl_ratio)
// w_K22_single[ m ][ itmp ] = wjj * secondt_const * (1- 1/epsl_ratio)
// w_K21_single[ m ][ itmp ] = w21 * secondt_const * 0.5 * kappa_val^2
// w_K12_single[ m ][ itmp ] = secondt_const * kappa_val # no actual weight computation necessary
    } // end for
}


