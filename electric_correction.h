//
// Created by Yimin on 5/8/22.
//

#ifndef LEVELSET_ELECTRIC_CORRECT_H
#define LEVELSET_ELECTRIC_CORRECT_H

#define CODIM 2
#define SQRT2 1.4142135623730951
#define NFOUR 45
#define BQSIZE 6

#include "bbfmm.h"
#include "gmres.h"
#include "levelset.h"
#include "Config.h"

using namespace bbfmm;

void electric_correction(Grid& g, levelset& ls, Surface& surf, Molecule& mol, scalar_t rescale, Config& cfg,
                         vector<vector<int>> &_contrib_id, 
                         vector<vector<scalar_t>> &K11_contrib_v,
                         vector<vector<scalar_t>> &K21_contrib_v,
                         vector<vector<scalar_t>> &K22_contrib_v
);


/* ***********************************
 *
 *
 *        auxiliary functions 
 *        and classes  
 *
 * ***********************************/

// biquintic interpolation
scalar_t biquintic_interpolation(Vector& X, Vector& Y, Matrix &wXY, Vector& xy);

// write file into vector
std::istream& operator>> (std::istream& in, std::vector<double>& v);

// 2 x 2 matrices COL MAJOR.
// data[0] = a11
// data[1] = a21
// data[2] = a12
// data[3] = a22.
Matrix D0(scalar_t eta, scalar_t kappa_1, scalar_t kappa_2);
Matrix M0(scalar_t kappa_1, scalar_t kappa_2);
Matrix A0(Vector& v1, Vector& v2);

Vector auxiliary_func_0(scalar_t t_, Matrix& D_, Matrix& A_);
scalar_t auxiliary_func_1(scalar_t t_, Matrix& D_, Matrix& A_, Matrix& M_);
scalar_t auxiliary_func_2(scalar_t t_, Matrix& D_, Matrix& A_);

scalar_t w_k0_ptilde1(Vector& alpha_beta, Matrix& AmatFI2, Vector& tmatFI2, 
std::vector<scalar_t>& alpha_all, scalar_t h_alpha, 
std::vector<scalar_t>& beta_all,  scalar_t h_beta,
std::vector<scalar_t>& weight_all,
std::function<scalar_t(scalar_t)> func);


Vector weightsinterp_s2biquintic_ab_plus(Vector &alpha_beta, 
std::vector<scalar_t>& alpha_all, scalar_t h_alpha, 
std::vector<scalar_t>& beta_all,  scalar_t h_beta,
std::vector<scalar_t>& weight_all);


#endif //LEVELSET_ELECTRIC_CORRECT_H
