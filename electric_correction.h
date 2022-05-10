//
// Created by Yimin on 5/8/22.
//

#ifndef LEVELSET_ELECTRIC_CORRECT_H
#define LEVELSET_ELECTRIC_CORRECT_H

#include "bbfmm.h"
#include "gmres.h"
#include "levelset.h"
#include "Config.h"

using namespace bbfmm;

void electric_correction(Grid& g, levelset& ls, Surface& surf, Molecule& mol, scalar_t rescale, Config& cfg);


/* ***********************************
 *
 *
 *        auxiliary functions
 *
 *
 * ***********************************/

// biquintic interpolation
scalar_t biquintic_interpolation(Vector& X, Vector& Y, Matrix &wXY, Vector& xy);

// write file into vector
std::istream& operator>> (std::istream& in, std::vector<double>& v);

#endif //LEVELSET_ELECTRIC_CORRECT_H
