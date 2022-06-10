//
// Created by lurker on 3/23/17.
//

#ifndef LEVELSET_ELECTRIC_H
#define LEVELSET_ELECTRIC_H

#include "bbfmm.h"
#include "gmres.h"
#include "levelset.h"
#include "Config.h"

void electric(Grid& g, levelset& ls, Surface& surf, Molecule& mol, scalar_t rescale, Config& cfg, 
              vector<vector<int>> &_contrib_id, 
              vector<vector<scalar_t>> &K11_contrib_v,
              vector<vector<scalar_t>> &K21_contrib_v,
              vector<vector<scalar_t>> &K22_contrib_v
);

#endif //LEVELSET_ELECTRIC_H
