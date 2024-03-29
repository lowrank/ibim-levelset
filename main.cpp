#define DISP
#define VIEW
#define GRID
#include <iostream>
#include "electric.h"
#include "electric_correction.h"

#undef VIEW
#undef GRID

#ifdef VIEW
#include "view.h"
#endif

int main(int argc, char* argv[]) {


#ifdef RUN_OMP
    omp_set_num_threads(omp_get_max_threads());
#endif

    if (argc <= 1) {
        std::cout << "USE " <<argv[0] << " PATH_OF_CONFIG_FILE " << std::endl;
        exit(0);
    }

    Config cfg;
    std::ifstream cfgFile;cfgFile.open(argv[1], std::ifstream::in);
    cfg.parse(cfgFile);cfgFile.close();

    cfg.print();

    Molecule mol; mol.load(cfg.options["pqr_file"]);

    scalar_t s = mol.centralize(200.0);
    mol.getCenter();
    scalar_t pr = 1.4 * s;

    scalar_t grid_lo = -300.0, grid_hi = 300;
    index_t size = atoi(cfg.options["grid_size"].c_str());
    scalar_t dx = (grid_hi - grid_lo) / scalar_t(size);

    std::cout << std::setw(15) << "h" << " " << std::setw(8) << dx / s << " Angstroms" << std::endl;
    std::cout << std::setw(15) << "s" << " " << std::setw(8) << s << " Rescale" << std::endl;

    levelset ls(size, size, size, 11, grid_lo, grid_lo, grid_lo, dx, cfg);

    Grid g(-2.0 * size, ls.Nx, ls.Ny, ls.Nz);
    Grid phi0(0., ls.Nx, ls.Ny, ls.Nz);

    RUN("EXPAND", ls.expand(mol, g, pr));
    RUN("INWARD", ls.evolve(g, 1.4, s, 0.2));

    // flip the sign. Outside is positive, then gradients are pointing outside.
    g *= -1.0;
    phi0 = g;

    // if a point is far, the distance is fixed explicitly.
    // RUN("PRESET", ls.setExterior(mol, g, pr));

    // only reinitialization on nearby points?
    RUN("REINIT", ls.reinitialize(g, phi0, atoi(cfg.options["reinit_step"].c_str()), 1, 0.5, pr));

    ls.setInclusion(g);


    /*
     *  surface information.
     *
     */
    Surface surf(g, ls, ls.thickness * ls.dx);

#ifdef GRID
    g.output("../data/output.grid");
    surf.output("../data/output.node");
#endif


#ifdef VIEW
    view& v = view::getInstance(2);
    v.loadLevelSet(g, ls, surf);
    v.run();
#endif
    
    vector<vector<int>> _contrib_id;
    vector<vector<scalar_t>> K11_contrib_v;
    vector<vector<scalar_t>> K21_contrib_v;
    vector<vector<scalar_t>> K22_contrib_v;
    
    electric_correction(g, ls, surf, mol, s, cfg, 
                            _contrib_id, 
                            K11_contrib_v, 
                            K21_contrib_v,
                            K22_contrib_v);
    
     
    electric(g, ls, surf, mol, s, cfg,
                        _contrib_id, 
                        K11_contrib_v, 
                        K21_contrib_v,
                        K22_contrib_v);
    return 0;
}