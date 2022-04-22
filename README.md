# Implicit Boundary Integral Method with FMM

## Install OpenBlas
- First, in directory ``contrib``, run ``bash openblas.sh`` to install OpenBLAS (with OPENMP support).
## Compile IBIM-LEVELSET
- ``mkdir build && cd build``
- ``cmake ..`` 
- ``make``

## Run Experiments
There is an example config file in ``./data`` directory. Setting the numeric values and run with
- ./levelset ../data/input.cfg

