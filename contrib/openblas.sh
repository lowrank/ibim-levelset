git clone https://github.com/xianyi/OpenBLAS.git

cd OpenBLAS && make USE_OPENMP=1 && make install PREFIX=../../contrib

cd .. && rm -rf OpenBLAS
