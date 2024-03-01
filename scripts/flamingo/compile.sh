module purge
module load intel_comp/2019
module load hdf5

make with_hdf5
make parents

