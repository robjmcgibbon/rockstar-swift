#!/bin/bash -l
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=128
#SBATCH --job-name rockstar_flamingo
#SBATCH --partition cosma8
#SBATCH --output /cosma8/data/dp004/dc-mcgi1/logs/%x.%j
#SBATCH --account dp004
#SBATCH -t 04:00:00
#

module purge
#module load gnu_comp/11.1.0
module load intel_comp/2021.3.0 compiler
module load openmpi/4.1.1
rockstar=/cosma/home/dp004/dc-mcgi1/flamingo/code/rockstar-swift/rockstar-galaxies

# TODO: Make sure snapshot_names.txt exists
#outbase="/snap8/scratch/dp004/dc-mcgi1/rockstar-flamingo/L1000N1800/HYDRO_FIDUCIAL"
outbase="/snap8/scratch/dp004/dc-mcgi1/rockstar-flamingo/L1000N1800/HYDRO_default"
mkdir -p $outbase
if [ -f "${outbase}/auto-rockstar.cfg" ]; then
    echo "Directory already exists. Is this a restart?"
    exit
fi

logfile="${SLURM_JOB_NAME}.${SLURM_JOB_ID}"
mkdir -p logs
ln -s "/cosma8/data/dp004/dc-mcgi1/logs/${logfile}" "./logs/${SLURM_JOB_ID}"
ln -s "/cosma8/data/dp004/dc-mcgi1/logs/${logfile}" "${outbase}/slurm_log"

# Starting master process
$rockstar -c flamingo.cfg &
# If restarting use the following instead
#$rockstar -c "${outbase}restart.cfg" &

# Wait until master process is running
until [ -f "${outbase}/auto-rockstar.cfg" ]; do
    sleep 5
done

# Launch writer processes
# Reader processes will be spawned automatically
mpirun -np $SLURM_NTASKS $rockstar -c "${outbase}/auto-rockstar.cfg"

wait
