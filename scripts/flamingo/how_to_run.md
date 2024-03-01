### flamingo.cfg

Simulation params
- Get softening length from used_params.yml
- Set softening length after multiplying by h
- Repeat for max physical softening

Data
- Set INBASE to snapshot directory
- Set NUM_BLOCKS
- Change SNAPSHOT_NAMES to current directory
- Change OUTBASE

Parallel processing
- If NUM_BLOCKS > NUM_WRITERS then set NUM_READERS
- Set NUM_WRITERS to total number of cores


### snapshot_names.txt

List snapshots to be analysed


### flamingo.sbatch

Job params
- Set number of nodes
- Set runtime
- Set account

Data
- Change outbase


### create_parents.sh

