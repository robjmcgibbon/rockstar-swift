outbase="/snap8/scratch/dp004/dc-mcgi1/rockstar-flamingo/L1000N0900/HYDRO_FIDUCIAL/"

# Consistent trees doesn't like leading zeros
for i in {0..77}; do
    j=$(printf "%04d" "$i")
    ln -s "${outbase}out_${j}.list" "${outbase}out_${i}.list"
done

perl ../../../scripts/gen_merger_cfg.pl "flamingo.cfg"

