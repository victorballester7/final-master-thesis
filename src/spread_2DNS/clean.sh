rm -f *.out
rm -rf data/ images/
mkdir -p data images
sed -n '3,6p' jobscriptMPI_IDRIS.slurm
cd data
mkdir -p kspectrum output vectrans EnergyProf EnstrophyProf
