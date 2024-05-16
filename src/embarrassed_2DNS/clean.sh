rm -f *.out
rm -rf data/ images/
mkdir -p data images data/average
sed -n '3,6p' jobscriptMPI*
cd data
mkdir -p kspectrum output vectrans EnergyProf EnstrophyProf average/kspectrum average/EnergyProf average/EnstrophyProf
