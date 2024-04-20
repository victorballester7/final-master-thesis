DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd $DIR

make clean dist
make hd2D

cp clean.sh RUNFILES/
cp hd2D RUNFILES/
cp jobscriptMPI.slurm RUNFILES/
cp input.prm RUNFILES/
cp status.prm RUNFILES/
cp Vis2Db.py RUNFILES/

cp RUNFILES/* $WORK/embarrassed_2DNS/

cd $WORK/embarrassed_2DNS
./clean.sh
