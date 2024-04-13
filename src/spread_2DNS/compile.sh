# get the pwd of the script file
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd $DIR

make clean dist
make hd2D

cp hd2D RUNFILES/
cp jobscriptMPI.slurm RUNFILES/
cp input.prm RUNFILES/
cp status.prm RUNFILES/
cp Vis2Db.py RUNFILES/

cp RUNFILES/* $WORK/spread_2DNS/

cd $WORK/spread_2DNS
./clean.sh

