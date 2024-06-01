# get the pwd of the script file
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

WORKDIR=/travail/vballester

cd $DIR

make clean dist
make hd2D

cp clean.sh RUNFILES/
cp clean_tests.sh RUNFILES/
cp copy2tests.sh RUNFILES/
cp create_test.sh RUNFILES/
cp hd2D RUNFILES/
cp jobscriptMPI_MESOPSL.slurm RUNFILES/
cp input.prm RUNFILES/
cp status.prm RUNFILES/
cp Vis2Db.py RUNFILES/
cp Vis2Db_bigsize.py RUNFILES/
cp averageEprof.py RUNFILES/
cp movie.py RUNFILES/
cp animate.sh RUNFILES/


cp RUNFILES/* $WORKDIR/spread_2DNS/

cd $WORKDIR/spread_2DNS
./clean.sh
