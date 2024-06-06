DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

WORKDIR=$WORK

cd $DIR

make clean dist
make hd2D

cp clean.sh RUNFILES/
cp clean_tests.sh RUNFILES/
cp average_tests.sh RUNFILES/
cp copy2tests.sh RUNFILES/
cp create_test.sh RUNFILES/
cp hd2D RUNFILES/
cp jobscriptMPI_IDRIS.slurm RUNFILES/
cp input.prm RUNFILES/
cp status.prm RUNFILES/
cp Vis2Db.py RUNFILES/
cp average.py RUNFILES/
cp averageEprof.py RUNFILES/

cp RUNFILES/* $WORKDIR/embarrassed_2DNS

cd $WORKDIR/embarrassed_2DNS

./clean.sh
