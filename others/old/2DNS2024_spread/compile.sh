make sclean clean dist
make hd2D

cp hd2D RUNFILES/
cp jobscriptMPI.slurm RUNFILES/
cp parameter.inp RUNFILES/
cp status.inp RUNFILES/
cp Vis2Db.py RUNFILES/

cp RUNFILES/* $WORK/TEST1/
