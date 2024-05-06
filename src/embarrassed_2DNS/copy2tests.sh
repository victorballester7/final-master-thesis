# copy all the necessary files to the test directories
files="clean.sh hd2D jobscriptMPI.slurm status.prm images/ data/ input.prm Vis2Db.py average.py"
for dir in test*; do
  cp -r $files $dir
done
