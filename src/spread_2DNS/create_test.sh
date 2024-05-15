# creates a test folder with the necessary files for the test to run

if [ -z "$1" ]; then
  echo "Usage: $0 <test_foldername>"
  exit 1
fi

mkdir $1

# copy all the necessary files to the test directories
files="clean.sh hd2D jobscriptMPI.slurm status.prm images/ data/ input.prm Vis2Db.py"
cp -r $files $1
