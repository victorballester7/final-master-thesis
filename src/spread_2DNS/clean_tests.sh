# execute ./clean.sh inside all the directories of the form test** to clean the test directories
for dir in test*; do
    cd $dir
    ./clean.sh
    cd ..
done