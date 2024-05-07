# execute ./average.py inside all the directories of the form test** to average the data in each test directories
for dir in test*; do
    cd $dir
    echo "Running average.py in $dir"
    python average.py
    cd ..
done
