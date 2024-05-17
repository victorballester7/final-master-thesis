sshfs vballester@mesopsl.obspm.fr:/obs/vballester ~/Desktop/MESOPSL/HOME/
if [ $? -eq 0 ]; then
    echo "MESOPSL HOME mounted successfully!"
else
    echo "MESOPSL HOME mount failed!"
fi

sshfs vballester@mesopsl.obspm.fr:/travail/vballester/ ~/Desktop/MESOPSL/WORK/
if [ $? -eq 0 ]; then
    echo "MESOPSL WORK mounted successfully!"
else
    echo "MESOPSL WORK mount failed!"
fi
