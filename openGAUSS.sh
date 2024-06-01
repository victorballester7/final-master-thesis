sshfs victor@gauss:/home/victor/Desktop/ ~/Desktop/GAUSS/Desktop/
if [ $? -eq 0 ]; then
    echo "Gauss Desktop mounted successfully!"
else
    echo "Gauss Desktop mount failed!"
fi
