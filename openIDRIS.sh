sshfs uft62hz@jean-zay.idris.fr:/linkhome/rech/genlpe01/uft62hz ~/Desktop/IDRIS/HOME/
if [ $? -eq 0 ]; then
    echo "IDRIS HOME mounted successfully!"
else
    echo "IDRIS HOME mount failed!"
fi

sshfs uft62hz@jean-zay.idris.fr:/gpfswork/rech/gsv/uft62hz ~/Desktop/IDRIS/WORK/
if [ $? -eq 0 ]; then
    echo "IDRIS WORK mounted successfully!"
else
    echo "IDRIS WORK mount failed!"
fi

sshfs uft62hz@jean-zay.idris.fr:/gpfsstore/rech/gsv/uft62hz ~/Desktop/IDRIS/STORE/
if [ $? -eq 0 ]; then
    echo "IDRIS STORE mounted successfully!"
else
    echo "IDRIS STORE mount failed!"
fi
