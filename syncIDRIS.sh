# # copy the files to IDRIS
date=$(date +"%Y%m%d")
file1="simple_2DNS"
file2="spread_2DNS"
file3="embarrassed_2DNS"
file4="pointvortices"

WORKIDRIS="/gpfswork/rech/gsv/uft62hz"

# mv src/${file1}** src/${file1}_${date}
# mv src/${file2}** src/${file2}_${date}
# mv src/${file3}** src/${file3}_${date}
# mv src/${file4}** src/${file4}_${date}

mkdir src/${date}
cp -r src/${file1} src/${date}
cp -r src/${file2} src/${date}
cp -r src/${file3} src/${date}
cp -r src/${file4} src/${date}

cd src
rm *.tar
tar -cvzf ${date}.tar ${date}/*
rm -r ${date}

# # destination in gauss
dest1="~/Desktop/"
scp ${date}.tar victor@gauss:${dest1}

dest2="~/CODES/"
# ssh victor@gauss "scp ${dest1}/${date}.tar uft62hz@jean-zay.idris.fr:${dest2}"

ssh victor@gauss << EOF
    scp ${dest1}/${date}.tar uft62hz@jean-zay.idris.fr:${dest2}
    ssh uft62hz@jean-zay.idris.fr << INNER_EOF
        cd CODES
        # remove current directories
        rm -R -- */
        tar -xvzf ${date}.tar
INNER_EOF
EOF