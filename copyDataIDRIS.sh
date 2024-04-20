# # copy the files to IDRIS
file1="simple_2DNS"
file2="spread_2DNS"
file3="embarrassed_2DNS"
file4="pointvortices"

# working directory in IDRIS supercomputer
WORKIDRIS="/gpfswork/rech/gsv/uft62hz"

# we enter in gauss and then in IDRIS
ssh victor@gauss << EOF
    ssh uft62hz@jean-zay.idris.fr << INNER_EOF
        cd ${WORKIDRIS}/${file2}
        tar -cvzf ${file2}.tar data/*
        cd ${WORKIDRIS}/${file3}
        tar -cvzf ${file3}.tar data/*
INNER_EOF
    scp uft62hz@jean-zay.idris.fr:${WORKIDRIS}/${file2}/${file2}.tar ~/Desktop/${file2}.tar
    scp uft62hz@jean-zay.idris.fr:${WORKIDRIS}/${file3}/${file3}.tar ~/Desktop/${file3}.tar
EOF

rm -r data/${file2} data/${file3}
mkdir -p data/${file2} data/${file3}

scp -r victor@gauss:~/Desktop/${file2}/data/* data/${file2}/
scp -r victor@gauss:~/Desktop/${file3}/data/* data/${file3}/

cd data/${file2}
tar -xvzf ${file2}.tar
cd ../${file3}
tar -xvzf ${file3}.tar


# remove the tar files
ssh victor@gauss << EOF
    ssh uft62hz@jean-zay.idris.fr << INNER_EOF
        rm ${WORKIDRIS}/${file2}/${file2}.tar
        rm ${WORKIDRIS}/${file3}/${file3}.tar
INNER_EOF
    rm ~/Desktop/${file2}.tar ~/Desktop/${file3}.tar
EOF