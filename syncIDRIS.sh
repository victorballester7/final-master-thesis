# # copy the files to IDRIS
date=$(date +"%Y%m%d") # current date
file1="simple_2DNS"
file2="spread_2DNS"
file3="embarrassed_2DNS"
file4="pointvortices"

# working directory in IDRIS supercomputer
WORKIDRIS="/gpfswork/rech/gsv/uft62hz"

# create a tmp directory to store the files
mkdir src/${date}
cp -r src/${file1} src/${date}
cp -r src/${file2} src/${date}
cp -r src/${file3} src/${date}
cp -r src/${file4} src/${date}

cd src
# remove the previous tar files
rm *.tar
# create a tar file
tar -cvzf ${date}.tar ${date}/*
# remove the tmp directory
rm -r ${date}

# destination in gauss computer
dest1="~/Desktop/"
scp ${date}.tar victor@gauss:${dest1}

# destination in IDRIS supercomputer
dest2="~/CODES/"

# we enter in gauss and then in IDRIS
ssh victor@gauss << EOF
    scp ${dest1}/${date}.tar uft62hz@jean-zay.idris.fr:${dest2}
    ssh uft62hz@jean-zay.idris.fr << INNER_EOF
        cd CODES
        # remove current directories
        rm -R -- */
        tar -xvzf ${date}.tar
INNER_EOF
EOF