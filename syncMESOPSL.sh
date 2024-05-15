# # copy the files to IDRIS
date=$(date +"%Y%m%d") # current date
file1="simple_2DNS"
file2="spread_2DNS"
file3="embarrassed_2DNS"
file4="pointvortices"

# create a tmp directory to store the files
mkdir -p src/${date} src/${file2}/RUNFILES src/${file3}/RUNFILES
cp -r src/${file1} src/${date}
cp -r src/${file2} src/${date}
cp -r src/${file3} src/${date}
cp -r src/${file4} src/${date}

cd src
# create a tar file
tar -cvzf ${date}.tar ${date}/*
# remove the tmp directory
rm -r ${date}

# destination in MESOPSL supercomputer
dest2="~/CODES/"
scp ${date}.tar vballester@mesopsl.obspm.fr:${dest2}

# 'ssh -t vballester@mesopsl.obspm.fr /bin/sh' prevents the loading of the .bash_profile, and so we avoid loading the modules each time
ssh -t vballester@mesopsl.obspm.fr /bin/sh <<EOF
  cd CODES
  # remove current directories
  rm -R -- */
  tar -xvzf ${date}.tar
EOF

# remove the tar file
rm *.tar
