#add colors to echo's
YELLOW='\033[1;33m'
BLUE='\033[1;34m'
GREEN='\033[1;32m'
RED='\033[1;31m'
RESET='\033[0m'

pwd=$(pwd)
data=$pwd/data/simple_2DNS

myNS="simple_2DNS"

echo -e "${YELLOW}Preparing...${RESET}"

cd src/$myNS

# check if we start a new run from 0 or not
# Get the first line of the file
first_line=$(head -n 1 status.prm)
# Extract the first string
first_string=$(echo "$first_line" | awk '{print $1}') 
# Convert the string to an integer
n=$(printf "%.0f" "$first_string") # = last output file

if [ $n -eq 0 ]; then
  # remove old files
  make sclean clean dist > /dev/null 2>&1
else
  make clean dist > /dev/null 2>&1
fi

# create necessary folders
mkdir -p $data $pwd/images $data/kspectrum $data/output $data/vectrans $data/EnergyProf $data/EnstrophyProf


# compile the code
echo -e "${YELLOW}Compiling...${RESET}"
make hd2D
if [ $? -ne 0 ]; then
  echo -e "${RED}Compilation failed!${RESET}"
  exit 1
fi
echo -e "${GREEN}Compilation done!${RESET}"
# run the code

echo -e "${YELLOW}Running...${RESET}"
../../bin/$myNS/hd2D
if [ $? -ne 0 ]; then
  echo -e "${RED}Running failed!${RESET}"
  exit 1
fi
echo -e "${GREEN}Running done!${RESET}"