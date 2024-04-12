#add colors to echo's
YELLOW='\033[1;33m'
BLUE='\033[1;34m'
GREEN='\033[1;32m'
RED='\033[1;31m'
RESET='\033[0m'

pwd=$(pwd)

myNS="2DNS2024_simple"

echo -e "${YELLOW}Preparing...${RESET}"
cd src/$myNS
# remove old files
make sclean clean dist > /dev/null 2>&1
# create necessary folders
mkdir -p $pwd/data $pwd/images $pwd/data/kspectrum $pwd/data/mspectrum $pwd/data/output $pwd/data/vectrans $pwd/data/EnergyProf $pwd/data/EnstrophyProf


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