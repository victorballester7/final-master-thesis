#add colors to echo's
YELLOW='\033[1;33m'
BLUE='\033[1;34m'
GREEN='\033[1;32m'
RED='\033[1;31m'
RESET='\033[0m'

echo -e "${YELLOW}Preparing...${RESET}"

CODE=disk
DATA_DIR=data/pointvortices/$CODE

# create necessary folders
rm -rf $DATA_DIR
mkdir -p $DATA_DIR/positions $DATA_DIR/EnergyProf $DATA_DIR/EnergyFlux $DATA_DIR/NumVortices

cd src/pointvortices/$CODE

# compile the code
echo -e "${YELLOW}Compiling...${RESET}"
make
if [ $? -ne 0 ]; then
  echo -e "${RED}Compilation failed!${RESET}"
  exit 1
fi
echo -e "${GREEN}Compilation done!${RESET}"
# run the code

echo -e "${YELLOW}Running...${RESET}"
cd ../../..
./bin/pointvortices/$CODE/main
if [ $? -ne 0 ]; then
  echo -e "${RED}Running failed!${RESET}"
  exit 1
fi
echo -e "${GREEN}Running done!${RESET}"
