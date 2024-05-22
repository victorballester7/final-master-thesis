#add colors to echo's
YELLOW='\033[1;33m'
BLUE='\033[1;34m'
GREEN='\033[1;32m'
RED='\033[1;31m'
RESET='\033[0m'

echo -e "${YELLOW}Preparing...${RESET}"

# create necessary folders
rm -rf data/pointvorticesDipoles
mkdir -p data data/pointvorticesDipoles data/pointvorticesDipoles/positions data/pointvorticesDipoles/EnergyProf data/pointvorticesDipoles/EnergyFlux data/pointvorticesDipoles/NumVortices

cd src/pointvorticesDipoles

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
cd ../..
./bin/pointvorticesDipoles/main
if [ $? -ne 0 ]; then
  echo -e "${RED}Running failed!${RESET}"
  exit 1
fi
echo -e "${GREEN}Running done!${RESET}"