# animate a series of images
#
# absolute path to the script
# Get the absolute path of the script
SCRIPT_PATH="$(realpath "${BASH_SOURCE[0]}")"
DIR_SCRIPT="$(dirname "$SCRIPT_PATH")"

if [ $1 -eq 0 ]; then
  echo "Usage: $0 <type>"
  echo "type: disk, dipoles"
  exit 1
fi

DIR_IMAGES="$DIR_SCRIPT/../../images/pointvortices"

if [ $1 == "disk" ]; then
  DIR_IMAGES="$DIR_IMAGES/disk"
elif [ $1 == "dipoles" ]; then
  DIR_IMAGES="$DIR_IMAGES/dipoles"
else
  echo "Invalid type: $1"
  exit 1
fi

cd $DIR_IMAGES

rm *.jpg

cd $DIR_SCRIPT

echo "Generating images..."
python plot.py $1

cd $DIR_IMAGES

# animate the images
echo "Generating animation..."
convert -delay 10 -loop 0 *.jpg ../../videos/pointvortices/$1/animation.mp4