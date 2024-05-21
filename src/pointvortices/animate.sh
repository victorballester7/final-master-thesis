# animate a series of images
#
# absolute path to the script
# Get the absolute path of the script
SCRIPT_PATH="$(realpath "${BASH_SOURCE[0]}")"
DIR_SCRIPT="$(dirname "$SCRIPT_PATH")"

DIR_IMAGES="$DIR_SCRIPT/../../images/pointvortices"

cd $DIR_IMAGES

rm *.jpg

cd $DIR_SCRIPT

echo "Generating images..."
python plot.py

cd $DIR_IMAGES

# animate the images
echo "Generating animation..."
convert -delay 10 -loop 0 *.jpg ../../videos/pointvortices/animation.mp4
