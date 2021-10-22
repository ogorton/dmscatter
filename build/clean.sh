shopt -s extglob
echo "Cleaning build directory..."
rm -r !(clean.sh|build.sh|build-serial.sh) && echo "Done."
