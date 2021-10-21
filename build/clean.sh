mkdir tmp
cp clean.sh tmp/
cp install.sh tmp/
rm ./* 
rm -r CMakeFiles
cp tmp/* ./
rm -r tmp
