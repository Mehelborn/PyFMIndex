git clone https://github.com/mehelborn/PyFMIndex.git
git submodule update --init --recursive --remote


cd lib/AvxWindowFmIndex/
cmake .
make
