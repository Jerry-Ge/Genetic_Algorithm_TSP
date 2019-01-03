#!/bin/bash
# Cmake into build directory
curl -O http://chalonverse.com/435/pa2.tar.gz
tar xzf pa2.tar.gz
echo "Compiling..."
mkdir build
cd build
RELEASE=ON CC=clang CXX=clang++ cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON .. || exit 1
make || exit 1
echo "Running clang-tidy..."
../run-clang-tidy.py -quiet -header-filter=.*/src/.* || exit 1
# Return to root folder (so cwd is correct)
cd ..
# Run student tests
echo "Running student tests..."
timeout 10 build/tests/tests [student]
# Run graded tests
echo "Running graded tests..."
timeout 10 build/tests/tests [graded] || exit 1
# Now try in address sanitizer...
#- mkdir buildAsan
#- cd buildAsan
#- CC=clang-4.0 CXX=clang++-4.0 SANITIZE=ON cmake ..
#- make
#- cd ..
#- echo "Running address sanitizer..."
#- timeout 10 buildAsan/tests/tests [graded]
