Compile cppHD.cpp with command: g++ -O3 -fPIC -std=c++20 -pipe -D_FILE_OFFSET_BITS=64 -o cppHD cppHD.cpp
Compile matrixMethod or matrixMethodVer2 with command: g++ -O3 -fPIC -std=c++14 -pipe -D_FILE_OFFSET_BITS=64 -o matrix matrixMethod.cpp -I .\eigen-5.0.0

To run cppHD: cppHD [FASTA FILE].fna [sequence length] [k]

To run matrixMethod or matrixMethodVer2: matrix [FASTA FILE].fna [sequence length] [k]
