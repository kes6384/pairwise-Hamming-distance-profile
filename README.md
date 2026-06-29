Compile cppHD.cpp with command: g++ -O3 -fPIC -std=c++20 -pipe -D_FILE_OFFSET_BITS=64 -o cppHD cppHD.cpp  
Compile matrixMethod or matrixMethodVer2 with command: g++ -O3 -fPIC -std=c++14 -pipe -D_FILE_OFFSET_BITS=64 -o matrix matrixMethod.cpp -I .\eigen-5.0.0  
Compile cppHD_openMP.cpp with command: g++ -O3 -fPIC -fopenmp -std=c++20 -pipe -D_FILE_OFFSET_BITS=64 -o cppHDOMP cppHD_openMP.cpp

To run cppHD: cppHD [FASTA FILE].fna [sequence length] [k]  
To run matrixMethod or matrixMethodVer2: matrix [FASTA FILE].fna [sequence length] [k]
To run cppHDOMP: cppHDOMP [FASTA FILE].fna [sequence length] [k] [numThreads]
