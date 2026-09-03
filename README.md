# PAIRWISE HAMMING DISTANCE PROFILE  
  
This tool offers several methods to compute the pairwise Hamming distance profile of a given sequence. It is designed for use with FASTA files and $k \leq 128$. The tool can be compiled and run using the command line. The results are output to a user-specified text file.  
  
The tool can be used to compute the full Hamming distance profile of a sequence, or to compute a sketch of the profile for faster (but less accurate) results. Support for multithreading is also provided.  
  
All of the provided methods are contained within profile.cpp. Required open source libraries are provided in the repo. MurmurHash repo: https://github.com/aappleby/smhasher/tree/master  

## REQUIREMENTS  

$C++ \geq 20$  
$Windows OS \geq 11$  
  
The tool has not been successfully tested with other OS or compilers.  
  
## INSTALL AND COMPILING  

```
git clone https://github.com/kes6384/pairwise-Hamming-distance-profile.git  
cd pairwise-Hamming-distance-profile  
g++ -O3 -fPIC -fopenmp -std=c++20 -pipe -D_FILE_OFFSET_BITS=64 -o profile profile.cpp
```
  
## RUNNING  

Once compiled, the tool can be used via the command line using the general format below:  

```
cd pairwise-Hamming-distance-profile  
./profile [INPUT FILE].fna [OUTPUT FILE].txt [SEQUENCE LENGTH] [k] [ROW SAMPLE RATE] [COLUMN SAMPLE RATE] [NUMBER OF THREADS]  
```

The sampling rate and number of threads arguments are optional and can be used together. Omit the sampling rate arguments to generate the full profile. To generate a sketch instead of the full profile, you must include both a row and column sampling rate.  

`-L`: length of the input sequence ... 
INPUT FILE - FASTA file that contains the sequence you want to analyze  
OUTPUT FILE - text file where results will be written to  
SEQUENCE LENGTH - length of the provided sequence to analyze (can be smaller than the length of the entire sequence in INPUT FILE)  
k - length of k-mers, must be $k 1 \leq k \geq 128$ 
ROW SAMPLE RATE - optional argument, rate at which to sample first k-mer; default is 1 (full profile)  
COLUMN SAMPLE RATE - optional argument, rate at which to sample second k-mer; default is 1 (full profile)  
NUMBER OF THREADS - optional argument, use this to run a multithreaded version of the tool; default is 1  
  
### EXAMPLES  
The FASTA file used for the following examples and the results of running the following commands can be found in the folder "test". 

To use one thread to calculate the full Hamming distance profile of 32-mers for the first 10,000 characters of the sequence found in test.fna and store the results in results.txt, any of the following three commands will work:
```
./profile test.fna results.txt 10000 32
```
```
./profile test.fna results.txt 10000 32 1.0 1.0
```
```
./profile test.fna results.txt 10000 32 1
```
```
./profile test.fna results.txt 10000 32 1.0 1.0 1
```
  
To use one thread to sample 32-mer pairs at a rate of 0.25, any of the following commands will work:
```
./profile test.fna results.txt 10000 32 0.5 0.5
```
```
./profile test.fna results.txt 10000 32 0.5 0.5 1
```
  
To calculate the full profile using 4 threads, any of the following commands will work:
```
./profile test.fna results.txt 10000 32 4
```
```
./profile test.fna results.txt 10000 32 1.0 1.0 4
```
  
To use 4 threads to sample 32-mer pairs at a rate of 0.25, use the following command:
```
./profile test.fna results.txt 10000 32 0.5 0.5 4
```
  
## OUTPUT  
  
The results of running the tool indicate the number of k-mer pairs in the given sequence that had a Hamming distance of 0, 1, ..., k. The results will be formatted in 2 columns, in the following manner:  
  
Hamming Distance : Number of Pairs
  
The following is an example of the output found in results.txt after running the command `./profile test.fna results.txt 10000 32`. For more example output, please see the test folder found in the repo.
  
Hamming Distance : Number of Pairs
0 : 0
1 : 0
2 : 0
3 : 0
4 : 0
5 : 0
6 : 0
7 : 0
8 : 1
9 : 7
10 : 26
11 : 122
12 : 677
13 : 2590
14 : 9719
15 : 31298
16 : 89884
17 : 235987
18 : 555450
19 : 1170238
20 : 2200547
21 : 3672056
22 : 5406069
23 : 6981408
24 : 7828013
25 : 7542293
26 : 6158943
27 : 4177891
28 : 2285129
29 : 969730
30 : 302018
31 : 59794
32 : 5606