# PAIRWISE HAMMING DISTANCE PROFILE  
  
This tool offers several methods to compute the pairwise Hamming distance profile of a given sequence. It is designed for use with FASTA files and $k \leq 128$. The tool can be compiled and run using the command line. The results are output to a user-specified text file.  
  
The tool can be used to compute the full Hamming distance profile of a sequence, or to compute a sketch of the profile for faster (but less accurate) results. Support for multithreading is also provided.  
  
All of the provided methods are contained within profile.cpp. Required open source libraries are provided in the repo. MurmurHash repo: https://github.com/aappleby/smhasher/tree/master  

## REQUIREMENTS  
  
C++ $\geq 20$  
Windows OS $\geq 11$  
  
  
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
./profile -i [INPUT FILE] -o [OUTPUT FILE].txt -l [SEQUENCE LENGTH] -k [k] -s [ROW SAMPLE RATE] [COLUMN SAMPLE RATE] -t [NUMBER OF THREADS]  
```

The arguments can be specified in any order using the options outlined below. The sampling rate and number of threads arguments are optional and can be used together. All other arguments are required. Omit the sampling rate arguments to generate the full profile. To generate a sketch instead of the full profile, you must include both a row and column sampling rate.  

`-i` : FASTA file that contains the sequence you want to analyze  
`-o` : text file where results will be written to  
`-l` : length of the input sequence to analyze (can be smaller than the length of the entire sequence in INPUT FILE)  
`-k` : length of k-mers, must be $k 1 \leq k \leq 128$ 
`-s` : generate a sketch instead of the full profile. Requires both row and column sampling rates  
&nbsp;&nbsp;&nbsp;&nbsp;- ROW SAMPLE RATE - rate at which to sample first k-mer in each pair  
&nbsp;&nbsp;&nbsp;&nbsp;- COLUMN SAMPLE RATE - rate at which to sample second k-mer in each pair  
`-t` : run multithreaded version of the tool. Otherwise, the tool will run with one thread  
  
### EXAMPLES  
The FASTA file used for the following examples and the results of running the following commands can be found in the [test folder](/test/). 

To use one thread to calculate the full Hamming distance profile of 32-mers for the first 10,000 characters of the sequence from the file test.fna in the folder test and store the results in results.txt, any of the following three commands will work:
```
./profile -i test/test.fna -o results.txt -l 10000 -k 32
```
```
./profile -i test/test.fna -o results.txt -l 10000 -k 32 -s 1.0 1.0
```
```
./profile -i test/test.fna -o results.txt -l 10000 -k 32 -t 1
```
```
./profile -i test/test.fna -o results.txt -l 10000 -k 32 -s 1.0 1.0 -t 1
```
  
To use one thread to sample 32-mer pairs at a rate of 0.25, any of the following commands will work:
```
./profile -i test/test.fna -o results.txt -l 10000 -k 32 -s 0.5 0.5
```
```
./profile -i test/test.fna -o results.txt -l 10000 -k 32 -s 0.5 0.5 -t 1
```
  
To calculate the full profile using 4 threads, any of the following commands will work:
```
./profile -i test/test.fna -o results.txt -l 10000 -k 32 -t 4
```
```
./profile -i test/test.fna -o results.txt -l 10000 -k 32 -s 1.0 1.0 -t 4
```
  
To use 4 threads to sample 32-mer pairs at a rate of 0.25, use the following command:
```
./profile -i test/test.fna -o results.txt -l 10000 -k 32 -s 0.5 0.5 -t 4
```
  
## OUTPUT  
  
The results of running the tool indicate the number of k-mer pairs in the given sequence that had a Hamming distance of 0, 1, ..., k. The results will be formatted in 2 columns, in the following manner:  
  
Hamming Distance : Number of Pairs
  
For example outputs, please see the [test folder](/test/) found in the repo. All example outputs were generated using k = 32 and sequence length = 10,000. [results_fullProfile.txt](/test/results_fullProfile.txt) was generated using the default method. [results_sampleRate0.5.txt](/test/results_sampleRate0.5.txt) was generated using a row sampling rate of 0.5 and a column sampling rate of 0.5. Changing the number of threads has no effect on the output generated.  