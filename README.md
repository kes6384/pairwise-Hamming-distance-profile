# Pairwise Hamming Distance Profile  
  
This tool offers several methods to compute the pairwise Hamming distance profile of a given sequence. It is designed for use with FASTA files and $k \leq 128$. The tool can be compiled and run using the command line. The results are output to a user-specified text file.  
  
The tool can be used to compute the full Hamming distance profile of a sequence, or to compute a sketch of the profile for faster (but less accurate) results. Support for multithreading is also provided.  
  
All of the provided methods are contained within profile.cpp. Required open source libraries are provided in the repo. MurmurHash repo: https://github.com/aappleby/smhasher/tree/master  

## Requirements 

C++ >= 17

## INSTALLATIONCOMPILING  

```
git clone ... //
cd ./xxx
g++ -O3 -fPIC -fopenmp -std=c++20 -pipe -D_FILE_OFFSET_BITS=64 -o profile profile.cpp
```
  
## RUNNING  

Once compiled, the tool can be used via the command line using the general format below:  

```
cd ./xxx
./profile [INPUT FILE].fna [OUTPUT FILE].txt [SEQUENCE LENGTH] [k] [ROW SAMPLE RATE] [COLUMN SAMPLE RATE] [NUMBER OF THREADS]  
```

The sampling rate and number of threads arguments are optional and can be used together. Omit the sampling rate arguments to generate the full profile. To generate a sketch instead of the full profile, you must include both a row and column sampling rate.  

`-L`: length of the input sequence ... 
INPUT FILE - FASTA file that contains the sequence you want to analyze  
OUTPUT FILE - text file where results will be written to  
SEQUENCE LENGTH - length of the provided sequence to analyze (can be smaller than the length of the entire sequence in INPUT FILE)  
k - length of k-mers, must be a value between 1 and 128 (inclusive)  
ROW SAMPLE RATE - optional argument, rate at which to sample first k-mer; default is 1 (full profile)  
COLUMN SAMPLE RATE - optional argument, rate at which to sample second k-mer; default is 1 (full profile)  
NUMBER OF THREADS - optional argument, use this to run a multithreaded version of the tool; default is 1  
  
### EXAMPLES  

./test/test.fa

```
profile test.fna results.txt 10000 32
```

This command will analyze the first 10000 characters of the sequence found in test.fna. The Hamming distance profile of 32-mers will be calculated and then stored in results.txt. Since no sampling rate is specified, the full profile will be calculated. Since no thread number is specified, the tool will run with one thread.  

```
profile test.fna results.txt 10000 32 0.5 0.5
```
This command will sample 32-mer pairs at a rate of 0.25 (0.5*0.5).  

```
profile test.fna results.txt 10000 32 4
```
This command will run the multithreaded version with 4 threads.  

```
profile test.fna results.txt 10000 32 0.5 0.5 4
```

This command will sample 32-mer pairs at a rate of 0.25 (0.5*0.5) with multithreading using 4 threads. 

## OUTPUT  
  
The results of running the tool indicate the number of k-mer pairs in the given sequence that had a Hamming distance of 0, 1, ..., k. The results will be formatted in 2 columns, in the following manner:  
  
Hamming Distance : Number of Pairs

Here is an example of output (also check `./test/out.txt`)

