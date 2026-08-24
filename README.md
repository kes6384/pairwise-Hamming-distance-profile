# Pairwise Hamming Distance Profile  
  
This tool offers several methods to compute the pairwise Hamming distance profile of a given sequence. It is designed for use with FASTA files and k <= 128. The results are output to text files. The tool can be compiled and run using the command line.  
  
The tool can be used to compute the full Hamming distance profile of a sequence, or to compute a sketch of the profile for faster (but less accurate) results. An option for multithreading is also provided.  
  
All of the provided methods are contained within profile.cpp. Required open source libraries are provided in the repo.  
  
## COMPILING  
  
g++ -O3 -fPIC -fopenmp -std=c++20 -pipe -D_FILE_OFFSET_BITS=64 -o profile profile.cpp  
  
## RUNNING  

Once compiled, the tool can be used via the command line using the general format below:  
  
profile [INPUT FILE].fna [OUTPUT FILE].txt [SEQUENCE LENGTH] [k] [ROW SAMPLE RATE] [COLUMN SAMPLE RATE] [NUMBER OF THREADS]  

The sampling rate and number of threads arguments are optional.  
  
INPUT FILE - FASTA file that contains the sequence you want to analyze  
OUTPUT FILE - text file where results will be written to  
SEQUENCE LENGTH - length of the provided sequence to analyze (may be smaller than the length of the entire sequence in INPUT FILE)  
k - length of k-mers  
ROW SAMPLE RATE - rate at which to sample first k-mer; default is 1 (full profile)  
COLUMN SAMPLE RATE - rate at which to sample second k-mer; default is 1 (full profile)  
NUMBER OF THREADS - Use this option to run a multithreaded version of the tool; default is 1  
  
### EXAMPLES  
  
profile test.fna results.txt 10000 32  
This command will analyze the first 10000 characters of the sequence found in test.fna. The Hamming distance profile of 32-mers will be calculated and then stored in results.txt. Since no sampling rate is specified, the full profile will be calculated. Since no thread number is specified, the tool will run one thread.  
  
profile test.fna results.txt 10000 32 0.5 0.5  
This command will sample 32-mers at a rate of 0.25 (0.5*0.5).  
  
profile test.fna results.txt 10000 32 4  
This command will run the multithreaded version with 4 threads.  
  
profile test.fna results.txt 10000 32 0.5 0.5 4  
This command will sample 32-mers at a rate of 0.25 (0.5*0.5) with multithreading using 4 threads. 

## OUTPUT  
  
The results of running the tool indicate the number of k-mer pairs in the given sequence that had a Hamming distance of 0, 1, ..., k. The results will be formatted in 2 columns, in the following manner:  
  
Hamming Distance : Number of Pairs