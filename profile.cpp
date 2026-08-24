//Algorithm for finding the Hammond Distance profile of a given sequence
//Uses popcount method
//Includes sketching and multithreaded options
//Takes FASTA files as input
//Outputs a text file representing the histogram of Hammond Distances between k-mers

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <bit>
#include <bitset>
#include <random>
#include <omp.h>
#include <vector>
#include "MurmurHash3.cpp" // https://github.com/aappleby/smhasher/tree/master
#include "parseFASTA.cpp" // Simple FASTA parser

#define POPCOUNT
#if defined(POPCOUNT)
    #define HD_FUNC hdPC
#endif

using namespace std;
//using namespace TFA;

// Masks used for popcount method
uint64_t popMask;
uint64_t popMask2;

const unsigned char seq_nt4_table[256] = { // translate ACGT to 0123
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

// parse FASTA file using parseFASTA.cpp
// Sequence is stored in seq and is of length len
int getSeq(char* seq, char* file, int len)
{
  return getSequence(file, len, seq);
}

// Returns HD of kmer1 and kmer2
// Uses bit manipulation and popcount
int hdPC(uint64_t kmer1, uint64_t kmer2 , int k)
{
  int dist = 0;
  // Mask bits to separate first and second bits of each character for each subsequence
  uint64_t subseq1_oddBit = kmer1 & popMask;
  uint64_t subseq2_oddBit = kmer2 & popMask;
  uint64_t subseq1_evenBit = kmer1 & (popMask2);
  uint64_t subseq2_evenBit = kmer2 & (popMask2);

  // XOR each masked subsequence to compare first and second bits
  // XOR result is 0 if the bits match
  uint64_t xor_firstBits = subseq1_oddBit ^ subseq2_oddBit;
  uint64_t xor_secondBits = (subseq1_evenBit ^ subseq2_evenBit) << 1;

  // NOR = 1 if both bits for the character matched = the character matched
  // ~NOR = 1 if the characters did not match
  // popcount(~nor) = number of characters that did not match
  // ~nor = or
  uint64_t orResult = xor_firstBits | xor_secondBits;
  dist = std::popcount(orResult);
  return dist;
}

// Output array of Hamming distance counts to a text file
// dists - array of Hamming distance counts
// len - length of array of Hamming distance counts (equivalent to k-mer size, since HD max is k)
// rate - sampling rate (theta1 * theta2) = 1 if not sketch
// outFile - file results are written to
void output(unsigned int *dists , int len , float rate , char* outFile)
{
  ofstream out(outFile);

  out << "Hamming Distance : Number of Pairs" << endl;
  for(int i=0; i<=len; i++)
  {
      out << ("%d" , i) << (" : ") << ("%d" , (int)((float)dists[i]/(2.0 * rate))) << endl;
  }

  out.close();
}

// Output vector of Hamming distance counts to a text file
// Used with multithreaded methods
// dists - vector of Hamming distance counts
// len - length of vector of Hamming distance counts (equivalent to k-mer size, since HD max is k)
// rate - sampling rate (theta1 * theta2) = 1 if not sketch
// outFile - file results are written to
void output_multithread(std::vector<uint64_t> &dists , int len , float rate , char* outFile)
{
  ofstream out(outFile);

  out << "Hamming Distance : Number of Pairs" << endl;
  for(int i=0; i<=len; i++)
  {
      out << ("%d" , i) << (" : ") << ("%d" , (int)((float)dists[i]/(2.0 * rate))) << endl;
  }

  out.close();
}

// No sketching, no multithreading
void regularVer(int kVal, int seqLen , int* seq , unsigned int* dists)
{
    // Iterate through k-mers and check HD for each pair
    uint64_t mask = (kVal == 32) ? ~0ULL : ((1ULL << (2 * kVal)) - 1);
    uint64_t kmer1 = 0;
    for (int i=0; i<kVal-1; i++)
    {
        int c = seq[i];
        kmer1 = ((kmer1 << 2) | (uint64_t)c) & mask;
    }
    for (int i=kVal-1; i<seqLen; i++)
    {
        int c = seq[i];
        kmer1 = ((kmer1 << 2) | (uint64_t)c) & mask;
        uint64_t kmer2 = kmer1;
        for (int j=i+1; j<seqLen; j++)
        {
            int c2 = seq[j];
            kmer2 = ((kmer2 << 2) | (uint64_t)c2) & mask;
            dists[HD_FUNC(kmer1  , kmer2 , kVal)] ++;
        }
    }
}

// Sketching
void sketch(int kVal , int seqLen , int* seq , double theta1 , double theta2 , unsigned int* dists)
{
    std::random_device rd;
    uint32_t seed1 = rd();
    uint32_t seed2 = rd();
    
    // Sample k-mers
    int numKmers = (seqLen - kVal) + 1;
    uint64_t *kmersRow = (uint64_t*)calloc(numKmers , sizeof(uint64_t));
    uint64_t *kmersCol = (uint64_t*)calloc(numKmers , sizeof(uint64_t));

    uint64_t mask = (kVal == 32) ? ~0ULL : ((1ULL << (2 * kVal)) - 1);
    uint64_t kmer1 = 0;
    for (int i=0; i<kVal-1; i++)
    {
        int c = seq[i];
        kmer1 = ((kmer1 << 2) | (uint64_t)c) & mask;
    }
    int countRow = 0;
    int countCol = 0;
    uint32_t *hashed = (uint32_t *)malloc(sizeof(uint32_t));
    for (int i=kVal-1; i<seqLen; i++)
    {
        int c = seq[i];
        kmer1 = ((kmer1 << 2) | (uint64_t)c) & mask;
        MurmurHash3_x86_32(&kmer1 , sizeof(kmer1) , seed1 , hashed);
        float hash = (float)(*hashed) / (float)(UINT32_MAX);
        if(hash < theta1)
        {
            kmersRow[countRow] = kmer1;
            countRow ++;
        }
        // Different seeds for row and column for independence
        MurmurHash3_x86_32(&kmer1 , sizeof(kmer1) , seed2 , hashed);
        hash = (float)(*hashed) / (float)(UINT32_MAX);
        if(hash < theta2)
        {
            kmersCol[countCol] = kmer1;
            countCol ++;
        }
        }
        free(hashed);

        // Calculate HD for sampled k-mers
        for(int i=0; i<countRow; i++)
        {
        uint64_t kmer1 = kmersRow[i];
        for(int j=0; j<countCol; j++)
        {
            uint64_t kmer2 = kmersCol[j];
            dists[HD_FUNC(kmer1 , kmer2 , kVal)] ++;
        }
    }
    
    free(kmersRow);
    free(kmersCol);
}

// Multithreading with no sketching
void multithread(int kVal , int seqLen , int* seq , int numThreads , char* outFile)
{
    int numKmers = (seqLen - kVal) + 1;
    uint64_t *kmers = (uint64_t*)calloc(numKmers , sizeof(uint64_t));

    // Iterate through k-mers and check HD for each pair
    uint64_t mask = (kVal == 32) ? ~0ULL : ((1ULL << (2 * kVal)) - 1);
    uint64_t kmer1 = 0;
    for (int i=0; i<kVal-1; i++)
    {
        int c = seq[i];
        kmer1 = ((kmer1 << 2) | (uint64_t)c) & mask;
    }
    int count = 0;
    for (int i=kVal-1; i<seqLen; i++)
    {
        int c = seq[i];
        kmer1 = ((kmer1 << 2) | (uint64_t)c) & mask;
        kmers[count] = kmer1;
        count ++;
    }

    std::vector<vector<uint64_t>> local_dists(numThreads , vector<uint64_t>(kVal + 1 , 0));

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        auto& hist = local_dists[tid];
        #pragma omp for schedule(dynamic , 1)
        for (int i = 0; i < numKmers; i++) {
            const uint64_t ki = kmers[i];
            for (int j = i + 1; j < numKmers; j++) {
                int hd = hdPC(ki, kmers[j], kVal);
                hist[hd]++;
            }
        }
    }

    std::vector<uint64_t> dists(kVal + 1 , 0);
    for(int t=0; t<numThreads; t++)
    {
        for(int h=0; h<=kVal; h++)
        {
        dists[h] += local_dists[t][h];
        }
    }
    
    free(kmers);
    output_multithread(dists , kVal , 0.5 , outFile);
}

void sketch_multithread(int kVal , int seqLen , int* seq , double theta1 , double theta2 , int numThreads , char* outFile)
{
    // Sample k-mers
    int numKmers = (seqLen - kVal) + 1;
    uint64_t *kmersRow = (uint64_t*)calloc(numKmers , sizeof(uint64_t));
    uint64_t *kmersCol = (uint64_t*)calloc(numKmers , sizeof(uint64_t));

    uint64_t mask = (kVal == 32) ? ~0ULL : ((1ULL << (2 * kVal)) - 1);
    uint64_t kmer1 = 0;
    for (int i=0; i<kVal-1; i++)
    {
        int c = seq[i];
        kmer1 = ((kmer1 << 2) | (uint64_t)c) & mask;
    }
    int countRow = 0;
    int countCol = 0;
    uint32_t *hashed = (uint32_t *)malloc(sizeof(uint32_t));
    for (int i=kVal-1; i<seqLen; i++)
    {
        int c = seq[i];
        kmer1 = ((kmer1 << 2) | (uint64_t)c) & mask;
        MurmurHash3_x86_32(&kmer1 , sizeof(kmer1) , 1 , hashed);
        float hash = (float)(*hashed) / (float)(UINT32_MAX);
        if(hash < theta1)
        {
            kmersRow[countRow] = kmer1;
            countRow ++;
        }
        MurmurHash3_x86_32(&kmer1 , sizeof(kmer1) , 2 , hashed);
        hash = (float)(*hashed) / (float)(UINT32_MAX);
        if(hash < theta2)
        {
            kmersCol[countCol] = kmer1;
            countCol ++;
        }
    }
    free(hashed);

    // Calculate HD for sampled k-mers
    std::vector<vector<uint64_t>> local_dists(numThreads , vector<uint64_t>(kVal + 1 , 0));

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        auto& hist = local_dists[tid];
        #pragma omp for schedule(dynamic , 1)
        for (int i = 0; i < countRow; i++) {
            const uint64_t ki = kmersRow[i];
            for (int j = 0; j < countCol; j++) {
                int hd = hdPC(ki, kmersCol[j], kVal);
                hist[hd]++;
            }
        }
    }

    std::vector<uint64_t> dists(kVal + 1 , 0);
    for(int t=0; t<numThreads; t++)
    {
        for(int h=0; h<=kVal; h++)
        {
            dists[h] += local_dists[t][h];
        }
    }

    free(kmersRow);
    free(kmersCol);
    output_multithread(dists , kVal , (theta1*theta2) , outFile);
}

// Command line args: input file name, output file name, seqLen, k, sampling rate 1, sampling rate 2, number of threads
// Sampling rates and number of threads are optional
// Default is full profile and one thread (no multithreading)
int main(int argc, char* argv[]) {
    char* inptFile = argv[1];
    char* outptFile = argv[2];
    int seqLen = atoi(argv[3]); // Sequence length
    int kVal = atoi(argv[4]); // k-mer length
    double theta1 = 1;
    double theta2 = 1;
    int numThreads = 1;

    popMask = (2.0)*((pow((long double)4 , (long double)kVal) - 1)/3.0); // Keep odd bits
    popMask2 = popMask >> 1; // Keep even bits

    // Parse FASTA file to get sequence
    char *charSeq = (char *)malloc(seqLen + 1);
    int retrievedLen = getSeq(charSeq , inptFile , seqLen);

    // 2-bit encoding
    int* seq = (int *)calloc(retrievedLen , sizeof(int));
    // Store as array of integers, each int is one encoded character
    for (int i=0; i<retrievedLen; i++)
    {
        seq[i] = seq_nt4_table[(uint8_t)charSeq[i]]; 
    }

    // No sketch, multithreading
    if(argc == 6)
    {
        int numThreads = atoi(argv[5]);
        omp_set_num_threads(numThreads);
        multithread(kVal , retrievedLen , seq , numThreads , outptFile);
    }
    // Sketch, no multithreading
    else if(argc == 7)
    {
        double theta1 = atof(argv[5]); // row sampling rate
        double theta2 = atof(argv[6]); // column sampling rate
        // Tracks how many pairs had a Hamming distance of i, where i is an index of the array
        unsigned int *dists = (unsigned int *)calloc(kVal+1 , sizeof(int));
        sketch(kVal , retrievedLen , seq , theta1 , theta2 , dists);
        output(dists , kVal , (theta1*theta2) , outptFile);
        free(dists);
    }
    // Sketch with multithreading
    else if(argc == 8)
    {
        double theta1 = atof(argv[5]); // row sampling rate
        double theta2 = atof(argv[6]); // column sampling rate
        int numThreads = atoi(argv[7]);
        omp_set_num_threads(numThreads);
        sketch_multithread(kVal , retrievedLen , seq , theta1 , theta2 , numThreads , outptFile);
    }
    // No sketch, no multithreading
    else
    {
        // Tracks how many pairs had a Hamming distance of i, where i is an index of the array
        unsigned int *dists = (unsigned int *)calloc(kVal+1 , sizeof(int));
        regularVer(kVal , retrievedLen , seq , dists);
        output(dists , kVal , 0.5 , outptFile);
        free(dists);
    }

    free(seq);
    free(charSeq);

    return 0;
}