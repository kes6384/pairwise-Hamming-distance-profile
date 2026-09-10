//Algorithm for finding the Hammond Distance profile of a given sequence
//Uses popcount method
//Includes sketching and multithreaded options
//Takes FASTA files as input and supports k<=128
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

// Struct for storing k-mers with k <= 128
struct uint128 {
    // Each k-mer is 2-bit encoded, so need 256 bits for 128 characters
    uint64_t fullKmer[4] = {0,0,0,0};
};

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
int hdPC(uint128 kmer1, uint128 kmer2 , int k)
{
    //cout << ("%d" , kmer1.fullKmer[3]) << endl;
    //cout << ("%d" , kmer2.fullKmer[3]) << endl;
  int numLoops = (k / 33) + 1;
  //int numLoops = 4;
  int dist = 0;
  // Mask bits to separate first and second bits of each character for each subsequence
  uint128 subseq1_oddBit;
  uint128 subseq1_evenBit;
  uint128 subseq2_oddBit;
  uint128 subseq2_evenBit;

  for(int i=0; i < numLoops; i++)
  {
    subseq1_oddBit.fullKmer[3-i] = kmer1.fullKmer[3-i] & popMask;
  }

  for(int i=0; i < numLoops; i++)
  {
    subseq2_oddBit.fullKmer[3-i] = kmer2.fullKmer[3-i] & popMask;
  }

  for(int i=0; i < numLoops; i++)
  {
    subseq1_evenBit.fullKmer[3-i] = kmer1.fullKmer[3-i] & popMask2;
  }

  for(int i=0; i < numLoops; i++)
  {
    subseq2_evenBit.fullKmer[3-i] = kmer2.fullKmer[3-i] & popMask2;
  }

  // XOR each masked subsequence to compare first and second bits
  // XOR result is 0 if the bits match
  uint128 xor_firstBits;
  uint128 xor_secondBits;

  for(int i=0; i < numLoops; i++)
  {
    xor_firstBits.fullKmer[3-i] = subseq1_oddBit.fullKmer[3-i] ^ subseq2_oddBit.fullKmer[3-i];
  }

  for(int i=0; i < numLoops; i++)
  {
    xor_secondBits.fullKmer[3-i] = (subseq1_evenBit.fullKmer[3-i] ^ subseq2_evenBit.fullKmer[3-i]) << 1;
  }

  // NOR = 1 if both bits for the character matched = the character matched
  // ~NOR = 1 if the characters did not match
  // popcount(~nor) = number of characters that did not match
  // ~nor = or
  uint128 orResult;
  for(int i=0; i < numLoops; i++)
  {
    orResult.fullKmer[3-i] = xor_firstBits.fullKmer[3-i] | xor_secondBits.fullKmer[3-i];
  }

  for(int i=0; i < numLoops; i++)
  {
    dist += std::popcount(orResult.fullKmer[3-i]);
  }

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
    int numLoops = (kVal / 33) + 1;
    //int numLoops = 4;
    // Iterate through k-mers and check HD for each pair
    uint128 mask;
    for(int i=0; i < numLoops; i++)
    {
        mask.fullKmer[i] = (kVal >= ((i+1)*32)) ? ~0ULL : ((1ULL << (2 * (kVal-(32*i)))) - 1);
    }
    uint128 kmer1;
    for (int i=0; i<kVal-1; i++)
    {
        int c = seq[i];
        for(int i=0; i < 4; i++)
        {
            if(i == 3)
            {
                kmer1.fullKmer[3] = ((kmer1.fullKmer[3] << 2) | (uint64_t)c) & mask.fullKmer[0];
            }
            else if(3-i < numLoops)
            {
                kmer1.fullKmer[i] = ((kmer1.fullKmer[i] << 2) | (kmer1.fullKmer[i+1] >> 62)) & mask.fullKmer[3-i];
            }
        }
    }
    //cout << ("%d" , kmer1.fullKmer[3]) << endl;
    for (int i=kVal-1; i<seqLen; i++)
    {
        int c = seq[i];
        for(int i=0; i < 4; i++)
        {
            if(i == 3)
            {
                kmer1.fullKmer[3] = ((kmer1.fullKmer[3] << 2) | (uint64_t)c) & mask.fullKmer[0];
            }
            else if(3-i < numLoops)
            {
                kmer1.fullKmer[i] = ((kmer1.fullKmer[i] << 2) | (kmer1.fullKmer[i+1] >> 62)) & mask.fullKmer[3-i];
            }
        }
        uint128 kmer2 = kmer1;
        for (int j=i+1; j<seqLen; j++)
        {
            int c2 = seq[j];
            for(int i=0; i < 4; i++)
            {
                if(i == 3)
                {
                    kmer2.fullKmer[3] = ((kmer2.fullKmer[3] << 2) | (uint64_t)c2) & mask.fullKmer[0];
                }
                else if(3-i < numLoops)
                {
                    kmer2.fullKmer[i] = ((kmer2.fullKmer[i] << 2) | (kmer2.fullKmer[i+1] >> 62)) & mask.fullKmer[3-i];
                }
            }
            dists[HD_FUNC(kmer1  , kmer2 , kVal)] ++;
        }
    }
}

// Sketching
void sketch(int kVal , int seqLen , int* seq , double theta1 , double theta2 , unsigned int* dists)
{
    //std::random_device rd;
    // arbitrary seeds for hashing
    uint32_t seed1 = 42;
    uint32_t seed2 = 24;
    
    int numLoops = (kVal / 33) + 1;
    // Sample k-mers
    int numKmers = (seqLen - kVal) + 1;
    uint128 *kmersRow = (uint128*)calloc(numKmers , sizeof(uint128));
    uint128 *kmersCol = (uint128*)calloc(numKmers , sizeof(uint128));

    uint128 mask;
    for(int i=0; i < numLoops; i++)
    {
        mask.fullKmer[i] = (kVal >= ((i+1)*32)) ? ~0ULL : ((1ULL << (2 * (kVal-(32*i)))) - 1);
    }
    uint128 kmer1;
    for (int i=0; i<kVal-1; i++)
    {
        int c = seq[i];
        for(int i=0; i < 4; i++)
        {
            if(i == 3)
            {
                kmer1.fullKmer[3] = ((kmer1.fullKmer[3] << 2) | (uint64_t)c) & mask.fullKmer[0];
            }
            else if(3-i < numLoops)
            {
                kmer1.fullKmer[i] = ((kmer1.fullKmer[i] << 2) | (kmer1.fullKmer[i+1] >> 62)) & mask.fullKmer[3-i];
            }
        }
    }
    int countRow = 0;
    int countCol = 0;
    uint32_t *hashed = (uint32_t *)malloc(sizeof(uint32_t));
    for (int i=kVal-1; i<seqLen; i++)
    {
        int c = seq[i];
        for(int i=0; i < 4; i++)
        {
            if(i == 3)
            {
                kmer1.fullKmer[3] = ((kmer1.fullKmer[3] << 2) | (uint64_t)c) & mask.fullKmer[0];
            }
            else if(3-i < numLoops)
            {
                kmer1.fullKmer[i] = ((kmer1.fullKmer[i] << 2) | (kmer1.fullKmer[i+1] >> 62)) & mask.fullKmer[3-i];
            }
        }
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
        uint128 kmer1 = kmersRow[i];
        for(int j=0; j<countCol; j++)
        {
            uint128 kmer2 = kmersCol[j];
            dists[HD_FUNC(kmer1 , kmer2 , kVal)] ++;
        }
    }
    
    free(kmersRow);
    free(kmersCol);
    // Account for kmers being in both the row and the column
}

// Multithreading with no sketching
void multithread(int kVal , int seqLen , int* seq , int numThreads , char* outFile)
{
    int numKmers = (seqLen - kVal) + 1;
    uint128 *kmers = (uint128*)calloc(numKmers , sizeof(uint128));

    int numLoops = (kVal / 33) + 1;
    // Iterate through k-mers and check HD for each pair
    uint128 mask;
    for(int i=0; i < numLoops; i++)
    {
        mask.fullKmer[i] = (kVal >= ((i+1)*32)) ? ~0ULL : ((1ULL << (2 * (kVal-(32*i)))) - 1);
    }
    uint128 kmer1;
    for (int i=0; i<kVal-1; i++)
    {
        int c = seq[i];
        for(int i=0; i < 4; i++)
        {
            if(i == 3)
            {
                kmer1.fullKmer[3] = ((kmer1.fullKmer[3] << 2) | (uint64_t)c) & mask.fullKmer[0];
            }
            else if(3-i < numLoops)
            {
                kmer1.fullKmer[i] = ((kmer1.fullKmer[i] << 2) | (kmer1.fullKmer[i+1] >> 62)) & mask.fullKmer[3-i];
            }
        }
    }
    int count = 0;
    for (int i=kVal-1; i<seqLen; i++)
    {
        int c = seq[i];
        for(int i=0; i < 4; i++)
        {
            if(i == 3)
            {
                kmer1.fullKmer[3] = ((kmer1.fullKmer[3] << 2) | (uint64_t)c) & mask.fullKmer[0];
            }
            else if(3-i < numLoops)
            {
                kmer1.fullKmer[i] = ((kmer1.fullKmer[i] << 2) | (kmer1.fullKmer[i+1] >> 62)) & mask.fullKmer[3-i];
            }
        }
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
            const uint128 ki = kmers[i];
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
    // arbitrary seeds for hashing
    uint32_t seed1 = 42;
    uint32_t seed2 = 24;

    int numLoops = (kVal / 33) + 1;

    // Sample k-mers
    int numKmers = (seqLen - kVal) + 1;
    uint128 *kmersRow = (uint128*)calloc(numKmers , sizeof(uint128));
    uint128 *kmersCol = (uint128*)calloc(numKmers , sizeof(uint128));

    uint128 mask;
    for(int i=0; i < numLoops; i++)
    {
        mask.fullKmer[i] = (kVal >= ((i+1)*32)) ? ~0ULL : ((1ULL << (2 * (kVal-(32*i)))) - 1);
    }
    uint128 kmer1;
    for (int i=0; i<kVal-1; i++)
    {
        int c = seq[i];
        for(int i=0; i < 4; i++)
        {
            if(i == 3)
            {
                kmer1.fullKmer[3] = ((kmer1.fullKmer[3] << 2) | (uint64_t)c) & mask.fullKmer[0];
            }
            else if(3-i < numLoops)
            {
                kmer1.fullKmer[i] = ((kmer1.fullKmer[i] << 2) | (kmer1.fullKmer[i+1] >> 62)) & mask.fullKmer[3-i];
            }
        }
    }
    int countRow = 0;
    int countCol = 0;
    uint32_t *hashed = (uint32_t *)malloc(sizeof(uint32_t));
    for (int i=kVal-1; i<seqLen; i++)
    {
        int c = seq[i];
        for(int i=0; i < 4; i++)
        {
            if(i == 3)
            {
                kmer1.fullKmer[3] = ((kmer1.fullKmer[3] << 2) | (uint64_t)c) & mask.fullKmer[0];
            }
            else if(3-i < numLoops)
            {
                kmer1.fullKmer[i] = ((kmer1.fullKmer[i] << 2) | (kmer1.fullKmer[i+1] >> 62)) & mask.fullKmer[3-i];
            }
        }
        MurmurHash3_x86_32(&kmer1 , sizeof(kmer1) , seed1 , hashed);
        float hash = (float)(*hashed) / (float)(UINT32_MAX);
        if(hash < theta1)
        {
            kmersRow[countRow] = kmer1;
            countRow ++;
        }
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
    std::vector<vector<uint64_t>> local_dists(numThreads , vector<uint64_t>(kVal + 1 , 0));

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        auto& hist = local_dists[tid];
        #pragma omp for schedule(dynamic , 1)
        for (int i = 0; i < countRow; i++) {
            const uint128 ki = kmersRow[i];
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
    char* inptFile;
    char* outptFile;
    int seqLen = 0;
    int kVal = 0;
    double theta1 = 1;
    double theta2 = 1;
    int numThreads = 1;

    // parse arguments
    for (int i = 1; i < argc; i++)
    {
        string arg = argv[i];

        if(arg == "-i") //input file
        {
            if(argc < (i + 1))
            {
                cout << "ARGUMENT ERROR" << endl;
                return 1;
            }
            inptFile = argv[++i];
        }
        else if(arg == "-o") //output file
        {
            if(argc < (i + 1))
            {
                cout << "ARGUMENT ERROR" << endl;
                return 1;
            }
            outptFile = argv[++i];
        }
        else if(arg == "-l") //sequence length
        {
            if(argc < (i + 1))
            {
                cout << "ARGUMENT ERROR" << endl;
                return 1;
            }
            seqLen = atoi(argv[++i]);
        }
        else if(arg == "-k") //k-mer length
        {
            if(argc < (i + 1))
            {
                cout << "ARGUMENT ERROR" << endl;
                return 1;
            }
            kVal = atoi(argv[++i]);
        }
        else if(arg == "-t") //multithreading
        {
            if(argc < (i + 1))
            {
                cout << "ARGUMENT ERROR" << endl;
                return 1;
            }
            numThreads = atoi(argv[++i]);
        }
        else if(arg == "-s") //sketch
        {
            if(argc < (i + 2))
            {
                cout << "ARGUMENT ERROR" << endl;
                return 1;
            }
            theta1 = atof(argv[++i]);
            theta2 = atof(argv[++i]);
        }
        else
        {
            cout << "ARGUMENT ERROR" << endl;
            return 1;
        }
    }
    // Check for input errors
    if(kVal > 128 || kVal < 1 || seqLen < 1 || seqLen < kVal || numThreads < 1 || theta1 <= 0 || theta2 <= 0 || theta1 > 1 || theta2 > 1 || inptFile == "" || outptFile == "")
    {
        cout << "ARGUMENT ERROR" << endl;
        return 1;
    }

    popMask = (2.0)*((pow((long double)4 , (long double)(min(kVal , 32))) - 1)/3.0); // Keep odd bits
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

    // Run correct method

    // multithreaded
    if(numThreads > 1)
    {
        omp_set_num_threads(numThreads);
        // sketching
        if(theta1*theta2 != 1)
        {
            sketch_multithread(kVal , retrievedLen , seq , theta1 , theta2 , numThreads , outptFile);
        }
        else
        {
            multithread(kVal , retrievedLen , seq , numThreads , outptFile);
        }
    }
    else if(theta1*theta2 != 1)
    {
        // Tracks how many pairs had a Hamming distance of i, where i is an index of the array
        unsigned int *dists = (unsigned int *)calloc(kVal+1 , sizeof(int));
        sketch(kVal , retrievedLen , seq , theta1 , theta2 , dists);
        output(dists , kVal , (theta1*theta2) , outptFile);
        free(dists);
    }
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