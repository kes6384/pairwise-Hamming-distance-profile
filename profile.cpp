//Algorithm for finding the Hammond Distance profile of a given sequence
//Uses popcount method
//Includes sketching and multithreaded options
//Takes FASTA files as input and supports k<=128
//Outputs a text file representing the histogram of Hammond Distances between k-mers
//Sketch versions make use of MurmurHash3 by Austin Appleby https://github.com/aappleby/smhasher/tree/master

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

#define HD_FUNC hdPC

using namespace std;

// Struct for storing k-mers with k <= 128
struct uint128 {
    // Each k-mer is 2-bit encoded, so need 256 bits for 128 characters
    uint64_t fullKmer[4] = {0,0,0,0};
};
// How many elements of fullKmer are actually needed
int maxNum;

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
// seq - extracted sequence will be stored here
// len - length of sequence to extract
int getSeq(char* seq, char* file, unsigned int len)
{
  return getSequence(file, len, seq);
}

// Separates out and stores k-mers found in sequence
// kmers - array to store kmers in
// kVal - k
// seqLen - length of sequence
// seq - sequence to analyze
void storeKmers(uint128 *kmers , int kVal , unsigned int seqLen , int* seq)
{
    uint128 mask;
    for(int i=0; i < maxNum+1; i++)
    {
        mask.fullKmer[i] = (kVal >= ((i+1)*32)) ? ~0ULL : ((1ULL << (2 * (kVal-(32*i)))) - 1);
    }
    uint128 kmer1;
    for (int i=0; i<kVal-1; i++)
    {
        int c = seq[i];
        for(int j=0; j < maxNum; j++)
        {
            kmer1.fullKmer[j] = ((kmer1.fullKmer[j] << 2) | (kmer1.fullKmer[j+1] >> 62)) & mask.fullKmer[maxNum-j];
        }
        kmer1.fullKmer[maxNum] = ((kmer1.fullKmer[maxNum] << 2) | (uint64_t)c) & mask.fullKmer[0];
    }
    int count = 0;
    for (unsigned int i=kVal-1; i<seqLen; i++)
    {
        int c = seq[i];
        for(int j=0; j < maxNum; j++)
        {
            kmer1.fullKmer[j] = ((kmer1.fullKmer[j] << 2) | (kmer1.fullKmer[j+1] >> 62)) & mask.fullKmer[maxNum-j];
        }
        kmer1.fullKmer[maxNum] = ((kmer1.fullKmer[maxNum] << 2) | (uint64_t)c) & mask.fullKmer[0];
        kmers[count] = kmer1;
        count ++;
    }
}

// Separates out and stores a sketch of k-mers in sequence
// kmersRow - vector to store row k-mers in
// kmersCol - vector to store column k-mers in
// countRow - address to store number of k-mers in row
// countCol - address to store number of k-mers in column
// kVal - k
// seqLen - length of sequence
// seq - sequence to analyze
// theta1 - row sampling rate
// theta2 - column sampling rate
void storeKmersSketch(std::vector<uint128> &kmersRow , std::vector<uint128> &kmersCol , unsigned int* countRow , unsigned int* countCol , int kVal , unsigned int seqLen , int* seq , double theta1 , double theta2)
{
    // arbitrary seeds for hashing
    uint32_t seed1 = 42;
    uint32_t seed2 = 24;
    
    uint128 mask;
    for(int i=0; i < maxNum+1; i++)
    {
        mask.fullKmer[i] = (kVal >= ((i+1)*32)) ? ~0ULL : ((1ULL << (2 * (kVal-(32*i)))) - 1);
    }
    uint128 kmer1;
    for (int i=0; i<kVal-1; i++)
    {
        int c = seq[i];
        for(int j=0; j < maxNum; j++)
        {
            kmer1.fullKmer[j] = ((kmer1.fullKmer[j] << 2) | (kmer1.fullKmer[j+1] >> 62)) & mask.fullKmer[maxNum-j];
        }
        kmer1.fullKmer[maxNum] = ((kmer1.fullKmer[maxNum] << 2) | (uint64_t)c) & mask.fullKmer[0];
    }
    *countRow = 0;
    *countCol = 0;
    uint32_t *hashed = (uint32_t *)malloc(sizeof(uint32_t));
    for (unsigned int i=kVal-1; i<seqLen; i++)
    {
        int c = seq[i];
        for(int j=0; j < maxNum; j++)
        {
            kmer1.fullKmer[j] = ((kmer1.fullKmer[j] << 2) | (kmer1.fullKmer[j+1] >> 62)) & mask.fullKmer[maxNum-j];
        }
        kmer1.fullKmer[maxNum] = ((kmer1.fullKmer[maxNum] << 2) | (uint64_t)c) & mask.fullKmer[0];
        // Decide if kmer should be sampled
        MurmurHash3_x86_32(&kmer1 , sizeof(kmer1) , seed1 , hashed);
        float hash = (float)(*hashed) / (float)(UINT32_MAX);
        if(hash < theta1)
        {
            if(*countRow < kmersRow.size())
            {
                kmersRow[*countRow] = kmer1;
            }
            else
            {
                kmersRow.push_back(kmer1);
            }
            (*countRow) ++;
        }
        // Different seeds for row and column for independence
        MurmurHash3_x86_32(&kmer1 , sizeof(kmer1) , seed2 , hashed);
        hash = (float)(*hashed) / (float)(UINT32_MAX);
        if(hash < theta2)
        {
            if(*countCol < kmersCol.size())
            {
                kmersCol[*countCol] = kmer1;
            }
            else
            {
                kmersCol.push_back(kmer1);
            }
            (*countCol) ++;
        }
    }
    free(hashed);
}

// Returns HD of kmer1 and kmer2
// Uses bit manipulation and popcount
int hdPC(uint128 kmer1, uint128 kmer2 , int k)
{
  int dist = 0;

  uint128 subseq1_oddBit;
  uint128 subseq1_evenBit;
  uint128 subseq2_oddBit;
  uint128 subseq2_evenBit;

  uint128 xor_firstBits;
  uint128 xor_secondBits;

  uint128 orResult;

  for(int i=0; i < maxNum+1; i++)
  {
    // Mask bits to separate first and second bits of each character for each subsequence
    subseq1_oddBit.fullKmer[i] = kmer1.fullKmer[i] & popMask;
    subseq2_oddBit.fullKmer[i] = kmer2.fullKmer[i] & popMask;
    subseq1_evenBit.fullKmer[i] = kmer1.fullKmer[i] & popMask2;
    subseq2_evenBit.fullKmer[i] = kmer2.fullKmer[i] & popMask2;

    // XOR each masked subsequence to compare first and second bits
    // XOR result is 0 if the bits match
    xor_firstBits.fullKmer[i] = subseq1_oddBit.fullKmer[i] ^ subseq2_oddBit.fullKmer[i];
    xor_secondBits.fullKmer[i] = (subseq1_evenBit.fullKmer[i] ^ subseq2_evenBit.fullKmer[i]) << 1;

    // NOR = 1 if both bits for the character matched = the character matched
    // ~NOR = 1 if the characters did not match
    // popcount(~nor) = number of characters that did not match
    // ~nor = or
    orResult.fullKmer[i] = xor_firstBits.fullKmer[i] | xor_secondBits.fullKmer[i];

    dist += std::popcount(orResult.fullKmer[i]);
  }

  return dist;
}

// Output array of Hamming distance counts to a text file
// dists - array of Hamming distance counts
// len - length of array of Hamming distance counts (equivalent to k-mer size, since max HD is k)
// rate - sampling rate (theta1 * theta2); equals 1 if not sketch
// outFile - file results are written to
void output(unsigned int *dists , int len , float rate , char* outFile)
{
  ofstream out(outFile);

  out << "Hamming Distance : Number of Pairs" << endl;
  for(int i=1; i<=len; i++)
  {
      out << ("%d" , i) << (" : ") << ("%d" , (int)((float)dists[i]/(2.0 * rate))) << endl;
  }

  out.close();
}

// Compute full HD profile of sequence with no sketching, no multithreading
// kVal - k
// seqLen - length of sequence
// seq - sequence to analyze
// dists - array to store HDs in
void regularVer(int kVal, unsigned int seqLen , int* seq , unsigned int* dists)
{
    unsigned int numKmers = (seqLen - kVal) + 1;
    uint128 *kmers = (uint128*)calloc(numKmers , sizeof(uint128));

    // Iterate through and store k-mers
    storeKmers(kmers , kVal , seqLen , seq);

    // Calculate HD for each pair
    for (unsigned int i=0; i<numKmers; i++)
    {
        for(unsigned int j=i+1; j<numKmers; j++)
        {
            dists[HD_FUNC(kmers[i]  , kmers[j] , kVal)] +=2;
        }
    }
    free(kmers);
}

// Compute sketch of HD profile of sequence with no multithreading
// kVal - k
// seqLen - length of sequence
// seq - sequence to analyze
// theta1 - row sampling rate
// theta2 - column sampling rate
// dists - array to store HDs in
void sketch(int kVal , unsigned int seqLen , int* seq , double theta1 , double theta2 , unsigned int* dists)
{   
    unsigned int numKmers = (seqLen - kVal) + 1;
    uint128 empty;
    // Initially assume needed space based on sampling rate
    std::vector<uint128> kmersRow((unsigned int)(numKmers*(theta1*theta2)) , empty);
    std::vector<uint128> kmersCol((unsigned int)(numKmers*(theta1*theta2)) , empty);

    unsigned int countRow;
    unsigned int countCol;

    // Sample k-mers
    storeKmersSketch(kmersRow , kmersCol , &countRow , &countCol , kVal , seqLen , seq , theta1 , theta2);

    // Calculate HD for sampled k-mers
    for(unsigned int i=0; i<countRow; i++)
    {
        uint128 kmer1 = kmersRow[i];
        for(unsigned int j=0; j<countCol; j++)
        {
            uint128 kmer2 = kmersCol[j];
            dists[HD_FUNC(kmer1 , kmer2 , kVal)] ++;
        }
    }
}

// Compute full HD profile of sequence with multithreading and no sketching
// kVal - k
// seqLen - length of sequence
// seq - sequence to analyze
// numThreads - number of threads to run
// dists - array to store HDs in
void multithread(int kVal , unsigned int seqLen , int* seq , int numThreads , unsigned int* dists)
{
    unsigned int numKmers = (seqLen - kVal) + 1;
    uint128 *kmers = (uint128*)calloc(numKmers , sizeof(uint128));

    // Iterate through and store k-mers
    storeKmers(kmers , kVal , seqLen , seq);

    // Calculate HD of each k-mer pair
    std::vector<vector<uint64_t>> local_dists(numThreads , vector<uint64_t>(kVal + 1 , 0));

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        auto& hist = local_dists[tid];
        #pragma omp for schedule(dynamic , 1)
        for (unsigned int i = 0; i < numKmers; i++) {
            const uint128 ki = kmers[i];
            for (unsigned int j = i + 1; j < numKmers; j++) {
                int hd = hdPC(ki, kmers[j], kVal);
                hist[hd]+=2;
            }
        }
    }

    for(int t=0; t<numThreads; t++)
    {
        for(int h=0; h<=kVal; h++)
        {
            dists[h] += local_dists[t][h];
        }
    }
    
    free(kmers);
}

// Compute sketch of HD profile of sequence with multithreading
// kVal - k
// seqLen - length of sequence
// seq - sequence to analyze
// theta1 - row sampling rate
// theta2 - column sampling rate
// numThreads - number of threads to run
// dists - array to store HDs in
void sketch_multithread(int kVal , unsigned int seqLen , int* seq , double theta1 , double theta2 , int numThreads , unsigned int* dists)
{
    unsigned int numKmers = (seqLen - kVal) + 1;
    uint128 empty;
    // Initially assume needed space based on sampling rate
    std::vector<uint128> kmersRow((unsigned int)(numKmers*(theta1*theta2)) , empty);
    std::vector<uint128> kmersCol((unsigned int)(numKmers*(theta1*theta2)) , empty);

    unsigned int countRow;
    unsigned int countCol;

    // Sample k-mers
    storeKmersSketch(kmersRow , kmersCol , &countRow , &countCol , kVal , seqLen , seq , theta1 , theta2);

    // Calculate HD for sampled k-mers
    std::vector<vector<uint64_t>> local_dists(numThreads , vector<uint64_t>(kVal + 1 , 0));

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        auto& hist = local_dists[tid];
        #pragma omp for schedule(dynamic , 1)
        for (unsigned int i = 0; i < countRow; i++) {
            const uint128 ki = kmersRow[i];
            for (unsigned int j = 0; j < countCol; j++) {
                int hd = hdPC(ki, kmersCol[j], kVal);
                hist[hd]++;
            }
        }
    }

    for(int t=0; t<numThreads; t++)
    {
        for(int h=0; h<=kVal; h++)
        {
            dists[h] += local_dists[t][h];
        }
    }
}

int main(int argc, char* argv[]) {
    char* inptFile;
    char* outptFile;
    unsigned int seqLen = 0;
    int kVal = 0;
    // Default is full profile with no multithreading
    double theta1 = 1.0;
    double theta2 = 1.0;
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
    // Check for argument errors
    if(kVal > 128 || kVal < 1 || seqLen < kVal || numThreads < 1 || theta1 <= 0 || theta2 <= 0 || theta1 > 1 || theta2 > 1 || inptFile == NULL || outptFile == NULL || inptFile == "" || outptFile == "")
    {
        cout << "ARGUMENT ERROR" << endl;
        return 1;
    }

    maxNum = ((kVal-1) / 32);

    popMask = (2.0)*((pow((long double)4 , (long double)(min(kVal , 32))) - 1)/3.0); // Keep odd bits
    popMask2 = popMask >> 1; // Keep even bits

    // Parse FASTA file to get sequence
    char *charSeq = (char *)malloc(seqLen + 1);
    unsigned int retrievedLen = getSeq(charSeq , inptFile , seqLen);

    if(retrievedLen == 0)
    {
        cout << "INPUT FILE ERROR" << endl;
        return 1;
    }

    // 2-bit encoding
    int* seq = (int *)calloc(retrievedLen , sizeof(int));
    // Store as array of integers, each int is one encoded character
    for (unsigned int i=0; i<retrievedLen; i++)
    {
        seq[i] = seq_nt4_table[(uint8_t)charSeq[i]]; 
    }

    // Tracks how many pairs had a Hamming distance of i, where i is an index of the array
    unsigned int *dists = (unsigned int *)calloc(kVal+1 , sizeof(unsigned int));

    // Run correct method

    // multithreaded
    if(numThreads > 1)
    {
        omp_set_num_threads(numThreads);
        // sketching with multithreading
        if(theta1*theta2 != 1)
        {
            sketch_multithread(kVal , retrievedLen , seq , theta1 , theta2 , numThreads , dists);
        }
        // full profile
        else
        {
            multithread(kVal , retrievedLen , seq , numThreads , dists);
        }
    }
    // sketching
    else if(theta1*theta2 != 1)
    {
        sketch(kVal , retrievedLen , seq , theta1 , theta2 , dists);
    }
    // full profile
    else
    {
        regularVer(kVal , retrievedLen , seq , dists);
    }
    
    output(dists , kVal , (theta1*theta2) , outptFile);

    free(dists);
    free(seq);
    free(charSeq);

    return 0;
}