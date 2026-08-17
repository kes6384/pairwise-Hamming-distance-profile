//Algorithm for finding the Hammond Distance profile of a given sequence
//Uses popcount method and creates a sketch
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
#include "MurmurHash3.cpp"
#include "parseFASTA.cpp" // Simple FASTA parser
//#include "tinyFA.hpp"  // FASTA file parser: https://github.com/edawson/tinyFA
//#include "pliib.hpp"

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
  //cout << ("%s" , seq);
}

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
// len - length of array of Hamming distance counts (equivalent to k-mer size + 1, since HD ranges from 0 to k)
// rate - sampling rate (theta1 * theta2)
void output(unsigned int *dists , int len , float rate , char* outFile)
{
  ofstream out(outFile);

  //cout << ("%f" , rate) << endl;
  out << "Hamming Distance : Number of Pairs" << endl;
  for(int i=0; i<=len; i++)
  {
      out << ("%d" , i) << (" : ") << ("%d" , (int)((float)dists[i]/(2.0 * rate))) << endl;
  }

  out.close();
}

//command line args: input file name, seqLen, k, sampling rate 1, 2, output file
int main(int argc, char* argv[]) {
  char* file = argv[1];
  int seqLen = atoi(argv[2]); // Sequence length
  int kVal = atoi(argv[3]); // k-mer length
  double theta1 = atof(argv[4]); // sampling rate
  double theta2 = atof(argv[5]); // sampling rate
  std::random_device rd;
  uint32_t seed1 = rd();
  uint32_t seed2 = rd();

  popMask = (2.0)*((pow((long double)4 , (long double)kVal) - 1)/3.0); // Keep odd bits
  popMask2 = popMask >> 1; // Keep even bits

  // Tracks how many pairs had a Hamming distance of i, where i is an index of the array
  unsigned int *dists = (unsigned int *)calloc(kVal+1 , sizeof(int));

  // Parse FASTA file to get sequence
  char *charSeq = (char *)malloc(seqLen + 1);
  int retrievedLen = getSeq(charSeq , file , seqLen);

  // 2-bit encoding
  int* seq = (int *)calloc(retrievedLen , sizeof(int));
  // Store as array of integers, each int is one encoded character
    for (int i=0; i<retrievedLen; i++)
    {
      seq[i] = seq_nt4_table[(uint8_t)charSeq[i]]; 
    }
  
  // Sample k-mers
  int numKmers = (retrievedLen - kVal) + 1;
  uint64_t *kmersRow = (uint64_t*)calloc(numKmers , sizeof(uint64_t));
  uint64_t *kmersCol = (uint64_t*)calloc(numKmers , sizeof(uint64_t));

  uint64_t mask = (kVal == 32) ? ~0ULL : ((1ULL << (2 * kVal)) - 1);
    //cout << ("%s" , charSeq) << endl;
    uint64_t kmer1 = 0;
    for (int i=0; i<kVal-1; i++)
    {
      int c = seq[i];
      kmer1 = ((kmer1 << 2) | (uint64_t)c) & mask;
    }
    int countRow = 0;
    int countCol = 0;
    uint32_t *hashed = (uint32_t *)malloc(sizeof(uint32_t));
    for (int i=kVal-1; i<retrievedLen; i++)
    {
      int c = seq[i];
      kmer1 = ((kmer1 << 2) | (uint64_t)c) & mask;
      MurmurHash3_x86_32(&kmer1 , sizeof(kmer1) , seed1 , hashed);
      float hash = (float)(*hashed) / (float)(UINT32_MAX);
      //float hashed = (float)(std::hash<uint64_t>{}(kmer1))/(float)UINT64_MAX;
      //cout << ("%f" , hash) << endl;
      if(hash < theta1)
      {
        kmersRow[countRow] = kmer1;
        countRow ++;
      }
      // Different seeds for row and column for independence
      MurmurHash3_x86_32(&kmer1 , sizeof(kmer1) , seed2 , hashed);
      hash = (float)(hash) / (float)(UINT32_MAX);
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
  
  free(seq);
  free(charSeq);
  free(kmersRow);
  free(kmersCol);
  output(dists , kVal , (theta1*theta2) , argv[6]);
  free(dists);
  return 0;
}