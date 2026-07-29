//Algorithm for finding the Hammond Distance profile of a given sequence
//Multithreaded popcount method
//Takes FASTA files as input
//Outputs a text file representing the histogram of Hammond Distances between k-mers

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <bit>
#include <bitset>
#include <omp.h>
#include <vector>
#include "parseFASTA.cpp" // Simple FASTA parser
//#include "tinyFA.hpp"  // FASTA file parser: https://github.com/edawson/tinyFA
//#include "pliib.hpp"

#define POPCOUNT // XOR or popcount, choose HD calculation method
#if defined(XOR)
    #define HD_FUNC hdXOR
#elif defined(POPCOUNT)
    #define HD_FUNC hdPC
#else
    #define HD_FUNC HD
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

// Compute Hamming distance between kmer1 and kmer2
// Method simply uses XOR on the two bit-encoded k-mers
int hdXOR(uint64_t kmer1 , uint64_t kmer2 , int k)
{
  uint64_t result = kmer1 ^ kmer2;

  int dist = 0;
  for(int i=0; i<k; i++)
  {
      // Check character
      int nxtChar = result % 4;
      if(nxtChar != 0)
          dist ++;
      // Move to next character
      result = result >> 2;
  }
  return dist;
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
void output(std::vector<uint64_t> &dists , int len)
{
  ofstream out("HD_cpp_OMP.txt");

  out << "Hamming Distance : Number of Pairs" << endl;
  for(int i=0; i<=len; i++)
  {
      out << ("%d" , i) << (" : ") << ("%d" , dists[i]) << endl;
  }

  out.close();
}

//command line args: file name, seqLen, k
int main(int argc, char* argv[]) {
  char* file = argv[1];
  int seqLen = atoi(argv[2]); // Sequence length
  int kVal = atoi(argv[3]); // k-mer length
  char kmerList = 'y'; // If sequence should be bit-encoded before k-mers are extracted

  int numThreads = atoi(argv[4]);
  omp_set_num_threads(numThreads);

  popMask = (2.0)*((pow((long double)4 , (long double)kVal) - 1)/3.0); // Keep odd bits
  popMask2 = popMask >> 1; // Keep even bits

  // Tracks how many pairs had a Hamming distance of i, where i is an index of the array
  //unsigned int *dists = (unsigned int *)calloc(kVal+1 , sizeof(int));

  // Parse FASTA file to get sequence
  char *charSeq = (char *)malloc(seqLen + 1);
  int retrievedLen = getSeq(charSeq , file , seqLen);

  // 2-bit encoding
  int* seq = (int *)calloc(retrievedLen , sizeof(int));
  // Store as array of integers, each int is one encoded character
  if (kmerList == 'y')
  {
    for (int i=0; i<retrievedLen; i++)
    {
      seq[i] = seq_nt4_table[(uint8_t)charSeq[i]]; 
    }
  }

  int numKmers = (retrievedLen - kVal) + 1;
  uint64_t *kmers = (uint64_t*)calloc(numKmers , sizeof(uint64_t));

  // Iterate through k-mers and check HD for each pair
  uint64_t mask = (kVal == 32) ? ~0ULL : ((1ULL << (2 * kVal)) - 1);
  if(kmerList == 'y')
  {
    //cout << ("%s" , charSeq) << endl;
    uint64_t kmer1 = 0;
    for (int i=0; i<kVal-1; i++)
    {
      int c = seq[i];
      kmer1 = ((kmer1 << 2) | (uint64_t)c) & mask;
    }
    int count = 0;
    for (int i=kVal-1; i<retrievedLen; i++)
    {
      int c = seq[i];
      kmer1 = ((kmer1 << 2) | (uint64_t)c) & mask;
      kmers[count] = kmer1;
      count ++;

      /*uint64_t kmer2 = kmer1;
      for (int j=i+1; j<retrievedLen; j++)
      {
        int c2 = seq[j];
        kmer2 = ((kmer2 << 2) | (uint64_t)c2) & mask;
        dists[HD_FUNC(kmer1  , kmer2 , kVal)] ++;
      }*/
    }
  }

  /*for(int i=0; i<numKmers; i++)
  {
    #pragma omp parallel for
    for(int j=i+1; j<numKmers; j++)
    {
      dists[HD_FUNC(kmers[i] , kmers[j] , kVal)] ++;
    }
  }*/

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
  
  free(seq);
  free(charSeq);
  free(kmers);
  output(dists , kVal);
  return 0;
}