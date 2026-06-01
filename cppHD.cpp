//Algorithm for finding the Hammond Distance profile of a given sequence
//Includes XOR and popcount methods
//Takes FASTA files as input
//Outputs a text file representing the histogram of Hammond Distances between k-mers

#include <iostream>
#include <fstream>
#include <cmath>
#include "tinyFA.hpp"  // FASTA file parser: https://github.com/edawson/tinyFA
#include "pliib.hpp"

using namespace std;
using namespace TFA;

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

// Use tinyFA library to parse FASTA file
// Sequence is stored in seq and is of length len
void getSeq(char*& seq, char* file, int len, char* contigName)
{
  tiny_faidx_t tf;
  // Check if an index exists, and create one if not.
  if (!checkFAIndexFileExists(file)){
      createFAIndex(file, tf);
  }
  else{
      // Parses an FAI file when passed a FASTA file name.
      parseFAIndex(file, tf);
  }
  getSequence(tf, contigName, seq, 0, len-1);
  //cout << "!" << endl;
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
  int dist = __builtin_popcount(orResult);
  return dist;
}

// Call the requested method to compute HD on kmer1 and kmer2
int HD(uint64_t kmer1 , uint64_t kmer2 , int k, string method)
{
  int hd = 0;
  //cout << ("%s" , method) << endl;

  if(method == "XOR")
  {
    //cout << "!" << endl;
    hd = hdXOR(kmer1 , kmer2 , k);
  }
  else if(method == "popcount")
  {
    hd = hdPC(kmer1 , kmer2 , k);
  }
  
  return hd;
}

// Output array of Hamming distance counts to a text file
// dists - array of Hamming distance counts
// len - length of array of Hamming distance counts (equivalent to k-mer size + 1, since HD ranges from 0 to k)
void output(int *dists , int len)
{
  ofstream out("HD_cpp_out.txt");

  out << "Hamming Distance : Number of Pairs" << endl;
  for(int i=0; i<=len; i++)
  {
      //out << ("%d : %d" , i , dists[i]) << endl;
      out << ("%d" , i) << (" : ") << ("%d" , dists[i]) << endl;
      //cout << ("%d" , dists[i]) << endl;
  }

  out.close();
}

//command line args: file name, seqLen, k, method (XOR, popcount), kmerList (y, n), sequence name
int main(int argc, char* argv[]) {
  char* file = argv[1];
  int seqLen = atoi(argv[2]); // Sequence length
  int kVal = atoi(argv[3]); // k-mer length
  char* method = argv[4]; // Which method to use to compute HD (XOR, popcount)
  char* kmerList = argv[5]; // If k-mers should be extracted once and put into a list or extracted in the loop
  //cout << ("%s" , method) << endl;

  char* contigName = argv[6]; // Sequence name

  popMask = (2.0)*((pow(4 , kVal) - 1)/3.0); // Keep odd bits
  popMask2 = popMask >> 1; // Keep even bits

  // Tracks how many pairs had a Hamming distance of i, where i is an index of the array
  int *dists = (int *)calloc(kVal+1 , sizeof(int));

  //cout << "hi";

  // Parse FASTA file to get sequence
  char *seq;
  getSeq(seq , file , seqLen , contigName);
  //cout << "hi";
  
  // Iterate through k-mers and check HD for each pair
  uint64_t mask = (kVal == 32) ? ~0ULL : ((1ULL << (2 * kVal)) - 1);
  //cout << ("%d" , seqLen) << endl;
  if(*kmerList == 'n')
  {
    uint64_t kmer1 = 0;
    for (int i=0; i<kVal-1; i++)
    {
      int c = seq_nt4_table[(uint8_t)seq[i]];
      //int c = seq[i];
      kmer1 = ((kmer1 << 2) | (uint64_t)c) & mask;
    }
    //cout << ("%c" , seq) << endl;
    //int c = seq_nt4_table[(uint8_t)seq[0]];
    //cout << ("%c" , seq) << endl;
    for (int i=kVal-1; i<seqLen; i++)
    {
      //cout << "hi";
      int c = seq_nt4_table[(uint8_t)seq[i]];
      //int c = seq[i];
      //int c = (uint8_t)seq[i];
      //cout << ("%d" , c) << endl;
      //cout << "?" << endl;
      kmer1 = ((kmer1 << 2) | (uint64_t)c) & mask;
      //cout << "hi";
      uint64_t kmer2 = 0;
      for (int j=i-kVal+1; j<i+1; j++)
      {
        int c = seq_nt4_table[(uint8_t)seq[j]];
        kmer2 = ((kmer2 << 2) | (uint64_t)c) & mask;
      }
      for (int j=i+1; j<seqLen; j++)
      {
        //cout << "hi";
        int c2 = seq_nt4_table[(uint8_t)seq[j]];
        kmer2 = ((kmer2 << 2) | (uint64_t)c2) & mask;
        //cout << (kmer1) << "," << (kmer2) << endl;
        dists[HD(kmer1  , kmer2 , kVal , method)] ++;
        //cout << ("%d" , HD(kmer1  , kmer2 , kVal , method)) << endl;
      }
    }
  }
  else
  {
    int numKmers = (seqLen-kVal)+1;
    uint64_t *kmers = (uint64_t *)calloc(numKmers , kVal);
    //cout << ("%d" , (seqLen-kVal)+1) << endl;
    uint64_t kmer1 = 0;
    for (int i=0; i<kVal-1; i++)
    {
      int c = seq_nt4_table[(uint8_t)seq[i]];
      //int c = seq[i];
      kmer1 = ((kmer1 << 2) | (uint64_t)c) & mask;
    }
    // extract all k-mers
    int count = 0;
    for(int i=kVal-1; i<seqLen; i++)
    {
      int c = seq_nt4_table[(uint8_t)seq[i]];
      kmer1 = ((kmer1 << 2) | (uint64_t)c) & mask;
      kmers[count] = kmer1;
      //cout << ("%d" , kmers[numKmers]) << endl;
      count ++;
    }
    // Check HD of all k-mer pairs
    for (int i=0; i<numKmers-1; i++)
    {
      kmer1 = kmers[i];
      uint64_t kmer2 = 0;
      for (int j=i+1; j<numKmers; j++)
      {
        kmer2 = kmers[j];
        //cout << (kmer1) << "," << (kmer2) << endl;
        dists[HD(kmer1  , kmer2 , kVal , method)] ++;
        //cout << ("%d" , HD(kmer1  , kmer2 , kVal , method)) << endl;
      }
    }
  }
  
  delete [] seq;
  output(dists , kVal);
  free(dists);
  return 0;
}