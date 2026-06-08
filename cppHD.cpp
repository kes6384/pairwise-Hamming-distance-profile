//Algorithm for finding the Hammond Distance profile of a given sequence
//Includes XOR and popcount methods
//Takes FASTA files as input
//Outputs a text file representing the histogram of Hammond Distances between k-mers

#include <iostream>
#include <fstream>
#include <cmath>
#include "tinyFA.hpp"  // FASTA file parser: https://github.com/edawson/tinyFA
#include "pliib.hpp"

#define popcount // XOR or popcount, choose HD calculation method
//#define compress

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
  getSequence(tf, contigName, seq, 0, len);
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
  dist = __builtin_popcount(orResult);
  return dist;
}

/*
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
}*/

// Output array of Hamming distance counts to a text file
// dists - array of Hamming distance counts
// len - length of array of Hamming distance counts (equivalent to k-mer size + 1, since HD ranges from 0 to k)
void output(unsigned long long *dists , int len)
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

//command line args: file name, seqLen, k, kmerList (y, n), sequence name
int main(int argc, char* argv[]) {
  char* file = argv[1];
  int seqLen = atoi(argv[2]); // Sequence length
  int kVal = atoi(argv[3]); // k-mer length
  char* kmerList = argv[4]; // If k-mers should be extracted once and put into a list or extracted in the loop
  //cout << ("%s" , method) << endl;

  char* contigName = argv[5]; // Sequence name

  popMask = (2.0)*((pow(4 , kVal) - 1)/3.0); // Keep odd bits
  popMask2 = popMask >> 1; // Keep even bits

  // Tracks how many pairs had a Hamming distance of i, where i is an index of the array
  unsigned long long *dists = (unsigned long long *)calloc(kVal+1 , sizeof(int));

  //cout << "hi";

  // Parse FASTA file to get sequence
  char *charSeq;
  getSeq(charSeq , file , seqLen , contigName);
  cout << ("%d" , strlen(charSeq)) << endl;
  //cout << "hi";

  //cout << ("%d" , (int)strlen(charSeq)) << endl;
  // 2-bit encoding
  #ifndef compress
  int* seq = (int *)calloc(seqLen , sizeof(int));
  // Store as array of integers, each int is one encoded character
  if (*kmerList == 'y')
  {
    for (int i=0; i<seqLen; i++)
    {
      seq[i] = seq_nt4_table[(uint8_t)charSeq[i]]; 
    }
  }
  #endif
  #ifdef compress
  // Store as array of bytes, each byte contains multiple encoded characters
  // Could be faster to use a larger type instead of 1 byte characters
  char* seq = (char *)calloc(((seqLen)/(sizeof(char)*4))+1 , sizeof(char));
  for (int i=0; i<((seqLen)/(sizeof(char)*4))+1; i++)
  {
    //cout << ("%d" , i) << endl;
    char entry = 0;
    //cout << ("%c" , entry) << endl;
    // One byte can fit 4 2-bit encoded characters
    int count = 0;
    for (int j=i*(sizeof(char)*4); j<seqLen; j++)
    {
      //cout << ("%d" , j) << endl;
      //cout << ("%d" , (uint8_t)charSeq[j]) << endl;
      //cout << ("%f" , seq_nt4_table[(uint8_t)charSeq[j]] ) << endl;
      entry = entry || (seq_nt4_table[(uint8_t)charSeq[j]] << count);
      count ++;
      if(count >= sizeof(char)*4)
        break;
    }
    //cout << ("%s" , entry) << endl;
    seq[i] = entry;
  }
  //cout << ("%s" , *seq) << endl;
  #endif
  //delete [] charSeq;
  
  // Iterate through k-mers and check HD for each pair
  uint64_t mask = (kVal == 32) ? ~0ULL : ((1ULL << (2 * kVal)) - 1);
  //cout << ("%d" , seqLen) << endl;
  if(*kmerList == 'y')
  {
    uint64_t kmer1 = 0;
    for (int i=0; i<kVal-1; i++)
    {
      //int c = seq_nt4_table[(uint8_t)charSeq[i]];
      #ifndef compress
      int c = seq[i];
      #endif
      #ifdef compress
      int c = seq[i/(sizeof(char)*4)];
      c = (c >> (i*2)) % 4;
      //c = c % 2*((i+1)-(i/(sizeof(char)*4)));
      #endif
      kmer1 = ((kmer1 << 2) | (uint64_t)c) & mask;
    }
    //cout << ("%c" , seq) << endl;
    //int c = seq_nt4_table[(uint8_t)seq[0]];
    for (int i=kVal-1; i<seqLen; i++)
    {
      //cout << "hi";
      //int c = seq_nt4_table[(uint8_t)charSeq[i]];
      #ifndef compress
      int c = seq[i];
      #endif
      #ifdef compress
      int c = seq[i/(sizeof(char)*4)];
      c = (c >> (i*2)) % 4;
      #endif
      //int c = (uint8_t)seq[i];
      //cout << ("%d" , c) << endl;
      //cout << "?" << endl;
      kmer1 = ((kmer1 << 2) | (uint64_t)c) & mask;
      //cout << "hi";
      //uint64_t kmer2 = 0;
      uint64_t kmer2 = kmer1;
      /*for (int j=i-kVal+1; j<i+1; j++)
      {
        int c = seq_nt4_table[(uint8_t)seq[j]];
        kmer2 = ((kmer2 << 2) | (uint64_t)c) & mask;
      }*/
      //cout << ("%d" , kmer1) << endl;
      //cout << ("%s" , charSeq[i]) << endl;
      for (int j=i+1; j<seqLen; j++)
      {
        //cout << "hi";
        //int c2 = seq_nt4_table[(uint8_t)charSeq[j]];
        #ifndef compress
        int c2 = seq[j];
        #endif
        #ifdef compress
        int c2 = seq[j/(sizeof(char)*4)];
        c2 = (c >> (j*2)) % 4;
        #endif
        kmer2 = ((kmer2 << 2) | (uint64_t)c2) & mask;
        //cout << (kmer1) << "," << (kmer2) << endl;
        #ifdef XOR
          dists[hdXOR(kmer1  , kmer2 , kVal)] ++;
        #endif
        #ifdef popcount
          dists[hdPC(kmer1  , kmer2 , kVal)] ++;
        #endif
        //cout << ("%d" , HD(kmer1  , kmer2 , kVal)) << endl;
      }
    }
  }
  else
  {
    uint64_t kmer1 = 0;
    for (int i=0; i<kVal-1; i++)
    {
      int c = seq_nt4_table[(uint8_t)charSeq[i]];
      kmer1 = ((kmer1 << 2) | (uint64_t)c) & mask;
    }
    //cout << ("%c" , seq) << endl;
    //int c = seq_nt4_table[(uint8_t)seq[0]];
    for (int i=kVal-1; i<seqLen; i++)
    {
      //cout << "hi";
      int c = seq_nt4_table[(uint8_t)charSeq[i]];
      //int c = (uint8_t)seq[i];
      //cout << ("%d" , c) << endl;
      //cout << "?" << endl;
      kmer1 = ((kmer1 << 2) | (uint64_t)c) & mask;
      //cout << "hi";
      //uint64_t kmer2 = 0;
      uint64_t kmer2 = kmer1;
      /*for (int j=i-kVal+1; j<i+1; j++)
      {
        int c = seq_nt4_table[(uint8_t)seq[j]];
        kmer2 = ((kmer2 << 2) | (uint64_t)c) & mask;
      }*/
      //cout << ("%d" , kmer1) << endl;
      //cout << ("%s" , charSeq[i]) << endl;
      for (int j=i+1; j<seqLen; j++)
      {
        //cout << "hi";
        int c2 = seq_nt4_table[(uint8_t)charSeq[j]];
        kmer2 = ((kmer2 << 2) | (uint64_t)c2) & mask;
        //cout << (kmer1) << "," << (kmer2) << endl;
        #ifdef XOR
          dists[hdXOR(kmer1  , kmer2 , kVal)] ++;
        #endif
        #ifdef popcount
          dists[hdPC(kmer1  , kmer2 , kVal)] ++;
        #endif
        //cout << ("%d" , HD(kmer1  , kmer2 , kVal)) << endl;
      }
    }
  }
  
  free(seq);
  delete [] charSeq;
  output(dists , kVal);
  free(dists);
  return 0;
}