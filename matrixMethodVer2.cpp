//Algorithm for finding the Hammond Distance profile of a given sequence
//Uses matrix multiplication in an attempt to inmprove speed
//Takes FASTA files as input
//Outputs a text file representing the histogram of Hammond Distances between k-mers

#include <iostream>
#include <fstream>
#include <cmath>
#include "parseFASTA.cpp" // Simple FASTA parser
// https://libeigen.gitlab.io/
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

using namespace std;

const unsigned char a_table[256] = { // translate A to 1
	1, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 1, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 1, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
};

const unsigned char c_table[256] = { // translate C to 1
	0, 1, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 1,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 1,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
};

const unsigned char g_table[256] = { // translate G to 1
	0, 0, 1, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 1,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 1,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
};

const unsigned char tu_table[256] = { // translate T/U to 1
	0, 0, 0, 1,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  1, 1, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  1, 1, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
};

// parse FASTA file using parseFASTA.cpp
// Sequence is stored in seq and is of length len
int getSeq(char* seq, char* file, int len)
{
  return getSequence(file, len, seq);
}

// Output array of Hamming distance counts to a text file
// dists - array of Hamming distance counts
// len - length of array of Hamming distance counts (equivalent to k-mer size + 1, since HD ranges from 0 to k)
void output(unsigned int *dists , int len)
{
  ofstream out("matrix_hd_out.txt");

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

    // Tracks how many pairs had a Hamming distance of i, where i is an index of the array
    //unsigned int *dists = (unsigned int *)calloc(kVal+1 , sizeof(int));

    // Tracks how many pairs had a Hamming distance of i, where i is an index of the array
    unsigned int *dists = (unsigned int *)calloc(kVal+1 , sizeof(int));

     // Parse FASTA file to get sequence
    char *charSeq = (char *)malloc(seqLen + 1);
    int retrievedLen = getSeq(charSeq , file , seqLen);
    int numKmers = (retrievedLen - kVal) + 1;

    // One large matrix containing masked k-mers for each base
    // stored in order: a, c, g, t
    Eigen::MatrixXf kmers(numKmers , kVal*4);

    // First row = first k-mer
    for (int i=0; i<kVal; i++)
    {
      uint8_t character = (uint8_t)charSeq[i];
      kmers(0,i) = (float)a_table[character]; //a
      kmers(0,(i+kVal)) = (float)c_table[character]; //c
      kmers(0,(i+(2*kVal))) = (float)g_table[character]; //g
      kmers(0,(i+(3*kVal))) = (float)tu_table[character]; //t
    }

    // Next k-mer is previous row shifted to the left with next character in last column
    for (int i=1; i<numKmers; i++)
    {
      uint8_t character = (uint8_t)charSeq[kVal+i-1];

      kmers.block(i,0,1,((kVal*4)-1)) = kmers.block(i-1,1,1,((kVal*4)-1));

      //kmers.block(i,0,1,kVal-1) = kmers.block(i-1,1,1,kVal-1);
      kmers(i,kVal-1) = (float)a_table[character];

      //kmers.block(i,kVal,1,kVal-1) = kmers.block(i-1,kVal+1,1,kVal-1);
      kmers(i,(2*kVal)-1) = (float)c_table[character];

      //kmers.block(i,2*kVal,1,kVal-1) = kmers.block(i-1,(2*kVal)+1,1,kVal-1);
      kmers(i,(3*kVal)-1) = (float)g_table[character];

      //kmers.block(i,3*kVal,1,kVal-1) = kmers.block(i-1,(3*kVal)+1,1,kVal-1);
      kmers(i,(4*kVal)-1) = (float)tu_table[character];
    }

    delete [] charSeq;

    Eigen::MatrixXf m(numKmers , numKmers);

    //m = m + AeBe for all e in sigma, where B = A transpose
    m = (kmers*kmers.transpose()).triangularView<Eigen::Lower>();

    // Use number of matching characters to get number of mismatched characters
    for(int i=0; i<numKmers; i++)
    {
        // Only look at lower triangle for results
        for(int j=i+1; j<numKmers; j++)
        {
          int dist = kVal - m(j,i);
          dists[dist] ++;
        }
    }

  output(dists , kVal);
  free(dists);
  return 0;
}