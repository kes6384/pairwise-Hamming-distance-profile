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
void getSeq(char* seq, char* file, int len)
{
  getSequence(file, len, seq);
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

    int numKmers = (seqLen - kVal) + 1;

    // Tracks how many pairs had a Hamming distance of i, where i is an index of the array
    unsigned int *dists = (unsigned int *)calloc(kVal+1 , sizeof(int));

    // Matrices for each base
    Eigen::MatrixXf a(numKmers , kVal);
    Eigen::MatrixXf t(numKmers , kVal);
    Eigen::MatrixXf c(numKmers , kVal);
    Eigen::MatrixXf g(numKmers , kVal);

     // Parse FASTA file to get sequence
    char *charSeq = (char *)malloc(seqLen + 1);
    getSeq(charSeq , file , seqLen);
    
    // Initialize matrices
    /* naive method
    for (int i=0; i<numKmers; i++)
    {
        for (int j=i; j<kVal+i; j++)
        {
            int index = j-i; 
            char character = charSeq[j];
            a(i,index) = a_table[(uint8_t)character];
            c(i,index) = c_table[(uint8_t)character];
            g(i,index) = g_table[(uint8_t)character];
            t(i,index) = tu_table[(uint8_t)character];
        }
    }*/

    // First row = first k-mer
    for (int i=0; i<kVal; i++)
    {
      uint8_t character = (uint8_t)charSeq[i];
      a(0,i) = (uint64_t)a_table[character];
      c(0,i) = (uint64_t)c_table[character];
      g(0,i) = (uint64_t)g_table[character];
      t(0,i) = (uint64_t)tu_table[character];
    }
    //cout << a.row(0) << endl;
    /*cout << c.row(0) << endl;
    cout << g.row(0) << endl;
    cout << t.row(0) << endl;*/
    // Next k-mer is previous row shifted to the left with next character in last column
    for (int i=1; i<numKmers; i++)
    {
      uint8_t character = (uint8_t)charSeq[kVal+i-1];

      a.block(i,0,1,kVal-1) = a.block(i-1,1,1,kVal-1);
      a(i,kVal-1) = (uint64_t)a_table[character];
      //cout << a.row(i) << endl;

      c.block(i,0,1,kVal-1) = c.block(i-1,1,1,kVal-1);
      c(i,kVal-1) = (uint64_t)c_table[character];

      g.block(i,0,1,kVal-1) = g.block(i-1,1,1,kVal-1);
      g(i,kVal-1) = (uint64_t)g_table[character];

      t.block(i,0,1,kVal-1) = t.block(i-1,1,1,kVal-1);
      t(i,kVal-1) = (uint64_t)tu_table[character];
    }

    Eigen::MatrixXf m(numKmers , numKmers);
    
    //m = m + AeBe for all e in sigma, where B = A transpose
    m = m + a*a.transpose();
    m = m + t*t.transpose();
    m = m + c*c.transpose();
    m = m + g*g.transpose();

    // Use number of matching characters to get number of mismatched characters
    for(int i=1; i<numKmers; i++)
    {
        for(int j=1; j<numKmers; j++)
        {
            int dist = kVal - m(i,j);
            dists[dist] ++;
            //cout << ("%d" , dist) << endl;
        }
    }

  delete [] charSeq;
  output(dists , kVal);
  free(dists);
  return 0;
}