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

    // Matrices
    // Square matrices with k columns
    // n = r
    Eigen::MatrixXf a(numKmers , kVal);
    Eigen::MatrixXf t(numKmers , kVal);
    Eigen::MatrixXf c(numKmers , kVal);
    Eigen::MatrixXf g(numKmers , kVal);
    // Transpose matrices
    Eigen::MatrixXf at(kVal , numKmers);
    Eigen::MatrixXf tt(kVal , numKmers);
    Eigen::MatrixXf ct(kVal , numKmers);
    Eigen::MatrixXf gt(kVal , numKmers);
    /* This part is for when numKmers != kVal. Would also effect later code since have to decompose into squares
    // n > r
    // Change n to the next smallest multiple of r
    while (numKmers % kVal != 0)
    {
        numKmers ++;
    }
    // n < r
    if(kVal % numKmers != 0)
    {

    }*/

     // Parse FASTA file to get sequence
    char *charSeq = (char *)malloc(seqLen + 1);
    getSeq(charSeq , file , seqLen);
    
    // Initialize matrices
    for (int i=0; i<numKmers; i++)
    {
        for (int j=i; j<kVal+i; j++)
        {
            int index = j-i; 
            char character = charSeq[j];
            a(i,index) = 0;
            t(i,index) = 0;
            c(i,index) = 0;
            g(i,index) = 0;
            at(index,i) = 0;
            tt(index,i) = 0;
            ct(index,i) = 0;
            gt(index,i) = 0;
            if(character == 'A')
            {
                a(i,index) = 1;
                at(index,i) = 1;
            }
            if(character == 'T')
            {
                t(i,index) = 1;
                tt(index,i) = 1;
            }
            if(character == 'C')
            {
                c(i,index) = 1;
                ct(index,i) = 1;
            }
            if(character == 'G')
            {
                g(i,index) = 1;
                gt(index,i) = 1;
            }
        }
    }

    Eigen::MatrixXf m(numKmers , numKmers);
    
    //m = m + AeBe for all e in sigma
    m = m + a*at;
    m = m + t*tt;
    m = m + c*ct;
    m = m + g*gt;

    for(int i=1; i<numKmers; i++)
    {
        for(int j=1; j<numKmers; j++)
        {
            int dist = kVal - m(i,j);
            dists[dist] ++;
        }
    }

  delete [] charSeq;
  output(dists , seqLen);
  return 0;
}