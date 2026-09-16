#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdint>
#include <cstring>

using namespace std;

// Gets a sequence of length seqLen from the given file
// seqLen - length of sequence to get
// seq - char array to store sequence in
// returns the length of the sequence retrieved
int getSequence(char* file , int seqLen , char* seq)
{
    ifstream sequence(file);

    if(!sequence)
        return 0;

    int charNum = 0; // Number of characters read so far
    while ((charNum < seqLen) && (sequence.peek() != EOF))
    {
        char nxtChar = sequence.get();
        // check for header lines
        if(nxtChar == '>')
        {
            // move file pointer to next line
            while(nxtChar!=10 && nxtChar!=13 && (sequence.peek() != EOF))
            {
                nxtChar = sequence.get();
            }
            if(sequence.peek() == EOF)
                break;
        }
        if(nxtChar != 10 && nxtChar != 13)
        {
            seq[charNum] = nxtChar;
            charNum ++;
        }
    }
    seq[charNum] = '\0';

    sequence.close();
    return charNum;
}