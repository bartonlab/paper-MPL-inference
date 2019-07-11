#include <vector>
#include <assert.h>

#include "inf.h"    // inference declarations
#include "io.h"     // input/output



/**********  I N P U T  **********/

// Retrieve sequences, counts, and generation times from a file (format in main.cpp and examples/)
// FUTURE: Make robust to incorrect time ordering by replacing IntVectors with maps from double (time) to IntVector

void getSequences(FILE *input, IntVVector &sequences, Vector &counts, std::vector<double> &times) {

    char c;
    int i;
    double d;
    
    while (fscanf(input,"%lf",&d)==1) {

        if (times.size()==0 || times.back()!=d) {
        
            times.push_back(d);
            sequences.push_back(IntVector());
            counts.push_back(std::vector<double>());
            
        }
        
        fscanf(input,"%lf",&d);
        counts.back().push_back(d);
        sequences.back().push_back(std::vector<int>());
        
        while (fscanf(input,"%c",&c)==1) {
    
            if (c=='\n' || c=='\r') break;
            
            fscanf(input,"%d",&i);
            (sequences.back()).back().push_back(i);
            
        }
        
    }
    
    // normalize counts
    
    for (int i=0;i<counts.size();i++) {
    
        double total = 0;
        for (int j=0;j<counts[i].size();j++) total        += counts[i][j];
        for (int j=0;j<counts[i].size();j++) counts[i][j] /= total;
        
    }
    
    // sanity check
    
    int length = sequences[0][0].size();
    
    for (int i=0;i<sequences.size();i++) { for (int j=0;j<sequences[i].size();j++) {
    
        assert(length==sequences[i][j].size() && "io.cpp - getSequences: read in sequence vectors with unequal lengths");
        
    } }
    
    // debug
    assert(sequences.size()==counts.size() && "io.cpp - getSequences: unequal lengths for sequences and counts");
    assert(sequences.size()==times.size()  && "io.cpp - getSequences: unequal lengths for sequences and times");
    
}


// Expand q-nary sequences into binary form

void expandSequences(int q, IntVVector &sequences) {
    
    IntVVector expandedSequences(sequences.size(), IntVector());
    
    for (int i=0;i<sequences.size();i++) { for (int j=0;j<sequences[i].size();j++) {
    
        expandedSequences[i].push_back(std::vector<int>(sequences[i][j].size() * (q-1), 0));
//        expandedSequences[i].push_back(std::vector<int>(sequences[i][j].size() * q, 0));
        
        for (int k=0;k<sequences[i][j].size();k++) {
        
            if (sequences[i][j][k]!=0) expandedSequences[i].back()[(k * (q-1)) + sequences[i][j][k] - 1] = 1;
//            expandedSequences[i].back()[(k * q) + sequences[i][j][k]] = 1;
            
        }
    
    } }
    
    sequences = expandedSequences;

}


// Retrieve a mutation matrix from a file

void getMu(FILE *input, Vector &mu) {

    char c;
    int i;
    double d;
    
    while (fscanf(input,"%lf",&d)==1) {

        mu.push_back(std::vector<double>());
        mu.back().push_back(d);
        
        while (fscanf(input,"%c",&c)==1) {
    
            if (c=='\n' || c=='\r') break;
            
            fscanf(input,"%lf",&d);
            mu.back().push_back(d);
            
        }
        
    }
    
    // sanity check
    
    int length = mu.size();
    
    for (int i=0;i<mu.size();i++) {
    
        assert(length==mu[i].size() && "io.cpp - getMu: read in mutation matrix that is not square");
        
    }
    
}


/**********  O U T P U T  **********/


// Print selection coefficients to file

void printSelectionCoefficients(FILE *output, const std::vector<double> &s) {
    
    for (int i=0;i<s.size();i++) fprintf(output,"%.6e\n",s[i]);
    fflush(output);

}

void printCovariance(FILE* output, double cov[], int L){
    
    for (int i=0;i<L;i++) {
    
        fprintf(output,"%.4e",cov[(i * L)]);
        for (int j=1;j<L;j++) fprintf(output," %.4e",cov[(i * L) + j]);
        fprintf(output,"\n");
    
    }
    fflush(output);

}

void printNumerator(FILE* output, double dx[], int L){
    
    for (int i=0;i<L;i++) {
        
        fprintf(output,"%.4e\n",dx[i]);
        
    }
    fflush(output);
    
}


/**********  A U X I L I A R Y  **********/


// Insert an element into vector in order (smallest to largest)

void insertInPlace(std::vector<double> &list, double element) {

    if (list.empty()) list.push_back(element);
    
    else {
    
        bool insert=true;
        int  imin=-1;
        int  imax=(int)list.size();
    
        while (imax-imin>1) {
      
            int imid = (imin+imax)/2;
 
            if      (list[imid] <  element) imin = imid;
            else if (list[imid] >  element) imax = imid;
            else if (list[imid] == element) { insert=false; break; }
        
        }
        
        if (insert) list.insert(list.begin()+imax,element);
        
    }

}

