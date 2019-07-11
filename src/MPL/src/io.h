#ifndef IO_H
#define IO_H


#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdio.h>

#include "inf.h"    // inference declarations


// Input
void getSequences(FILE *input, IntVVector &sequences, Vector &counts, std::vector<double> &times);
void expandSequences(int q, IntVVector &sequences);
void getMu(FILE *input, Vector &mu);

// Output
void printSelectionCoefficients(FILE *output, const std::vector<double> &s);
void printCovariance(FILE* output, double cov[], int L);
void printNumerator(FILE* output, double dx[], int L);

// Auxiliary
void insertInPlace(std::vector<double> &, double);


// STRING MANIPULATION

// Converts generic to string

template <class T>

inline std::string tostring (const T& t) {
    
    std::stringstream ss;
    ss << t;
    return ss.str();
    
}


// Converts a string to an integer

inline int strtoint(const std::string &s) {
    
    std::istringstream i(s);
    int x;
    
    if (!(i >> x)) printf("String to integer conversion failed!");
    return x;
    
}


// Converts a string to a double

inline double strtodouble(const std::string &s) {
    
    std::istringstream i(s);
    double x;
    
    if (!(i >> x)) printf("String to double conversion failed!");
    return x;

}


#endif
