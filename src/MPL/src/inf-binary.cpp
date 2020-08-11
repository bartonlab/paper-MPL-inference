#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <chrono>

#include "inf.h"    // inference declarations
#include "io.h"     // input/output


typedef std::chrono::high_resolution_clock Clock;
bool useDebug = false;


// Compute single and pair allele frequencies from binary sequences and counts

void computeAlleleFrequencies(const IntVector &sequences, // vector of sequence vectors
                              const std::vector<double> &counts, // vector of sequence counts
                              int q, // number of states (e.g., number of nucleotides or amino acids)
                              std::vector<double> &p1, // single allele frequencies
                              std::vector<double> &p2 // pair allele frequencies
                              ) {
    
    // Set frequencies to zero
    
    for (int a=0;a<p1.size();a++) p1[a] = 0;
    for (int a=0;a<p2.size();a++) p2[a] = 0;

    int L = (int) p1.size();

    // Iterate through sequences and count the frequency of each state at each site,
    // and the frequency of each pair of states at each pair of sites
    
    for (int k=0;k<sequences.size();k++) {
    
        // NOTE: If we want to change convention so that mutation matrix is NOT shifted wrt frequencies,
        // then we need to count sequences[k][.]!=q-1 here, instead of !=0
        
        for (int i=0;i<sequences[k].size();i++) { if (sequences[k][i]!=0) {

            int a = (i * (q-1)) + sequences[k][i] - 1; // map from q-value sequence[k][.] to position in binary vector

            p1[a] += counts[k];

            for (int j=i+1;j<sequences[k].size();j++) { if (sequences[k][j]!=0) {

                int b = (j * (q-1)) + sequences[k][j] - 1;

                p2[(a * L) + b] += counts[k]; // note: pair frequencies are
                p2[(b * L) + a] += counts[k]; // symmetric so must update both

            } }

        } }
        
    }

}


// Interpolate single and pair allele frequencies between two points, weighted by (1-x) and x, respectively

void interpolateFrequencies(const std::vector<double> &p1_0, // single allele frequencies, 1st time point
                            const std::vector<double> &p2_0, // pair allele frequencies, 1st time point
                            const std::vector<double> &p1_1, // single allele frequencies, 2nd time point
                            const std::vector<double> &p2_1, // pair allele frequencies, 2nd time point
                            double x, // interpolation weight (0.5 = even split between 1st and 2nd times)
                            std::vector<double> &p1, // interpolated single allele frequencies
                            std::vector<double> &p2 // interpolated pair allele frequencies
                            ) {

    int L = (int) p1_0.size();

    // Iterate through states and perform simple linear interpolation of frequencies
    
    for (int a=0;a<L;a++) {
        
        p1[a] = ((1 - x) * p1_0[a]) + (x * p1_1[a]);
            
        for (int b=a+1;b<L;b++) {
        
            double xp2 = ((1 - x) * p2_0[(a * L) + b]) + (x * p2_1[(a * L) + b]);
            
            p2[(a * L) + b] = xp2;
            p2[(b * L) + a] = xp2;
            
        }
            
    }

}


// Update the summed covariance matrix

void updateCovariance(double dg, // time step
                      const std::vector<double> &p1, // single allele frequencies
                      const std::vector<double> &p2, // pair allele frequencies
                      double totalCov[] // integrated covariance matrix
                      ) {

    int L = (int) p1.size();
                          
    // Iterate through states and add contributions to covariance matrix

    for (int a=0;a<L;a++) {
        
        totalCov[(a * L) + a] += dg * p1[a] * (1 - p1[a]);
            
        for (int b=a+1;b<L;b++) {
        
            double dCov = dg * (p2[(a * L) + b] - (p1[a] * p1[b]));
            
            totalCov[(a * L) + b] += dCov;
            totalCov[(b * L) + a] += dCov;
            
        }
            
    }

}


// Update the summed covariance matrix

void updateCovarianceIntegrate(double dg, // time step
                               const std::vector<double> &p1_0, // single allele frequencies
                               const std::vector<double> &p2_0, // pair allele frequencies
                               const std::vector<double> &p1_1, // single allele frequencies
                               const std::vector<double> &p2_1, // pair allele frequencies
                               double totalCov[] // integrated covariance matrix
                               ) {

    int L = (int) p1_0.size();
                          
    // Iterate through states and add contributions to covariance matrix

    for (int a=0;a<L;a++) {
        
        totalCov[(a * L) + a] += dg * ( ((3 - (2 * p1_1[a])) * (p1_0[a] + p1_1[a])) - (2 * p1_0[a] * p1_0[a]) ) / 6;
        
        for (int b=a+1;b<L;b++) {
        
            double dCov1 = -dg * ((2 * p1_0[a] * p1_0[b]) + (2 * p1_1[a] * p1_1[b]) + (p1_0[a] * p1_1[b]) + (p1_1[a] * p1_0[b])) / 6;
            double dCov2 = dg * 0.5 * (p2_0[(a * L) + b] + p2_1[(a * L) + b]);

            totalCov[(a * L) + b] += dCov1 + dCov2;
            totalCov[(b * L) + a] += dCov1 + dCov2;
            
        }
        
    }

}


// Update the summed mutation vector (flux out minus flux in)

void updateMu(double dg, // time step
              const std::vector<double> &p1, // single allele frequencies
              std::vector<double> &totalMu // contribution to selection estimate from mutation
              ) {

    // Iterate through states and add contributions to mutation term
                  
    for (int a=0;a<p1.size();a++) totalMu[a] += dg * (p1[a] - (1 - p1[a]));

}


// Update the summed mutation vector (flux out minus flux in)
// Note: since first row of mutation matrix is the reference, the mutation matrix is SHIFTED wrt frequencies,
// because the reference frequency is not explicitly tracked

void updateMu(double dg, // time step
              const Vector &muMatrix, // mutation matrix
              const std::vector<double> &p1, // single allele frequencies
              std::vector<double> &totalMu // contribution to selection estimate from mutation
              ) {

    int q = (int) muMatrix.size(); // number of tracked alleles (states)
    int L = (int) p1.size() / (q-1); // number of sites in the sequence

    // Iterate through sites and add contributions to mutation term
    
    for (int i=0;i<L;i++) {
        
        double ptot = 0; // sum of frequencies of all explicitly tracked alleles at each site
        for (int a=0;a<(q-1);a++) ptot += p1[(i * (q-1)) + a];

        for (int a=0;a<(q-1);a++) {

            double fluxIn  = 0; // mutational flux from other alleles into the present one
            double fluxOut = 0; // mutational flux from the present alleles out to others

            // Iterate through other states to compute the mutational flux
            
            for (int b=0;b<a;b++) {

                fluxIn  += p1[(i * (q-1)) + b] * muMatrix[b+1][a+1];
                fluxOut += p1[(i * (q-1)) + a] * muMatrix[a+1][b+1];

            }
            for (int b=a+1;b<(q-1);b++) {

                fluxIn  += p1[(i * (q-1)) + b] * muMatrix[b+1][a+1];
                fluxOut += p1[(i * (q-1)) + a] * muMatrix[a+1][b+1];

            }
            
            // Add contribution from the implicit reference state
            
            fluxIn  +=          (1 - ptot) * muMatrix[0][a+1];
            fluxOut += p1[(i * (q-1)) + a] * muMatrix[a+1][0];
            
            // Get total contribution to the mutation term from fluxes

            totalMu[(i * (q-1)) + a] += dg * (fluxOut - fluxIn);

        }

    }

}


// Update the summed mutation vector (flux out minus flux in)
// Note: since first row of mutation matrix is the reference, the mutation matrix is SHIFTED wrt frequencies,
// because the reference frequency is not explicitly tracked

void updateMuIntegrate(double dg, // time step
                       const Vector &muMatrix, // mutation matrix
                       const std::vector<double> &p1_0, // single allele frequencies
                       const std::vector<double> &p1_1, // single allele frequencies
                       std::vector<double> &totalMu // contribution to selection estimate from mutation
                       ) {

    int q = (int) muMatrix.size(); // number of tracked alleles (states)
    int L = (int) p1_0.size() / (q-1); // number of sites in the sequence

    // Iterate through sites and add contributions to mutation term

    for (int i=0;i<L;i++) {

        double ptot = 0; // sum of frequencies of all explicitly tracked alleles at each site
        for (int a=0;a<(q-1);a++) ptot += 0.5 * (p1_0[(i * (q-1)) + a] + p1_1[(i * (q-1)) + a]);

        for (int a=0;a<(q-1);a++) {

            double fluxIn  = 0; // mutational flux from other alleles into the present one
            double fluxOut = 0; // mutational flux from the present alleles out to others

            // Iterate through other states to compute the mutational flux

            for (int b=0;b<a;b++) {

                fluxIn  += 0.5 * (p1_0[(i * (q-1)) + b] + p1_1[(i * (q-1)) + b]) * muMatrix[b+1][a+1];
                fluxOut += 0.5 * (p1_0[(i * (q-1)) + a] + p1_1[(i * (q-1)) + a]) * muMatrix[a+1][b+1];

            }
            for (int b=a+1;b<(q-1);b++) {

                fluxIn  += 0.5 * (p1_0[(i * (q-1)) + b] + p1_1[(i * (q-1)) + b]) * muMatrix[b+1][a+1];
                fluxOut += 0.5 * (p1_0[(i * (q-1)) + a] + p1_1[(i * (q-1)) + a]) * muMatrix[a+1][b+1];

            }

            // Add contribution from the implicit reference state

            fluxIn  +=                                            (1 - ptot) * muMatrix[0][a+1];
            fluxOut += 0.5 * (p1_0[(i * (q-1)) + a] + p1_1[(i * (q-1)) + a]) * muMatrix[a+1][0];

            // Get total contribution to the mutation term from fluxes

            totalMu[(i * (q-1)) + a] += dg * (fluxOut - fluxIn);

        }

    }

}


// Process asymptotic sequences (long time limit)

void processAsymptotic(const IntVVector &sequences, const Vector &counts, const Vector &muMatrix, int q, double totalCov[], double dx[]) {

    int    L     = ((int) sequences[0][0].size()) * (q-1);  // sequence length (i.e. number of tracked alleles)
    double cNorm = (double) counts.size();                  // total number of time points, needed to normalize overall frequencies
    
    IntVector cSequences;               // sequence vector collapsed in time
    std::vector<double> cCounts;        // count vector collapsed in time
    std::vector<double> totalMu(L,0);   // accumulated mutation term
    std::vector<double> p1(L,0);        // total allele frequency vector
    std::vector<double> p2(L*L,0);      // total allele pair frequencies
    
    // collapse sequence and count vectors
    
    for (int k=0;k<sequences.size();k++) {
    
        for (int i=0;i<sequences[k].size();i++) {
        
            cSequences.push_back(sequences[k][i]);
            cCounts.push_back(counts[k][i]/cNorm);
            
        }
        
    }
    
    // compute allele frequencies and covariance
    
    computeAlleleFrequencies(cSequences, cCounts, q, p1, p2);
    updateCovariance(1, p1, p2, totalCov);
    updateMu(1, muMatrix, p1, totalMu);
    
    // gather dx and totalMu terms
    
    for (int a=0;a<L;a++) dx[a] += totalMu[a];

}


// Process standard sequences (time series)

void processStandard(const IntVVector &sequences, // vector of sequence vectors
                     const Vector &counts, // vector of sequence counts
                     const std::vector<double> &times, // sequence sampling times
                     const Vector &muMatrix, // matrix of mutation rates
                     int q, // number of states (e.g., number of nucleotides or amino acids)
                     double totalCov[], // integrated covariance matrix
                     double dx[] // selection estimate numerator
                     ) {

    int L = ((int) sequences[0][0].size()) * (q-1); // sequence length (i.e. number of tracked alleles)
    std::vector<double> totalMu(L,0);               // accumulated mutation term
    std::vector<double> p1(L,0);                    // current allele frequency vector
    std::vector<double> p2(L*L,0);                  // current allele pair frequencies
    std::vector<double> lastp1(L,0);                // previous allele frequency vector
    std::vector<double> lastp2(L*L,0);              // previous allele pair frequencies
    std::vector<double> xp1(L,0);                   // processed (midpoint, etc) allele frequency vector
    std::vector<double> xp2(L*L,0);                 // processed (midpoint, etc) allele pair frequencies
    
    // set initial allele frequency and covariance then loop
    
    computeAlleleFrequencies(sequences[0], counts[0], q, lastp1, lastp2);
    for (int a=0;a<L;a++) dx[a] -= lastp1[a];
    
    //--- debug
    if (useDebug) {
        int lwidth = 5;
        printf("dx = {\t");
        for (int a=0;a<L;a++) { if (a%lwidth==0 && a>0) printf("\n\t"); printf("%.4e\t",dx[a]); }
        printf("}\n\n");
    }
    //---
    
    // Iterate through sets of sequences collected at all time points,
    // computing allele frequencies, interpolating, and
    // adding contributions to selection coefficient estimates
    
    //printf("t1\tt2\tt2-t1\titer\tnSteps\tdt\tweight\tdg\n");
    
    for (int k=1;k<sequences.size();k++) {
    
        computeAlleleFrequencies(sequences[k], counts[k], q, p1, p2);
//        interpolateFrequencies(lastp1, lastp2, p1, p2, 0.5, xp1, xp2);
//        updateCovariance(times[k] - times[k-1], xp1, xp2, totalCov);
//        updateMu(times[k] - times[k-1], muMatrix, xp1, totalMu);
        updateCovarianceIntegrate(times[k] - times[k-1], lastp1, lastp2, p1, p2, totalCov);
        updateMuIntegrate(times[k] - times[k-1], muMatrix, lastp1, p1, totalMu);
        
        if (k==sequences.size()-1) { for (int a=0;a<L;a++) dx[a] += p1[a]; }
        else { lastp1 = p1; lastp2 = p2; }
    
    }
    
    //--- debug
    if (useDebug) {
        int lwidth = 5;
        printf("dx = {\t");
        for (int a=0;a<L;a++) { if (a%lwidth==0 && a>0) printf("\n\t"); printf("%.4e\t",dx[a]); }
        printf("}\n\n");
    }
    //---
    
    // Gather dx and totalMu terms
    
    for (int a=0;a<L;a++) dx[a] += totalMu[a];

}


// Add Gaussian regularization for selection coefficients (modifies integrated covariance)

void regularizeCovariance(const IntVVector &sequences, // vector of sequence vectors
                          int q, // number of states (e.g., number of nucleotides or amino acids)
                          double gammaN, // normalized regularization strength
                          bool useCovariance, // if false, don't consider off-diagonal terms
                          double totalCov[] // integrated covariance matrix
                          ) {

    int L = ((int) sequences[0][0].size()) * (q-1);

    for (int a=0;a<L;a++) totalCov[(a * L) + a] += 2 * gammaN; // standard regularization + diagonal part of gauged state regularization

    // regularization including gauged state

    if (useCovariance) { for (int i=0;i<sequences[0][0].size();i++) {

            int ix = i * (q-1);

            for (int a=0;a<q-1;a++) { for (int b=a+1;b<q-1;b++) {

                totalCov[((ix + a) * L) + ix + b] += gammaN;
                totalCov[((ix + b) * L) + ix + a] += gammaN;

            } }

    } }

}


// MAIN PROGRAM

int run(RunParameters &r) {
    
    // READ IN SEQUENCES FROM DATA
    
    IntVVector sequences;       // set of integer sequences at each time point
    Vector counts;              // counts for each sequence at each time point
    std::vector<double> times;  // list of times of measurements
    
    if (FILE *datain = fopen(r.getSequenceInfile().c_str(),"r")) { getSequences(datain, sequences, counts, times); fclose(datain); }
    else { printf("Problem retrieving data from file %s! File may not exist or cannot be opened.\n",r.getSequenceInfile().c_str()); return EXIT_FAILURE; }
    
    if (r.useVerbose) {
        if (times.size()>1) printf("Got %d time points, %.2f difference between first two times.\n",(int)sequences.size(),times[1]-times[0]);
        printf("Total number of loci (sites) is %d, number of alleles (states) is %d.\n",(int)sequences[0][0].size(),(int)r.q);
        if (counts.size()>1) {
            double sum = 0;
            for (int i=0;i<counts[1].size();i++) sum += counts[1][i];
            printf("Sum of all counts at the second time point is %.2f (1.00 expected).\n\n",sum);
        }
        printf("Parameters: N = %.2e, mu = %.2e, gamma = %.2e.\n\n",r.N,r.mu,r.gamma);
        if (counts.size()>0) {
            printf("First sequence at the first time point:\n%.2f\t",counts[0][0]);
            for (int a=0;a<sequences[0][0].size();a++) printf("%d ",sequences[0][0][a]);
            printf("\n\n");
        }
    }
    
    Vector muMatrix;    // matrix of mutation rates
    
    if (r.useMatrix) {
    
        if (FILE *muin = fopen(r.getMuInfile().c_str(),"r")) { getMu(muin, muMatrix); fclose(muin); }
        else { printf("Problem retrieving data from file %s! File may not exist or cannot be opened.\n",r.getMuInfile().c_str()); return EXIT_FAILURE; }
        
        r.q = (int) muMatrix.size();
        
    }
    else {
        
        muMatrix.resize(r.q, std::vector<double>(r.q, r.mu));
        for (int i=0;i<r.q;i++) muMatrix[i][i] = 0;
        
    }
    
    // PROCESS SEQUENCES
    
    int    L         = (int) sequences[0][0].size() * (r.q-1);  // sequence length (i.e. number of tracked alleles)
    double tol       = r.tol;                                   // tolerance for changes in covariance between time points
    double gammaN    = r.gamma/r.N;                             // regularization strength divided by population size
    double *dx       = new double[L];                           // difference between start and end allele frequencies
    double *totalCov = new double[L*L];                         // accumulated allele covariance matrix
    
    for (int a=0;a<  L;a++) dx[a]       = 0;
    for (int a=0;a<L*L;a++) totalCov[a] = 0;
    
    // _ START TIMER
    auto t_start = Clock::now();
    
    if (r.useAsymptotic) processAsymptotic(sequences, counts, muMatrix, r.q, totalCov, dx);
    else                 processStandard(sequences, counts, times, muMatrix, r.q, totalCov, dx);
    
    // If there is more than one input trajectory, loop through all of them and add contributions
    // NOTE: CURRENT CODE ASSUMES THAT ALL VALUES OF N ARE EQUAL
    
    if (r.infiles.size()>1) { for (int k=1;k<r.infiles.size();k++) {
    
        // Reset trajectory variables and reload them with new data
        
        sequences.clear();
        counts.clear();
        times.clear();
        
        if (FILE *datain = fopen(r.getSequenceInfile(k).c_str(),"r")) { getSequences(datain, sequences, counts, times); fclose(datain); }
        else { printf("Problem retrieving data from file %s! File may not exist or cannot be opened.\n",r.getSequenceInfile().c_str()); return EXIT_FAILURE; }
        
        // Add contributions to dx and totalCov
        
        if (r.useAsymptotic) processAsymptotic(sequences, counts, muMatrix, r.q, totalCov, dx);
        else                 processStandard(sequences, counts, times, muMatrix, r.q, totalCov, dx);
    
    } }
    
    // REGULARIZE
    
    regularizeCovariance(sequences, r.q, gammaN, r.useCovariance, totalCov);
    
    // RECORD COVARIANCE (optional)
    
    if (r.saveCovariance) {
        if (FILE *dataout = fopen(r.getCovarianceOutfile().c_str(),"w")) { printCovariance(dataout, totalCov, L); fclose(dataout); }
        else { printf("Problem writing data to file %s! File may not be created or cannot be opened.\n",r.getCovarianceOutfile().c_str()); return EXIT_FAILURE; }
    }
    
    // RECORD NUMERATOR (optional)
    
    if (r.saveNumerator) {
        if (FILE *dataout = fopen(r.getNumeratorOutfile().c_str(),"w")) { printNumerator(dataout, dx, L); fclose(dataout); }
        else { printf("Problem writing data to file %s! File may not be created or cannot be opened.\n",r.getCovarianceOutfile().c_str()); return EXIT_FAILURE; }
    }
    
    // INFER THE SELECTION COEFFICIENTS -- solve Cov . sMAP = dx
    
    std::vector<double> sMAP(L,0);
    
    if (r.useCovariance) {
    
        int status;
    
        gsl_matrix_view _cov = gsl_matrix_view_array(totalCov, L, L);   // gsl covariance + Gaussian regularization
        gsl_vector_view  _dx = gsl_vector_view_array(dx, L);            // gsl dx vector
        gsl_vector    *_sMAP = gsl_vector_alloc(L);                     // maximum a posteriori selection coefficients for each allele
        gsl_permutation  *_p = gsl_permutation_alloc(L);
        
        gsl_linalg_LU_decomp(&_cov.matrix, _p, &status);
        gsl_linalg_LU_solve(&_cov.matrix, _p, &_dx.vector, _sMAP);
        
        for (int a=0;a<L;a++) sMAP[a] = gsl_vector_get(_sMAP, a);
        
        gsl_permutation_free(_p);
        gsl_vector_free(_sMAP);
        
        delete[] dx;
        delete[] totalCov;
        
    }
    
    else {
    
        for (int a=0;a<L;a++) sMAP[a] = dx[a] / totalCov[(a * L) + a];
    
    }
    
    auto t_end = Clock::now();
    // ^ END TIMER
    
    printf("%lld\n",std::chrono::duration_cast<std::chrono::nanoseconds>(t_end - t_start).count());
    
    // WRITE TO FILE
    
    if (FILE *dataout = fopen(r.getSelectionCoefficientOutfile().c_str(),"w")) { printSelectionCoefficients(dataout, sMAP); fclose(dataout); }
    else { printf("Problem writing data to file %s! File may not be created or cannot be opened.\n",r.getSelectionCoefficientOutfile().c_str()); return EXIT_FAILURE; }
    
    if (r.useVerbose) {
        int lwidth = 5;
        printf("s = {\t");
        for (int a=0;a<L;a++) { if (a%lwidth==0 && a>0) printf("\n\t"); printf("%.4e\t",sMAP[a]); }
        printf("}\n");
    }
    
    return EXIT_SUCCESS;
 
}
