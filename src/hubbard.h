#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <x86intrin.h>

#include <igraph.h>

#include "slater_condon.h"

#define MAX_Sz_BLOCKS 100

// Function to find a global ID in the hash table given an alpha and beta ID
size_t findGlobalID(size_t alphaID, size_t betaID, size_t nalpha) ;

// Functions to find alpha and beta IDs from globalID
size_t findAlphaID(size_t globalID, size_t nalpha, size_t nbeta) ;
size_t findBetaID(size_t globalID, size_t nalpha, size_t nbeta) ;

double solveQuad(double a, double b, double c) ;

size_t get_matelem(size_t deti, size_t detj) ;

void printBits(size_t num, size_t len) ;

void generateConfigurations(size_t norb, size_t nelec, size_t* configAll, size_t* size) ;

// Function to find the positions of a list of configurations in a sorted list
void findPositions(size_t* configList, size_t sizeList, size_t* configs, size_t sizeConfigs, size_t* positions) ;

// Function to compare two configurations for qsort and bsearch
int compare(const void* a, const void* b) ;

// Function to calculate binomial coefficient using lgamma function
long long binomialCoeff(size_t n, size_t k) ;

// A function to print a matrix
void print_matrix_d(double** matrix, int rows, int cols) ;

void printPositions(size_t* positions, size_t size) ;

int getPhase(size_t alphaConfig, size_t newAlphaConfig, size_t h, size_t p) ;

int getExecDegree(size_t detI, size_t detJ) ;

int getHoles_1ex(size_t detI, size_t detJ, size_t *holesOut) ;
int getPart_1ex(size_t detI, size_t detJ, size_t *particlesOut) ;

// Function to generate all possible alpha determinants
void generateDeterminants(size_t* configAlpha, size_t sizeAlpha, const igraph_t* graph, size_t alphaConfig, size_t betaConfig, igraph_vector_t* alphaDeterminants, igraph_vector_t* alphaMEs, int alphaBeta, double** wmat) ;

// Function to generate all possible alpha determinants given a list of alpha determinants
void generateAllDeterminants(size_t *configAlpha, size_t sizeAlpha, const igraph_t* graph, size_t* alphaConfigs, size_t numConfigs, igraph_vector_t* allAlphaDeterminants, double** wmat) ;

size_t* igraphVectorToIntArray(const igraph_vector_t* igraph_vector) ;

// Main function that calculates MEs
void getAllHubbardMEs(size_t Idet, igraph_vector_t* MElist, igraph_vector_t* Jdetlist, size_t *configAlpha, size_t sizeAlpha, size_t *configBeta, size_t sizeBeta, const igraph_t* graph, double** wmat) ;

// Get the diagonal part of the hubbard Hamiltonian
int getHubbardDiag(size_t Idet, size_t *configAlpha, size_t sizeAlpha, size_t *configBeta, size_t sizeBeta) ;

// Main function that generates S2 operator
void getS2Operator(size_t Idet, igraph_vector_t* MElist, igraph_vector_t* Jdetlist, size_t *configAlpha, size_t sizeAlpha, size_t *configBeta, size_t sizeBeta, const igraph_t* graph, int natom, int natomax) ;

// Get maximum neighbors
int getMaxNeighbors(const igraph_t* graph, size_t nsites) ;

// A function to declare a matrix of given size and initialize it to 0
double** declare_matrix(int rows, int cols) ;

// A function to save a matrix in a file in CSV format
void save_matrix(double** matrix, int rows, int cols, char* filename) ;

// A function to find the number of permutations and the phase for a given pair of bit strings
int calculate (int a, int b) ;

// Function to get the Sz  operator
void getSzOperator(size_t detIa, size_t detIb, double *szval, size_t* configAlpha, size_t sizeAlpha, size_t* configBeta, size_t sizeBeta, int nblk, size_t* SzBlock) ;

// Function to get the num  operator
void getNumOperator(size_t detIa, size_t detIb, double *numval, double *numvala, size_t* configAlpha, size_t sizeAlpha, size_t* configBeta, size_t sizeBeta, int nblk, size_t* NumBlock) ;
