#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <x86intrin.h>

#include <igraph.h>

#include "slater_condon.h"

// Function to find a global ID in the hash table given an alpha and beta ID
size_t findGlobalID(size_t alphaID, size_t betaID, size_t nalpha) ;

size_t get_matelem(size_t deti, size_t detj) ;

void printBits(size_t num, size_t len) ;

void generateConfigurations(size_t norb, size_t nelec, size_t* configAll, size_t* size) ;

// Function to find the positions of a list of configurations in a sorted list
void findPositions(size_t* configList, size_t sizeList, size_t* configs, size_t sizeConfigs, size_t* positions) ;

// Function to compare two configurations for qsort and bsearch
int compare(const void* a, const void* b) ;

// Function to calculate binomial coefficient using lgamma function
long long binomialCoeff(size_t n, size_t k) ;

void printPositions(size_t* positions, size_t size) ;

// Function to generate all possible alpha determinants
void generateDeterminants(size_t* configAlpha, size_t sizeAlpha, const igraph_t* graph, size_t alphaConfig, igraph_vector_t* alphaDeterminants) ;

// Function to generate all possible alpha determinants given a list of alpha determinants
void generateAllDeterminants(size_t *configAlpha, size_t sizeAlpha, const igraph_t* graph, size_t* alphaConfigs, size_t numConfigs, igraph_vector_t* allAlphaDeterminants) ;

size_t* igraphVectorToIntArray(const igraph_vector_t* igraph_vector) ;
