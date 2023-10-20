#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <x86intrin.h>

#include <igraph.h>

#include "slater_condon.h"

int get_matelem(size_t deti, size_t detj) ;

void printBits(int num, int len) ;

void generateConfigurations(int norb, int nelec, int* configAll, int* size) ;

// Function to find the positions of a list of configurations in a sorted list
void findPositions(int* configList, int sizeList, int* configs, int sizeConfigs, int* positions) ;

// Function to compare two configurations for qsort and bsearch
int compare(const void* a, const void* b) ;

// Function to calculate binomial coefficient using lgamma function
long long binomialCoeff(int n, int k) ;

void printPositions(int* positions, int size) ;

// Function to generate all possible alpha determinants
void generateAlphaDeterminants(int* configAlpha, int sizeAlpha, const igraph_t* graph, size_t alphaConfig, igraph_vector_t* alphaDeterminants) ;

int* igraphVectorToIntArray(const igraph_vector_t* igraph_vector) ;
