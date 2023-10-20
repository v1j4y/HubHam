#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <x86intrin.h>
#include "uthash.h"

#include <igraph.h>

#include "slater_condon.h"

// Structure for the combined ID
struct CombinedID {
    size_t alphaID;
    size_t betaID;
};

// Structure for the hash table
struct IDMap {
    struct CombinedID combinedID;  // key
    unsigned long long ID;
    UT_hash_handle hh;  // makes this structure hashable
};

// Function to add a pair of IDs to the hash table
void addID(struct IDMap** idMap, struct CombinedID combinedID, size_t globalID) ;

// Function to find a global ID in the hash table given an alpha and beta ID
unsigned long long findGlobalID(struct IDMap** idMap, size_t alphaID, size_t betaID) ;

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

// Function to generate all possible alpha determinants given a list of alpha determinants
void generateAllAlphaDeterminants(int *configAlpha, int sizeAlpha, const igraph_t* graph, size_t* alphaConfigs, int numConfigs, igraph_vector_t* allAlphaDeterminants) ;

int* igraphVectorToIntArray(const igraph_vector_t* igraph_vector) ;
