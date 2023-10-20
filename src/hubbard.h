#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <x86intrin.h>

#include "slater_condon.h"

int get_matelem(size_t deti, size_t detj) ;
void generateConfigurations(int norb, int nelec, int* configAll, int* size) ;

// Function to find the positions of a list of configurations in a sorted list
void findPositions(int* configList, int sizeList, int* configs, int sizeConfigs, int* positions) ;

// Function to compare two configurations for qsort and bsearch
int compare(const void* a, const void* b) ;

// Function to calculate binomial coefficient using lgamma function
long long binomialCoeff(int n, int k) ;

void printPositions(int* positions, int size) ;
