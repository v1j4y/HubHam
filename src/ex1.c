/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Standard symmetric eigenproblem corresponding to the Laplacian operator in 1 dimension.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions = matrix dimension.\n\n";

#include <slepceps.h>

#include "hubbard.h"
#include "readgraphmllib.h"

int main(int argc,char **argv)
{
  Mat            A;           /* problem matrix */
  EPS            eps;         /* eigenproblem solver context */
  EPSType        type;
  PetscReal      error,tol,re,im;
  PetscScalar    kr,ki;
  Vec            xr,xi;

  // Assume configAlpha and configBeta are sorted lists of all possible alpha and beta configurations
  size_t norb = 6;
  size_t nalpha = 3;
  size_t nbeta = 3;

  size_t sizeAlpha = binomialCoeff(norb, nalpha);
  size_t sizeBeta = binomialCoeff(norb, nbeta);

  size_t* configAlpha = malloc(sizeAlpha * sizeof(size_t));
  size_t* configBeta = malloc(sizeBeta * sizeof(size_t));

  generateConfigurations(norb, nalpha, configAlpha, &sizeAlpha);
  generateConfigurations(norb, nbeta, configBeta, &sizeBeta);
  printf("Number of configurations (alpha = %ld) beta = %ld:\n",sizeAlpha,sizeBeta);

  // Sort the lists for binary search
  qsort(configAlpha, sizeAlpha, sizeof(size_t), compare);
  qsort(configBeta, sizeBeta, sizeof(size_t), compare);

  // Find the positions of a list of specific configurations
  size_t alphaConfigs[] = {0b000111, 0b001110};  // Example alpha configurations
  size_t betaConfigs[] = {0b010011, 0b001101};   // Example beta configurations

  size_t sizeAlphaConfigs = sizeof(alphaConfigs) / sizeof(alphaConfigs[0]);
  size_t sizeBetaConfigs = sizeof(betaConfigs) / sizeof(betaConfigs[0]);

  size_t* posAlpha = malloc(sizeAlphaConfigs * sizeof(size_t));
  size_t* posBeta = malloc(sizeBetaConfigs * sizeof(size_t));

  findPositions(configAlpha, sizeAlpha, alphaConfigs, sizeAlphaConfigs, posAlpha);
  findPositions(configBeta, sizeBeta, betaConfigs, sizeBetaConfigs, posBeta);

  printf("Positions of alpha configurations:\n");
  printPositions(posAlpha, sizeAlphaConfigs);

  printf("\nPositions of beta configurations:\n");
  printPositions(posBeta, sizeBetaConfigs);

  const char* graphmlFileName = "/home/chilkuri/Documents/codes/c_codes/hubbard_slepc/data/graphm3.graphml";
  FILE* graphmlFile = fopen(graphmlFileName, "r");

  if (graphmlFile == NULL) {
    fprintf(stderr, "Error opening the file.\n");
    return 1;  // Return an error code or use another error handling method.
  }

  igraph_t graph;
  igraph_empty(&graph, 0, IGRAPH_DIRECTED);

  if (readGraphMLFile(graphmlFile, &graph)) {
    // Successfully read the graph, now you can work with 'graph'.
    igraph_integer_t num_vertices = igraph_vcount(&graph);
    printf("Number of vertices: %ld\n", (long)num_vertices);
  }

  // Example alpha configuration
  size_t alphaConfig = 0b001101;
  size_t betaConfig  = 0b010101;

  // Generate all possible alpha determinants
  igraph_vector_t alphaDeterminants;
  igraph_vector_init(&alphaDeterminants, 0);
  igraph_vector_t alphaMEs;
  igraph_vector_init(&alphaMEs, 0);
  generateDeterminants(configAlpha, sizeAlpha, &graph, alphaConfig, &alphaDeterminants, &alphaMEs);

  // Print the generated alpha determinants
  printf("\nNumber of generated alpha determinants:%ld\n",igraph_vector_size(&alphaDeterminants));
  size_t *int_alphaDeterminants = igraphVectorToIntArray(&alphaDeterminants);
  for (size_t i = 0; i < igraph_vector_size(&alphaDeterminants); ++i) {
    printBits(configAlpha[int_alphaDeterminants[i]], norb);
    unsigned long long foundGlobalID = findGlobalID(configAlpha[int_alphaDeterminants[i]], betaConfig, sizeAlpha);
    printf("%ld - %f (id = %llu)\n", i, VECTOR(alphaDeterminants)[i], foundGlobalID);
  }

  // Now you can find the global ID in the hash table given an alpha and beta ID
  // For example:
  size_t alphaID;
  size_t betaID;
  findPositions(configAlpha, sizeAlpha, &configAlpha[10], 1, &alphaID);
  findPositions(configBeta, sizeBeta, &configBeta[14], 1, &betaID);

  size_t foundGlobalID = findGlobalID(alphaID, betaID, sizeAlpha);

  if (foundGlobalID != 0) {
    printf("(%ld, %ld) Global ID: %ld\n", alphaID, betaID, foundGlobalID);
  }
  printf(" alphaID = %ld betaID = %ld\n",findAlphaID(foundGlobalID,sizeAlpha,sizeBeta),findBetaID(foundGlobalID,sizeAlpha,sizeBeta));
  printf(" Phase = %d Phase = %d\n",getPhase(14,7,1,3), getPhase(7,14,3,1));

  // Declare a matrix of size 3 x 4
  int rows = 400;
  int cols = 400;
  int** matrix = declare_matrix(rows, cols);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscInt       n=sizeAlpha*sizeBeta,i,Istart,Iend,nev,maxit,its,nconv;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,(char*)0,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n1-D Laplacian Eigenproblem, n=%" PetscInt_FMT "\n\n",n));
  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetFromOptions(A));
  PetscCall(MatSetUp(A));

  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
  printf(" Istart=%ld Iend=%ld\n",Istart,Iend);
  for (i=Istart;i<Iend;i++) {

    igraph_vector_t MElist;
    igraph_vector_init(&MElist, 0);
    igraph_vector_t Jdetlist;
    igraph_vector_init(&Jdetlist, 0);

    //if (i>0) PetscCall(MatSetValue(A,i,i-1,-1.0,INSERT_VALUES));
    //if (i<n-1) PetscCall(MatSetValue(A,i,i+1,-1.0,INSERT_VALUES));
    //PetscCall(MatSetValue(A,i,i,2.0,INSERT_VALUES));
    getAllHubbardMEs(i, &MElist, &Jdetlist, configAlpha, sizeAlpha, configBeta, sizeBeta, &graph);
    for (int j = 0; j < igraph_vector_size(&Jdetlist); ++j) {
      int Jid = VECTOR(Jdetlist)[j];
      matrix[i][Jid] = VECTOR(MElist)[j];
      PetscCall(MatSetValue(A,i,Jid,(double)VECTOR(MElist)[j],INSERT_VALUES));
      //printf(" i=%d j=%d \n",i,Jid);
    }

    igraph_vector_destroy(&MElist);
    igraph_vector_destroy(&Jdetlist);
  }
  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));

  PetscCall(MatCreateVecs(A,NULL,&xr));
  PetscCall(MatCreateVecs(A,NULL,&xi));

  // Save file
  save_matrix(matrix, rows, cols, "/tmp/benzene_c.csv");

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Create eigensolver context
  */
  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));

  /*
     Set operators. In this case, it is a standard eigenvalue problem
  */
  PetscCall(EPSSetOperators(eps,A,NULL));
  PetscCall(EPSSetProblemType(eps,EPS_HEP));

  /*
     Set solver parameters at runtime
  */
  PetscCall(EPSSetFromOptions(eps));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(EPSSolve(eps));
  /*
     Optional: Get some information from the solver and display it
  */
  PetscCall(EPSGetIterationNumber(eps,&its));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %" PetscInt_FMT "\n",its));
  PetscCall(EPSGetType(eps,&type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type));
  PetscCall(EPSGetDimensions(eps,&nev,NULL,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %" PetscInt_FMT "\n",nev));
  PetscCall(EPSGetTolerances(eps,&tol,&maxit));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%" PetscInt_FMT "\n",(double)tol,maxit));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Get number of converged approximate eigenpairs
  */
  PetscCall(EPSGetConverged(eps,&nconv));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %" PetscInt_FMT "\n\n",nconv));

  if (nconv>0) {
    /*
       Display eigenvalues and relative errors
    */
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
         "           k          ||Ax-kx||/||kx||\n"
         "   ----------------- ------------------\n"));

    for (i=0;i<nconv;i++) {
      /*
        Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and
        ki (imaginary part)
      */
      PetscCall(EPSGetEigenpair(eps,i,&kr,&ki,xr,xi));
      /*
         Compute the relative error associated to each eigenpair
      */
      PetscCall(EPSComputeError(eps,i,EPS_ERROR_RELATIVE,&error));

#if defined(PETSC_USE_COMPLEX)
      re = PetscRealPart(kr);
      im = PetscImaginaryPart(kr);
#else
      re = kr;
      im = ki;
#endif
      if (im!=0.0) PetscCall(PetscPrintf(PETSC_COMM_WORLD," %9f%+9fi %12g\n",(double)re,(double)im,(double)error));
      else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"   %12f       %12g\n",(double)re,(double)error));
    }
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n"));
  }

  igraph_vector_destroy(&alphaDeterminants);
  igraph_vector_destroy(&alphaMEs);

  // Don't forget to destroy the graph when you're done.
  igraph_destroy(&graph);

  // Close the file when you're done with it.
  fclose(graphmlFile);

  free(configAlpha);
  free(configBeta);

  // Free the memory allocated for the matrix
  for (int i = 0; i < rows; i++) {
    free(matrix[i]);
  }
  free(matrix);

  /*
     Free work space
  */
  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&A));
  PetscCall(VecDestroy(&xr));
  PetscCall(VecDestroy(&xi));
  PetscCall(SlepcFinalize());
  return 0;
}
