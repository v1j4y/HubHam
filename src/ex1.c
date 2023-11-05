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
#include "get_s2.h"

int main(int argc,char **argv)
{
  Mat            A;           /* problem matrix */
  EPS            eps;         /* eigenproblem solver context */
  EPSType        type;
  PetscReal      error,tol,re,im;
  PetscScalar    kr,ki;
  Vec            xr,xi;
  PetscLogDouble t1,t2,tt1,tt2;
  PetscReal normfin;
  PetscReal xymatfin = 0.0;

  const char* graphmlFileName = "/home/chilkuri/Documents/codes/c_codes/hubbard_slepc/data/c4h6.graphml";
  FILE* graphmlFile = fopen(graphmlFileName, "r");

  if (graphmlFile == NULL) {
    fprintf(stderr, "Error opening the file.\n");
    return 1;  // Return an error code or use another error handling method.
  }

  igraph_t graph;
  igraph_empty(&graph, 0, IGRAPH_DIRECTED);
  igraph_integer_t num_vertices;

  if (readGraphMLFile(graphmlFile, &graph)) {
    // Successfully read the graph, now you can work with 'graph'.
    num_vertices = igraph_vcount(&graph);
  }

  // Assume configAlpha and configBeta are sorted lists of all possible alpha and beta configurations
  size_t norb = num_vertices;
  size_t nalpha = norb/2;
  size_t nbeta = norb/2;

  size_t sizeAlpha = binomialCoeff(norb, nalpha);
  size_t sizeBeta = binomialCoeff(norb, nbeta);

  size_t* configAlpha = malloc(sizeAlpha * sizeof(size_t));
  size_t* configBeta = malloc(sizeBeta * sizeof(size_t));

  generateConfigurations(norb, nalpha, configAlpha, &sizeAlpha);
  generateConfigurations(norb, nbeta, configBeta, &sizeBeta);

  // Sort the lists for binary search
  qsort(configAlpha, sizeAlpha, sizeof(size_t), compare);
  qsort(configBeta, sizeBeta, sizeof(size_t), compare);

  // Declare a matrix of size 3 x 4
  int rows = sizeAlpha * sizeBeta;
  int cols = rows;
  int** matrix = declare_matrix(rows, cols);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Define Hamiltonian
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  double U =  1.0;
  double t = -1.0;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscInt       n=sizeAlpha*sizeBeta,i,Istart,Iend,nev,maxit,its,nconv;
  PetscInt		   ncv, mpd;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,(char*)0,help));

  PetscCall(PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n1-D Laplacian Eigenproblem, n=%" PetscInt_FMT "\n\n",n));
  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,n,n,3*num_vertices,NULL,3*num_vertices,NULL,&A));
  PetscCall(MatMPIAIJSetPreallocation(A,3*num_vertices,NULL,3*num_vertices,NULL));
  //PetscCall(MatSetFromOptions(A));
  //PetscCall(MatSetUp(A));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n Number of vertices: %ld",(long)num_vertices));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n Number of configurations alpha = %ld, beta = %ld\n",sizeAlpha,sizeBeta));
  PetscCall(PetscTime(&tt1));

  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {

    igraph_vector_t MElist;
    igraph_vector_init(&MElist, 0);
    igraph_vector_t Jdetlist;
    igraph_vector_init(&Jdetlist, 0);

    int diag = getHubbardDiag(i, configAlpha, sizeAlpha, configBeta, sizeBeta);
    PetscCall(MatSetValue(A,i,i,(double)diag,INSERT_VALUES));
    matrix[i][i] = (double)diag*U;
    getAllHubbardMEs(i, &MElist, &Jdetlist, configAlpha, sizeAlpha, configBeta, sizeBeta, &graph);
    for (int j = 0; j < igraph_vector_size(&Jdetlist); ++j) {
      int Jid = VECTOR(Jdetlist)[j];
      matrix[i][Jid] = t*VECTOR(MElist)[j];
      PetscCall(MatSetValue(A,i,Jid,t*(double)VECTOR(MElist)[j],INSERT_VALUES));
    }

    igraph_vector_destroy(&MElist);
    igraph_vector_destroy(&Jdetlist);
  }
  PetscCall(PetscTime(&tt2));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Time used to build the matrix: %f\n",tt2-tt1));


  PetscCall(PetscTime(&tt1));
  PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
  PetscCall(PetscTime(&tt2));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Time used to assemble the matrix: %f\n",tt2-tt1));

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
  PetscCall(EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL));

  /*
     Set solver parameters at runtime
  */
  PetscCall(EPSSetFromOptions(eps));
  tol = 1.e-9;
  maxit = 10000000;
  PetscCall(EPSSetTolerances(eps,tol,maxit));
  ncv  = 9;
  mpd  = 10;
  nev  = 80;
  PetscCall(EPSSetDimensions(eps,nev,PETSC_DECIDE,PETSC_DECIDE));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(PetscTime(&t1));
  PetscCall(EPSSolve(eps));
  PetscCall(PetscTime(&t2));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Time used to Solve EVP: %f\n",t2-t1));

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

      /*
       * Get eigenvector
       */
      Vec vi, vr;
      BV V;
      PetscReal norm;
      PetscReal *values;
      int sbx = 0, pos1 = 1, pos2 = 4;
      PetscReal xymat = 0.0;
      PetscReal xymat2 = 0.0;
      PetscReal xymat3 = 0.0;
      PetscReal xymat4 = 0.0;
      PetscReal weight3 = 0.0;
      PetscReal norm2 = 0.0;
      PetscReal norm3 = 0.0;
      PetscReal norm4 = 0.0;
      EPSGetBV(eps,&V);
      BVGetColumn(V,i,&vi);
      VecNorm(vi,NORM_2,&norm);
      VecGetArray(vi, &values);
      norm = 0.0;
      get_s2(values, &Istart, &Iend, values, &num_vertices, &norm, &norm2, &norm3, &norm4, &xymat, &xymat2, &xymat3, &xymat4, &weight3,
             &sbx, &sbx, &sbx, &sbx, &sbx, &sbx,
             &sbx, &sbx, &sbx, &sbx,
             &sbx, &sbx, &pos1, &pos2, &pos1, num_vertices, configAlpha, sizeAlpha, configBeta, sizeBeta);
      MPI_Reduce(&xymat, &xymatfin, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
      MPI_Reduce(&norm, &normfin, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,"norm = %10.5f xymat = %10.5f\n",normfin,xymatfin));
      VecRestoreArray(vi, &values);
      BVRestoreColumn(V,i,&vi);
    }
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n"));
  }
  //printf(" %d \n",calculate(7, 7));
  //int ideter[num_vertices], addr;
  //getdet(1, ideter, configAlpha, sizeAlpha, configBeta, sizeBeta, num_vertices);
  //printf("\n");
  //for( int i=0;i<num_vertices;++i ) {
  //  printf(" %d ",ideter[i]);
  //}
  //adr (ideter, &addr, configAlpha, sizeAlpha, configBeta, sizeBeta, num_vertices);
  //printf("\n\n %d \n",addr);

  //igraph_vector_destroy(&alphaDeterminants);
  //igraph_vector_destroy(&alphaMEs);

  // Don't forget to destroy the graph when you're done.
  igraph_destroy(&graph);

  // Close the file when you're done with it.
  fclose(graphmlFile);

  free(configAlpha);
  free(configBeta);

  // Free the memory allocated for the matrix
  //for (int i = 0; i < rows; i++) {
  //  free(matrix[i]);
  //}
  //free(matrix);

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
