/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Hubbard Hamiltonian - A high-performance Exact diagonalization program

   Copyright (c) 2023-, Vijay Gopal CHILKURI
   LICENCE : GNU GPL V2
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Standard Hubbard eigenproblem corresponding to the topology given in the file.\n\n"
  "The command line options are:\n"
  "  -f <file>, where <file> = name of the graphml file.\n\n";

#include <slepceps.h>

#include "hubbard.h"
#include "readgraphmllib.h"
#include "utils.h"

int main(int argc,char **argv)
{
  Mat            A;           /* problem matrix */
  Mat            S2;          /* S2 oper matrix */
  EPS            eps;         /* eigenproblem solver context */
  EPSType        type;
  PetscReal      error,tol,re,im;
  PetscScalar    kr,ki, dot, spin;
  Vec            xr,xi, vs2;
  PetscLogDouble t1,t2,tt1,tt2;
  PetscReal normfin;
  PetscReal xymatfin = 0.0;

  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&argc,&argv,(char*)0,help));

  //const char* graphmlFileName = "/home/chilkuri/Documents/codes/c_codes/hubbard_slepc/data/graphm3.graphml";
  char        graphmlFileName[PETSC_MAX_PATH_LEN]; /* input file name */
  char        szblk[PETSC_MAX_PATH_LEN]; /* input file name */
  PetscBool       flg;
  PetscCall(PetscOptionsGetString(NULL, NULL, "-f", graphmlFileName, sizeof(graphmlFileName), &flg));
  PetscCheck(flg, PETSC_COMM_WORLD, PETSC_ERR_USER, "Must indicate graphml file with the -f option");
  FILE* graphmlFile = fopen(graphmlFileName, "r");

  if (graphmlFile == NULL) {
    fprintf(stderr, "Error opening the file.\n");
    return 1;  // Return an error code or use another error handling method.
  }

  /*
    Enter the options for the number of holes and
    the total number of alpha electrons.
  */
  PetscInt  nelec;
  PetscInt  nalpha;
  PetscInt  DoS2 = 0;
  PetscInt  hasW = 0;
  PetscInt  hasV = 0;
  PetscInt  DoSz = 0;
  PetscInt  DoNum = 0;
  PetscInt  DBGPrinting = 0;
  PetscReal t_inp = -1.0;
  PetscReal Ut_inp = 10.0;
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-ht",&t_inp,&flg));
  PetscCheck(flg, PETSC_COMM_WORLD, PETSC_ERR_USER, "Must indicate t with the -ht option");
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-hUt",&Ut_inp,&flg));
  PetscCheck(flg, PETSC_COMM_WORLD, PETSC_ERR_USER, "Must indicate U/t with the -hUt option");
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-hs2",&DoS2,&flg));
  PetscCheck(flg, PETSC_COMM_WORLD, PETSC_ERR_USER, "Must indicate whether or not to S2 with the -hs2 option (1=true, 0=false)");
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-pd",&DBGPrinting,&flg));
  PetscCheck(flg, PETSC_COMM_WORLD, PETSC_ERR_USER, "Must indicate debug printing with the -pd option");
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-ne",&nelec,&flg));
  PetscCheck(flg, PETSC_COMM_WORLD, PETSC_ERR_USER, "Must indicate the total number of e- with the -ne option");
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-na",&nalpha,&flg));
  PetscCheck(flg, PETSC_COMM_WORLD, PETSC_ERR_USER, "Must indicate the number of alpha e- with the -na option");
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-hsz",&DoSz,NULL));
  PetscCheck(flg, PETSC_COMM_WORLD, PETSC_ERR_USER, "Must indicate whether or not to Sz with the -hSz option (1=true, 0=false)");
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-hnum",&DoNum,NULL));
  PetscCheck(flg, PETSC_COMM_WORLD, PETSC_ERR_USER, "Must indicate whether or not to <n> with the -hnum option (1=true, 0=false)");
  PetscCall(PetscOptionsGetString(NULL, NULL, "-hSzblk", szblk, sizeof(szblk), NULL));
  PetscCheck(flg, PETSC_COMM_WORLD, PETSC_ERR_USER, "Indicate sz blocks with -hSzblk option which inputs pairs of numbers separated by , (e.g. 1,2,3,4)");
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-hhasW",&hasW,NULL));
  PetscCheck(flg, PETSC_COMM_WORLD, PETSC_ERR_USER, "Indicate whether or not the graph has weight -hhasW option (1=true, 0=false)");
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-hhasV",&hasV,NULL));
  PetscCheck(flg, PETSC_COMM_WORLD, PETSC_ERR_USER, "Indicate whether or not the graph has nearest neighbor repulsion -hhasV option (1=true, 0=false)");

  /* 
   * Read the Sz Blocks into an array
   *
   */
  size_t* SzBlock = malloc(MAX_Sz_BLOCKS * sizeof(size_t));
  size_t* NumBlock = malloc(MAX_Sz_BLOCKS * sizeof(size_t));
  char *p = szblk;
  int nblk = 0;
  setSzBlockList(p, &nblk, SzBlock) ;
  for(size_t k=0;k<nblk;++k) NumBlock[k] = SzBlock[k];
  //printf(" doSz = %ld nblk=%d HASV=%ld ---- \n",DoSz, nblk,hasV);

  igraph_t graph;
  igraph_set_attribute_table(&igraph_cattribute_table);
  igraph_empty(&graph, 0, IGRAPH_DIRECTED);
  igraph_integer_t num_vertices;

  if (readGraphMLFile(graphmlFile, &graph)) {
    // Successfully read the graph, now you can work with 'graph'.
    num_vertices = igraph_vcount(&graph);
  }

  // Print all the edge weights
  int wrows = num_vertices;
  int wcols = wrows;
  double** wmatrix = declare_matrix(wrows, wcols);
  double** vmatrix = declare_matrix(wrows, wcols);

  getWeightMatrix(&graph, wmatrix, (size_t)hasW);
  getRepulsionMatrix(&graph, vmatrix, (size_t)hasV);
  //print_matrix_d(vmatrix, wrows, wcols);
  //for (size_t i = 0; i < num_vertices; ++i) {
  //  igraph_vector_int_t orbital_id_allowed;
  //  igraph_vector_int_init(&orbital_id_allowed, 0);
  //  getConnectedVertices(&graph, (igraph_integer_t)i, &orbital_id_allowed);
  //  for (size_t j = 0; j < igraph_vector_int_size(&orbital_id_allowed); ++j) {
  //      size_t orbital_id = VECTOR(orbital_id_allowed)[j];
  //      printf(" %ld %ld \n",i,orbital_id);
  //  }
  //}

  // Assume configAlpha and configBeta are sorted lists of all possible alpha and beta configurations
  size_t norb = num_vertices;
  size_t nbeta = nelec - nalpha;
  int natomax = 100;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Define Hamiltonian
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  double U =  Ut_inp;
  double t = t_inp;

  size_t sizeAlpha = binomialCoeff(norb, nalpha);
  size_t sizeBeta = binomialCoeff(norb, nbeta);

  size_t* configAlpha = malloc(sizeAlpha * sizeof(size_t));
  size_t* configBeta = malloc(sizeBeta * sizeof(size_t));
  size_t sizeTotal  = sizeAlpha * sizeBeta;

  int max_nbrs = getMaxNeighbors(&graph, norb);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"==========================================="));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nHubHam: Hubbard Eigenproblem, n=%" PetscInt_FMT "\n",sizeTotal));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"===========================================\n\n"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Info] Norbs      \t\t %" PetscInt_FMT "\n",(size_t)norb));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Info] Nalpha     \t\t %" PetscInt_FMT "\n",(size_t)nalpha));
  if(DBGPrinting) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Info] Nbeta      \t\t %" PetscInt_FMT "\n",(size_t)nbeta));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Info] N(betas)  \t\t %" PetscInt_FMT "\n",sizeBeta));
  }
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Info] Max Nbrs   \t\t %" PetscInt_FMT "\n",(size_t)max_nbrs));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Info] N(alphas) \t\t %" PetscInt_FMT "\n",sizeAlpha));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Info] t          \t\t %10.5f |t| \n",(double)t_inp));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Info] U          \t\t %10.5f |t| \n",(double)Ut_inp));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Wait] Generating CSFs... \t\t \n"));

  generateConfigurations(norb, nalpha, configAlpha, &sizeAlpha);
  generateConfigurations(norb, nbeta, configBeta, &sizeBeta);

  // Sort the lists for binary search
  qsort(configAlpha, sizeAlpha, sizeof(size_t), compare);
  qsort(configBeta, sizeBeta, sizeof(size_t), compare);
  
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," [Done] Generating CSFs !  \t\t \n"));

  // Declare a matrix of size 3 x 4
  int rows = sizeAlpha * sizeBeta;
  int cols = rows;
  //double** matrix = declare_matrix(rows, cols);


  PetscInt       n=sizeAlpha*sizeBeta,i,Istart,Iend,nev,maxit,its,nconv;
  PetscInt		   ncv, mpd;

  // General Matrix
  //PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  //PetscCall(MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,n,n,3*num_vertices,NULL,3*num_vertices,NULL,&A));
  //PetscCall(MatMPIAIJSetPreallocation(A,3*num_vertices,NULL,3*num_vertices,NULL));
  //PetscCall(MatSetFromOptions(A));
  //PetscCall(MatSetUp(A));
  // Symmetric Matrix
  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(MatCreateSBAIJ(PETSC_COMM_WORLD,1,PETSC_DECIDE,PETSC_DECIDE,n,n,4 + num_vertices,NULL,4 + num_vertices,NULL,&A));
  //PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n));
  PetscCall(MatSetType ( A, MATSBAIJ));
  PetscCall(MatMPIBAIJSetPreallocation(A,1,4 + num_vertices,NULL,4 + num_vertices,NULL));

  // To print the matrix
  //PetscViewer viewer; // declare a viewer object
  //PetscViewerASCIIGetStdout(PETSC_COMM_WORLD, &viewer); // get the standard output
  //
  /*
   * Matrix for the S2 operator
    */
  //PetscCall(MatCreate(PETSC_COMM_WORLD,&S2));
  //PetscCall(MatSetSizes(S2,PETSC_DECIDE,PETSC_DECIDE,n,n));
  //PetscCall(MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,n,n,4*num_vertices,NULL,4*num_vertices,NULL,&S2));
  //PetscCall(MatMPIAIJSetPreallocation(S2,4*num_vertices,NULL,4*num_vertices,NULL));
  // Symmetric Matrix
  //PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n));
  /*
   * Matrix for the S2 operator
    */
  if(DoS2){
    // Symmetric Matrix
    PetscCall(MatCreate(PETSC_COMM_WORLD,&S2));
    PetscCall(MatCreateSBAIJ(PETSC_COMM_WORLD,1,PETSC_DECIDE,PETSC_DECIDE,n,n,norb*num_vertices,NULL,norb*num_vertices,NULL,&S2));
    PetscCall(MatSetType (S2, MATSBAIJ));
    PetscCall(MatMPIBAIJSetPreallocation(S2,1,norb*num_vertices,NULL,norb*num_vertices,NULL));
  }

  //PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n Number of vertices: %ld",(long)num_vertices));
  //PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n Number of configurations alpha = %ld, beta = %ld\n",sizeAlpha,sizeBeta));
  PetscCall(PetscTime(&tt1));

  /*
   * Initialize Hamiltonian
    */
  PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));
  for (i=Istart;i<Iend;i++) {

    igraph_vector_t MElist;
    igraph_vector_init(&MElist, 0);
    igraph_vector_t Jdetlist;
    igraph_vector_init(&Jdetlist, 0);

    int diag = getHubbardDiag(i, configAlpha, sizeAlpha, configBeta, sizeBeta);
    
    if(hasV) {

      double MElistV = 0.0;

      double Udiag = diag*U;
      getHubbardDiagVij(i, configAlpha, sizeAlpha, configBeta, sizeBeta, &graph, &MElistV, vmatrix);
      //printf(" %ld U=%10.5f MEV = %10.5f tot=%10.5f\n",i,Udiag,MElistV,Udiag+MElistV);
      PetscCall(MatSetValue(A,i,i,(PetscReal)(Udiag + MElistV),INSERT_VALUES));
    }
    else{
      double Udiag = diag*U;
      //printf(" %ld U=%10.5f  tot=%10.5f\n",i,Udiag,Udiag);
      PetscCall(MatSetValue(A,i,i,(PetscReal)diag*U,INSERT_VALUES));
    }
    //matrix[i][i] = (PetscReal)diag*U;
    getAllHubbardMEs(i, &MElist, &Jdetlist, configAlpha, sizeAlpha, configBeta, sizeBeta, &graph, wmatrix);
    for (int j = 0; j < igraph_vector_size(&Jdetlist); ++j) {
      PetscInt Jid = VECTOR(Jdetlist)[j];
      //matrix[i][Jid] = t*(PetscReal)VECTOR(MElist)[j];
      //matrix[Jid][i] = t*(PetscReal)VECTOR(MElist)[j];
      if( i > Jid ) PetscCall(MatSetValue(A,Jid,i,t*(PetscReal)VECTOR(MElist)[j],INSERT_VALUES));
      else          PetscCall(MatSetValue(A,i,Jid,t*(PetscReal)VECTOR(MElist)[j],INSERT_VALUES));
      //PetscCall(MatSetValue(A,Jid,i,t*(PetscReal)VECTOR(MElist)[j],INSERT_VALUES));
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
  //PetscCall(MatView(A, viewer));

  if(DoS2) {
    /*
     * Initialize S2 operator
      */
    PetscCall(MatGetOwnershipRange(S2,&Istart,&Iend));
    for (i=Istart;i<Iend;i++) {

      igraph_vector_t MElist;
      igraph_vector_init(&MElist, 0);
      igraph_vector_t Jdetlist;
      igraph_vector_init(&Jdetlist, 0);

      getS2Operator(i, &MElist, &Jdetlist, configAlpha, sizeAlpha, configBeta, sizeBeta, &graph, num_vertices, natomax);
      for (int j = 0; j < igraph_vector_size(&Jdetlist); ++j) {
        PetscInt Jid = VECTOR(Jdetlist)[j];
        //matrix[i][Jid] = VECTOR(MElist)[j];
        //printf(" %d %10.5f \n",Jid, VECTOR(MElist)[j]);
        if( i > Jid) PetscCall(MatSetValue(S2,Jid,i,t*(PetscReal)VECTOR(MElist)[j],INSERT_VALUES));
        else         PetscCall(MatSetValue(S2,i,Jid,t*(PetscReal)VECTOR(MElist)[j],INSERT_VALUES));
        //PetscCall(MatSetValue(S2,Jid,i,t*(PetscReal)VECTOR(MElist)[j],INSERT_VALUES));
      }

      igraph_vector_destroy(&MElist);
      igraph_vector_destroy(&Jdetlist);
    }
    PetscCall(PetscTime(&tt2));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," Time used to build the S2 operator: %f\n",tt2-tt1));


    PetscCall(PetscTime(&tt1));
    PetscCall(MatAssemblyBegin(S2,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(S2,MAT_FINAL_ASSEMBLY));
    PetscCall(PetscTime(&tt2));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," Time used to assemble the matrix: %f\n",tt2-tt1));
  }

  PetscCall(MatCreateVecs(A,NULL,&xr));
  PetscCall(MatCreateVecs(A,NULL,&xi));
  PetscCall(MatCreateVecs(A,NULL,&vs2));

  // Save file
  //save_matrix(matrix, rows, cols, "/tmp/benzene_c.csv");

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
  tol = 1.e-16;
  maxit = 10000000;
  PetscCall(EPSSetTolerances(eps,tol,maxit));
  //ncv  = 9;
  //mpd  = 10;
  //nev  = 4;
  //PetscCall(EPSSetDimensions(eps,nev,PETSC_DECIDE,PETSC_DECIDE));
  PetscCall(EPSSetFromOptions(eps));

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
    if(DoSz){
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,
         "           k          ||Ax-kx||/||kx||         S2                 Sz   \n"
         "   ----------------- -----------------------------------------------------\n"));
    }
    else if(DoNum){
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,
         "           k          ||Ax-kx||/||kx||         S2                 <N>  \n"
         "   ----------------- -----------------------------------------------------\n"));
    }
    else{
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,
         "           k          ||Ax-kx||/||kx||         S2         \n"
         "   ----------------- -------------------------------------\n"));
    }

    for (i=0;i<nev;i++) {
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

      /*
       * Get Spin S2 value
       */
      if(DoS2) {
        /*
         * Get Spin S2 value
         */
        PetscCall(MatMult(S2, xr, vs2));
        PetscCall(VecDot(xr, vs2, &dot));

        spin = solveQuad(1.0, 1.0, -1.0*fabs(dot));
      }
      else spin = 100;

      //double** szvalAll = declare_matrix(num_vertices, num_vertices);
      //double* szvalAll = (double *)malloc(num_vertices * num_vertices * sizeof(double));
      double szvalAll[num_vertices * num_vertices];
      if(DoSz) {
        /*
         * Get Sz values
         */
        //double** szvalpsi = declare_matrix(num_vertices, num_vertices);
        //double* szvalpsi = (double *)malloc(num_vertices * num_vertices * sizeof(double));
        double szvalpsi[num_vertices * num_vertices];
        double szval[2]; 
        for(size_t i1=0;i1<num_vertices;++i1) {
        for(size_t i2=0;i2<num_vertices;++i2) {
          szvalpsi[i1*num_vertices + i2] = 0.0;
          szvalAll[i1*num_vertices + i2] = 0.0;
        }
        }
        for(size_t k=0;k<2;++k) {
          szval[k] = 0.0; 
        }
        PetscCall(VecGetOwnershipRange(xr,&Istart,&Iend));
        for (size_t j=Istart;j<Iend;j++) {
          PetscInt ix[1];
          PetscScalar y[1];
          ix[0] = (size_t)j;
          PetscCall(VecGetValues(xr, 1, ix, y));
          size_t posi = j;
          size_t alphaID = findAlphaID(posi, sizeAlpha, sizeBeta);
          size_t betaID  = findBetaID(posi, sizeAlpha, sizeBeta);
          size_t detIa[1];
          size_t detIb[1];
          detIa[0] = configAlpha[alphaID];
          detIb[0] = configBeta[betaID];
          for(size_t i1=0;i1<num_vertices;++i1) {
            for(size_t i2=0;i2<num_vertices;++i2) {
              SzBlock[0] = i1;
              SzBlock[1] = i2;
              getSzOperator(detIa[0], detIb[0], szval, configAlpha, sizeAlpha, configBeta, sizeBeta, 2, SzBlock) ;
              double szprodpsi=1.0;
              szprodpsi = szval[0]*szval[1];
              szvalpsi[i1*num_vertices + i2] += szprodpsi*y[0]*y[0];;
            }
          }

          //printf(" --> %d \n",isDiag);
          if(DBGPrinting) {
            PetscCall(PetscPrintf(PETSC_COMM_WORLD,"%12f \t",y[0]));
            printf("\t alpha: ");
						for(size_t l=0; l < nalpha; ++l) {
              printf(" %d ",(1<<(l)) & detIa[0]? 1 : 0);
            }
            printf("\t CSF: ");
						for(size_t l=0; l < nbeta; ++l) {
              printf(" %d ",(1<<(l)) & detIb[0]? 1 : 0);
            }
            PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n"));
          }
        }
        MPI_Reduce(&szvalpsi, &szvalAll, num_vertices*num_vertices, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
      }

      double numvalAll[nblk]; 
      double weightNAll;
      if(DoNum) {
        /*
         * Get electron Number values
         *
         * The code can be used to calculate the
         * total weight of the neutral configurations 
         * or the weight of the <N>==1 on any combination
         * of sites as follows:
         *
         * -> with -hSblk 0,1,2,3,4,5,6,7 we get the total
         * weight of the neutral determinants.
         *
         * -> with -hSblk 0 we get the probability of 
         *  finding 1 electron on the site labelled 0.
         *
         */
        double numvalpsi[nblk]; 
        double weightNpsi;
        double numval[nblk]; 
        double numvala[nblk]; 
        int    isN=0;
        for(size_t k=0;k<nblk;++k) {
          numvalAll[k] = 0.0;
          numvalpsi[k] = 0.0; 
          numval[k] = 0.0; 
          numvala[k] = 0.0; 
        }
        weightNpsi=0.0;
        PetscCall(VecGetOwnershipRange(xr,&Istart,&Iend));
        for (size_t j=Istart;j<Iend;j++) {
          PetscInt ix[1];
          PetscScalar y[1];
          ix[0] = (size_t)j;
          PetscCall(VecGetValues(xr, 1, ix, y));
          size_t posi = j;
          size_t alphaID = findAlphaID(posi, sizeAlpha, sizeBeta);
          size_t betaID  = findBetaID(posi, sizeAlpha, sizeBeta);
          size_t detIa[1];
          size_t detIb[1];
          detIa[0] = configAlpha[alphaID];
          detIb[0] = configBeta[betaID];
          isN = 0;

          getNumOperator(detIa[0], detIb[0], numval, numvala, configAlpha, sizeAlpha, configBeta, sizeBeta, nblk, NumBlock) ;
          for(size_t k=0;k<nblk;++k) {
             numvalpsi[k] += numvala[k]*y[0]*y[0];
             if(abs(numval[k] - 1.0) < 10E-12) isN+=1;
          }
          if(isN==nblk) weightNpsi += y[0]*y[0];
          //printf(" --> %d \n",isDiag);
          if(DBGPrinting) {
            PetscCall(PetscPrintf(PETSC_COMM_WORLD,"%12f \t",y[0]));
            printf("\t alpha: ");
						for(size_t l=0; l < nalpha; ++l) {
              printf(" %d ",(1<<(l)) & detIa[0]? 1 : 0);
            }
            printf("\t CSF: ");
						for(size_t l=0; l < nbeta; ++l) {
              printf(" %d ",(1<<(l)) & detIb[0]? 1 : 0);
            }
            PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n"));
          }
        }
        MPI_Reduce(&numvalpsi, &numvalAll, nblk, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
        MPI_Reduce(&weightNpsi, &weightNAll, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
      }


      if (im!=0.0) PetscCall(PetscPrintf(PETSC_COMM_WORLD," %9f%+9fi %12g\n",(double)re,(double)im,(double)error));
      else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"   %12f       %12g       %12f",(double)re,(double)error,(double)fabs(spin)));
      //PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n"));
      //if(DoSz) {
      //  int a=NumBlock[0];
      //  int b=NumBlock[1];
      //  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"       %8.5f ",(double)szvalAll[a][b]));
      //}
      if(DoSz) {
        for(size_t l=0;l<nblk;l=l+2) {
          size_t i1 = NumBlock[l];
          size_t i2 = NumBlock[l+1];
          PetscCall(PetscPrintf(PETSC_COMM_WORLD,"       %8.5f ",(double)szvalAll[i1*num_vertices + i2]));
        }
      }
      if(DoNum) {
        //for(size_t k=0;k<nblk;++k) {
        //  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"       %8.5f ",(double)numvalAll[k]));
        //}
        PetscCall(PetscPrintf(PETSC_COMM_WORLD,"     W = %8.5f ",weightNAll));
      }
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n"));

    }
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n"));
  }

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
  PetscCall(VecDestroy(&vs2));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"===========================================\n"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"HubHam: Success !!!                      \n"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"===========================================\n"));
  PetscCall(SlepcFinalize());
  return 0;
}
