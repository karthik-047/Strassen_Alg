#include <ctime>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
int size; // Size to get from user

template <typename T> T **Allocate2DArray(int nRows, int nCols) {
  //(step 1) allocate memory for array of elements of column
  T **ppi = new T *[nRows];

  //(step 2) allocate memory for array of elements of each row
  T *curPtr = new T[nRows * nCols];

  // Now point the pointers in the right place
  for (int i = 0; i < nRows; ++i) {
    *(ppi + i) = curPtr;
    curPtr += nCols;
  }
  return ppi;
}

// Print matrix of size n
void print_matrix(int **mat, int n) {
  for (int i = 0; i < n; i++) {
    printf("\n");
    for (int j = 0; j < n; j++) {
      printf("\t %d", mat[i][j]);
    }
  }
}

// Deallocate memory
template <typename T> void Free2DArray(T **Array) {
  delete[] * Array;
  delete[] Array;
}

// Final result matrix using Strassen
int **C_top;

// Terminal matrix size below which serial matrix mult takes place
int THRESHOLD;

// Conventional matrix multiplication
void seqMatMult(int m, int n, int p, int **A, int **B, int **C) {
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++) {
      printf("\n");
      for (int k = 0; k < p; k++)
        C[i][j] += A[i][k] * B[k][j];
    }
}

// To perform multiplication when you hit threshold
void matmultleaf(int mf, int ml, int nf, int nl, int pf, int pl, int **A,
                 int **B, int **C)
/*
  subroutine that uses the simple triple loop to multiply
  a submatrix from A with a submatrix from B and store the
  result in a submatrix of C.
*/
// mf, ml; /* first and last+1 i index */
// nf, nl; /* first and last+1 j index */
// pf, pl; /* first and last+1 k index */
{
  // printf("\n Leaf check");
  // printf("\n %d %d %d %d %d %d", mf, ml, nf, nl, pf, pl);
  for (int i = mf; i < ml; i++)
    for (int j = nf; j < nl; j++) {
      C[i][j] = 0.0;
      for (int k = pf; k < pl; k++)
        C[i][j] += A[i][k] * B[k][j];
    }
}

// Copy parent matrix into submatrices
void copyQtrMatrix(int **X, int m, int **Y, int mf, int nf) {
  printf("\n");
  for (int i = 0; i < m; i++) {
    X[i] = &Y[mf + i][nf];
  }
  // print_matrix(X, m);
}

// Add matrix blocks
void AddMatBlocks(int **T, int m, int n, int **X, int **Y) {
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      T[i][j] = X[i][j] + Y[i][j];
}

// Subtract Matrix blocks
void SubMatBlocks(int **T, int m, int n, int **X, int **Y) {
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      T[i][j] = X[i][j] - Y[i][j];
}

// Strassen Algorithm
void strassenMMult(int ml, int nl, int pl, int **A, int **B, int **X) {
  // printf("\nStrassencalls: %d %d %d", ml, nl, pl);

  // Check either of the matrix dims is <= threshold to perform og method (ml, nl, pl have same values)
  if (ml <= THRESHOLD) {
    // printf("\nMatmulleaf %d %d %d", ml, nl, pl);
    matmultleaf(0, ml, 0, nl, 0, pl, A, B, X);

  } else {
    int m2 = ml / 2;
    int n2 = nl / 2;
    int p2 = pl / 2;

    int **M1 = Allocate2DArray<int>(m2, n2);
    int **M2 = Allocate2DArray<int>(m2, n2);
    int **M3 = Allocate2DArray<int>(m2, n2);
    int **M4 = Allocate2DArray<int>(m2, n2);
    int **M5 = Allocate2DArray<int>(m2, n2);
    int **M6 = Allocate2DArray<int>(m2, n2);
    int **M7 = Allocate2DArray<int>(m2, n2);

    int **wAM1 = Allocate2DArray<int>(m2, p2);
    int **wBM1 = Allocate2DArray<int>(p2, n2);
    int **wAM2 = Allocate2DArray<int>(m2, p2);
    int **wBM3 = Allocate2DArray<int>(p2, n2);
    int **wBM4 = Allocate2DArray<int>(p2, n2);
    int **wAM5 = Allocate2DArray<int>(m2, p2);
    int **wAM6 = Allocate2DArray<int>(m2, p2);
    int **wBM6 = Allocate2DArray<int>(p2, n2);
    int **wAM7 = Allocate2DArray<int>(m2, p2);
    int **wBM7 = Allocate2DArray<int>(p2, n2);

    int **A11 = new int *[m2];
    int **A12 = new int *[m2];
    int **A21 = new int *[m2];
    int **A22 = new int *[m2];

    int **B11 = new int *[p2];
    int **B12 = new int *[p2];
    int **B21 = new int *[p2];
    int **B22 = new int *[p2];

    int **C11 = new int *[m2];
    int **C12 = new int *[m2];
    int **C21 = new int *[m2];
    int **C22 = new int *[m2];

    copyQtrMatrix(A11, m2, A, 0, 0);
    copyQtrMatrix(A12, m2, A, 0, p2);
    copyQtrMatrix(A21, m2, A, m2, 0);
    copyQtrMatrix(A22, m2, A, m2, p2);

    copyQtrMatrix(B11, p2, B, 0, 0);
    copyQtrMatrix(B12, p2, B, 0, n2);
    copyQtrMatrix(B21, p2, B, p2, 0);
    copyQtrMatrix(B22, p2, B, p2, n2);

    copyQtrMatrix(C11, m2, X, 0, 0);
    copyQtrMatrix(C12, m2, X, 0, n2);
    copyQtrMatrix(C21, m2, X, m2, 0);
    copyQtrMatrix(C22, m2, X, m2, n2);

    // M1 = (A11 + A22)*(B11 + B22)
    printf("\nM1");
    AddMatBlocks(wAM1, m2, p2, A11, A22);
    AddMatBlocks(wBM1, p2, n2, B11, B22);

    strassenMMult(m2, n2, p2, wAM1, wBM1, M1);

    // M2 = (A21 + A22)*B11
    AddMatBlocks(wAM2, m2, p2, A21, A22);
    printf("\nM2");
    strassenMMult(m2, n2, p2, wAM2, B11, M2);

    // M3 = A11*(B12 - B22)
    SubMatBlocks(wBM3, p2, n2, B12, B22);
    printf("\nM3");
    strassenMMult(m2, n2, p2, A11, wBM3, M3);

    // M4 = A22*(B21 - B11)
    SubMatBlocks(wBM4, p2, n2, B21, B11);
    printf("\nM4");
    strassenMMult(m2, n2, p2, A22, wBM4, M4);

    // M5 = (A11 + A12)*B22
    AddMatBlocks(wAM5, m2, p2, A11, A12);
    printf("\nM5");
    strassenMMult(m2, n2, p2, wAM5, B22, M5);

    // M6 = (A21 - A11)*(B11 + B12)
    SubMatBlocks(wAM6, m2, p2, A21, A11);
    AddMatBlocks(wBM6, p2, n2, B11, B12);
    printf("\nM6");
    strassenMMult(m2, n2, p2, wAM6, wBM6, M6);

    // M7 = (A12 - A22)*(B21 + B22)
    SubMatBlocks(wAM7, m2, p2, A12, A22);
    AddMatBlocks(wBM7, p2, n2, B21, B22);
    printf("\nM7");
    strassenMMult(m2, n2, p2, wAM7, wBM7, M7);

    for (int i = 0; i < m2; i++)
      for (int j = 0; j < n2; j++) {
        C11[i][j] = M1[i][j] + M4[i][j] - M5[i][j] + M7[i][j];
        C12[i][j] = M3[i][j] + M5[i][j];
        C21[i][j] = M2[i][j] + M4[i][j];
        C22[i][j] = M1[i][j] - M2[i][j] + M3[i][j] + M6[i][j];
      }

    Free2DArray<int>(M1);
    Free2DArray<int>(M2);
    Free2DArray<int>(M3);
    Free2DArray<int>(M4);
    Free2DArray<int>(M5);
    Free2DArray<int>(M6);
    Free2DArray<int>(M7);

    Free2DArray<int>(wAM1);
    Free2DArray<int>(wBM1);
    Free2DArray<int>(wAM2);
    Free2DArray<int>(wBM3);
    Free2DArray<int>(wBM4);
    Free2DArray<int>(wAM5);
    Free2DArray<int>(wAM6);
    Free2DArray<int>(wBM6);
    Free2DArray<int>(wAM7);
    Free2DArray<int>(wBM7);

    delete[] A11;
    delete[] A12;
    delete[] A21;
    delete[] A22;
    delete[] B11;
    delete[] B12;
    delete[] B21;
    delete[] B22;
    delete[] C11;
    delete[] C12;
    delete[] C21;
    delete[] C22;
  }
}
// Called from main and calls strassen alg
void matmultS(int m, int n, int p, int **A, int **B, int **C) {
  strassenMMult(m, n, p, A, B, C);
}

int main(int argc, char *argv[]) {

  size = atoi(argv[1]);
  THRESHOLD = atoi(argv[2]);
  double start, end;

  int **A = Allocate2DArray<int>(size, size);
  int **B = Allocate2DArray<int>(size, size);
  C_top = Allocate2DArray<int>(size, size);         // Strassen
  int **C_check = Allocate2DArray<int>(size, size); // Traditional method

  int i, j;

  for (int i = 0; i < size; ++i)
    for (int j = 0; j < size; ++j) {
      B[i][j] = rand() % 3 + 1;
      A[i][j] = rand() % 3 + 1;
      C_top[i][j] = 0;
      C_check[i][j] = 0;
    }

  matmultS(size, size, size, A, B, C_top);
  seqMatMult(size, size, size, A, B, C_check);
  printf("\nMatrix A\n");
  print_matrix(A, size);
  printf("\nMatrix B\n");
  print_matrix(B, size);
  printf("\nMatrix C\n");
  print_matrix(C_top, size);
  printf("\nMatrix C serial\n");
  print_matrix(C_check, size);

  bool correctness = true;

  for (int i = 0; i < size; ++i)
    for (int j = 0; j < size; ++j)
      if (C_check[i][j] != C_top[i][j])
        correctness = false;
  if (correctness) {
    printf("\nOpenMP - Correct: matrix size = %d", size);

  } else {
    printf("\nOpenMP - Incorrect: matrix size = %d", size);
  }
  Free2DArray<int>(A);
  Free2DArray<int>(B);
  Free2DArray<int>(C_top);
  Free2DArray<int>(C_check);

  return 0;
}
