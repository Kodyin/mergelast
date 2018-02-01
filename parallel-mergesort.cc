/**
 *  \file parallel-mergesort.cc
 *
 *  \brief Implement your parallel mergesort in this file.
 */

#include <assert.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
using namespace std;

#include "sort.hh"
void
parallelSort (int N, keytype* A);
void
mergeSort (keytype* A, int l, int r, keytype* B, int lvl);
// void
// mergeSort (keytype* A, int l, int r, int lvl);
void
merge (keytype* A, int l1, int r1, int l2, int r2, keytype* B, int l3, int lvl);
int
binaryS (int tar, keytype* A, int l, int r);
void
mergeSequential (keytype* A, int l1, int r1, int l2, int r2, keytype* B, int l3);

void
parallelSort (int N, keytype* A)
{
    //omp_set_nested(1);
    //omp_set_num_threads (4);
    keytype* M = newKeys(N);
    int depth = 0;
    #pragma omp parallel
    {
         #pragma omp single nowait 
        {
          	depth = omp_get_num_threads();
            mergeSort (A, 0, N-1, M, depth);
        }
    }
    free(M);
  /* Lucky you, you get to start from scratch */
}

void
mergeSort (keytype* A, int l, int r, keytype* B, int lvl){
        int n = r - l + 1;
        if(n<=1) 
            return;
        else {
            int mid = (l+r)>>1;
            // if(lvl>1)
            // {
            // 	#pragma omp task 
            //     mergeSort(A, l, mid, B, lvl/2);
            //     mergeSort(A, mid+1, r, B, lvl/2);
            //     #pragma omp taskwait
            // }
            // else {
            //     sequentialSort(mid-l+1,A+l);
            // 	sequentialSort(r-mid,A+mid+1);
            // }
            if(lvl>1){
            	#pragma omp task firstprivate(A,B)
            	mergeSort(A, l, mid, B, lvl/2);
            	mergeSort(A, mid+1, r, B, lvl/2);
            	#pragma omp taskwait
            }
            else {
            	mergeSort(A, l, mid, B, lvl);
            	mergeSort(A, mid+1, r, B, lvl);
            }
            merge(A, l, mid, mid+1, r, B, l, lvl);
            memcpy(A + l, B + l, (r - l + 1) * sizeof(keytype)); 
      }
    
}
int
binaryS (int tar, keytype* A, int l, int r){
    int low = l, high = max(l, r+1);
    while(low<high){
        int mid = (low+high)>>1;
        if(tar<=A[mid]) high = mid;
        else low = mid + 1;
    }
    return high;
}

void
merge (keytype* A, int l1, int r1, int l2, int r2, keytype* B, int l3, int lvl){
    //if(lvl<4) printf("Bhread %d is merging %d ~ %d with %d ~ %d \n", omp_get_thread_num(), l1, r1, l2, r2);        
        int n1 = r1 - l1 + 1, n2 = r2 - l2 + 1;
        if(lvl>1){
        	if(n1<n2){
	            swap(l1,l2);
	            swap(r1,r2);
	            swap(n1,n2);
	        }
	        if(n1<=0) return;
	        int mid1 = (l1+r1)>>1;
	        int mid2 = binaryS(A[mid1], A, l2, r2);
	        int mid3 = l3 + (mid1 - l1) + (mid2 - l2);
	        B[mid3] =  A[mid1];
	    
	    	#pragma omp task shared(A, B)
	        merge(A, l1, mid1-1, l2, mid2-1, B, l3, lvl/2);

	        merge(A, mid1+1, r1, mid2, r2, B, mid3+1, lvl/2);
	    	#pragma omp taskwait    
        } else {
        	mergeSequential (A, l1, r1, l2, r2, B, l3);
        }                                   
}

void
mergeSequential (keytype* A, int l1, int r1, int l2, int r2, keytype* B, int l3){
    int a = l1, b = l2;
    for(int i = 0; i<r1-l1+1+r2-l2+1; i++){
        if(a>r1) {
        	memcpy(B + l3 + i, A + b, (r2 - b + 1) * sizeof(keytype));
			break;
		}
        else if(b>r2) {
        	memcpy(B + l3 + i, A + a, (r1 - a + 1) * sizeof(keytype));
			break;
        }
        else if(A[a] < A[b]) B[l3+i] = A[a++];
        else B[l3+i] = A[b++];
    }
}

/* eof */
