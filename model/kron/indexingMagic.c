#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// A space optimized Dynamic Programming Solution
// A utility function to return minimum of two integers
int min(int a, int b)
{
    return (a<b)? a: b;
}
int binomial(int n, int k)
{
    int* C = (int*)calloc(k+1, sizeof(int));
    int i, j, res;
 
    C[0] = 1;
 
    for(i = 1; i <= n; i++)
    {
        for(j = min(i, k); j > 0; j--)
            C[j] = C[j] + C[j-1];
    }
 
    res = C[k];  // Store the result before freeing memory
 
    free(C);  // free dynamically allocated memory to avoid memory leak
 
    return res;
}

const int K = 4;
const int N = 100;

static int C[N+K+1][K+2];

void populateC()
{
    for(int k = 0; k < K+1; ++k)
    {
        for(int n = 1; n < N+K+1; ++n)
        {
            C[n][k] = binomial(n+k-1, k);
        }
    }
}

int binomialCached(int n, int k)
{
    return C[n][k];
}

void largest(int i, int n, int k, int* val_off)
{
   int x = binomial(n+k-1, k);
   //int x = binomialCached(n, k);
   if (x > i)
   {
        largest(i, n-1, k, val_off);
   }
   else
   {
        val_off[0] = n;
        val_off[1] = x;
    }
}

static int gs_val_off[2]; //We reuse this and rely on the fact that nothing is multithreaded!???

void id2combo_0(int *combo, int idx, int n, int k, int* val_off)
{
    //This is based on http://stackoverflow.com/questions/12146910/finding-the-lexicographic-index-of-a-permutation-of-a-given-array
    //idx is assumed to be 0 based indexing. Combo will get populated assuming 1 based indexing!
    
    largest(idx,n,k,val_off);
    combo[k-1] = val_off[0] + 1; //k-1 since this is C!
    int offset = val_off[1];
    //don't waste a recursion call to try to populate 0
    if (k > 1)
    { 
        id2combo_0(combo,idx-offset, n, k-1, val_off);    
    }
}

void id2combo(int *combo, int idx, int n, int k)
{
    //We assume combo that idx is 1 based indexing
    //so we subtract 1 and turn it over to function which deals with 0
    //based indexing    
    id2combo_0(combo,idx-1, n, k, gs_val_off); //Ooo the joys of 1 based indexing ...
}
int main()
{
    populateC();

    int combo[5] = {0};
    id2combo(combo, 2, 3, 5);
    for(int i = 0; i < 5; ++i)
        printf("%i ", combo[i]);
    printf("\n");

    //run timing
    clock_t start = clock(), diff;

    for(int i = 1; i < 1000000; ++i)
    {
        id2combo(combo, i, 3, 5);
    }

    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);

    for(int i = 0; i < 5; ++i)
        printf("%i ", combo[i]);
    printf("\n");


    return 0;
}