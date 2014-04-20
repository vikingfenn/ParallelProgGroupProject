#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <pthread.h>

#define SIZE 128

static __inline__ unsigned long long rdtsc(void){
  unsigned long long int result = 0;
  unsigned long int upper, lower, tmp;
  __asm__ volatile(
                "0:                  \n"
                "\tmftbu   %0           \n"
                "\tmftb    %1           \n"
                "\tmftbu   %2           \n"
                "\tcmpw    %2,%0        \n"
                "\tbne     0b         \n"
                : "=r"(upper),"=r"(lower),"=r"(tmp)
      );
  result = upper;
  result = result << 32;
  result = result|lower;
  return (result);
}

void * tm_run (void * arg);

int main(int argc, char** argv)
{
    int myrank, numprocs;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    int a[SIZE];

    int threads = 10;
    pthread_t t_id[threads];
    pthread_attr_t attr;
    pthread_attr_init(&attr );
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    int t;
    for (t = 0; t < threads; t++){
        //printf("New Thread!\n");
        if (pthread_create(&t_id[t],&attr,&tm_run,&a))\
            printf("\nError during thread creation!\n");
    }
    pthread_attr_destroy(&attr);//printf("%i: attr destroyed\n",myrank);
    //Threads join with original again
    for (t = 0; t < threads; t++){
        //if(!pthread_join(t_id[i],NULL)) printf("%i: Thread %i exited with error.\n",myrank,i);
        int rv = pthread_join(t_id[t],NULL);
        if (rv!=0){
            printf("%i: Thread %i exited with error %i .\n",myrank,t,rv);
        }
        else{
            //printf("%i: %i: thread joined successfully\n",myrank,i);
        }
    }
    if (myrank==0)printf("Transactional memory test completed\n");
    MPI_Finalize();
    return 0;
}

void * tm_run (void * arg){

    int v, w, z;
    int b[SIZE];
    int * a = arg;

    int r;
    for (r=0; r<3; r++)
    {
        unsigned long long t0, t1;

        for (v=0; v<SIZE; v++)
        {
            a[v] = v;
            b[v] = -v;
        }

        t0 = rdtsc();
        for (v=0; v<SIZE; v++)
          for (w=0; w<SIZE; w++)
            for (z=0; z<SIZE; z++)
            {
              #pragma tm_atomic
              {
                a[z] = a[w] + b[v];
              }
            }

        t1 = rdtsc();
        printf("tm_atomic    ran in %f seconds \n", (float)(t1-t0)/1600000000);

        for (v=0; v<SIZE; v++)
        {
            a[v] = v;
            b[v] = -v;
        }

        t0 = rdtsc();
        for (v=0; v<SIZE; v++)
          for (w=0; w<SIZE; w++)
            for (z=0; z<SIZE; z++)
            {
              {
                a[z] = a[w] + b[v];
              }
            }

        t1 = rdtsc();
        printf("non-atomic   ran in %f cycles \n", (float)(t1-t0)/1600000000);
    }
    return NULL;
}