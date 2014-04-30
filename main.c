#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <pthread.h>

#define SIZE 128

pthread_mutex_t lock;

inline unsigned long long rdtsc(void)
{
        unsigned long long result=0;
        unsigned a, d;

        do {
                __asm__ __volatile__("rdtsc" : "=a" (a), "=d" (d));
                result = ((unsigned long long)a) | (((unsigned long long)d) << 32);
        } while (__builtin_expect ((int) result == -1, 0));
        return result;
}

void * tm_run (void * arg);
void * mutex_run (void * arg);
void * basic_run (void * arg);

int threads = 2;


int main(int argc, char** argv)
{
	
    int myrank, numprocs, subgroup_rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm mu_comm, tm_comm, ba_comm, copy_comm;
	MPI_Group original, mu_group, tm_group, ba_group;
	
	
    int *ranks1= malloc((numprocs/2)*(sizeof(int)));
    int *ranks2= malloc((numprocs/2)*(sizeof(int)));
    int *ranks3= malloc((numprocs/2)*(sizeof(int)));
	
	MPI_Comm_dup(MPI_COMM_WORLD, &copy_comm);
	MPI_Comm_group(copy_comm, &original);
	MPI_Comm_group(copy_comm, &tm_group);
	MPI_Comm_group(copy_comm, &mu_group);
	MPI_Comm_group(copy_comm, &ba_group);
	
    int a[SIZE], i;
	float tm_runtime = 0.0;
	float mutex_runtime = 0.0;
	float basic_runtime = 0.0;
    
    for (i = 0; i < numprocs; i++)
    {
		if (i < numprocs/2)
		{
			ranks2[i] = i;
		}
		else
		{
			ranks3[i - numprocs/2] = i;
		}
	}
    
	if (myrank < numprocs/2){
		MPI_Group_incl(original, numprocs/2, ranks2, &mu_group);
	}
	else{
		MPI_Group_incl(original, numprocs/2, ranks3, &ba_group);
    }
    
	MPI_Comm_create(copy_comm, tm_group, &tm_comm);
	MPI_Comm_create(copy_comm, mu_group, &mu_comm);
	MPI_Comm_create(copy_comm, ba_group, &ba_comm);
    	

    pthread_t t_id[threads];
    pthread_attr_t attr;
    pthread_attr_init(&attr );
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    int t;
    float res = 0.0;
    for (t = 0; t < threads; t++){
        //printf("New Thread!\n");
        /*if (myrank < numprocs/3)
          if (pthread_create(&t_id[t],&attr,&tm_run,(void*)&a))
              printf("\nError during thread creation!\n");
        else*/ if (myrank < numprocs/2){
          if (pthread_create(&t_id[t],&attr,&mutex_run,(void*)&a))
              printf("\nError during thread creation!\n");
		  }
        else
          if (pthread_create(&t_id[t],&attr,&basic_run,(void*)&a))
              printf("\nError during thread creation!\n");
    }
    pthread_attr_destroy(&attr);//printf("%i: attr destroyed\n",myrank);
    //Threads join with original again
    for (t = 0; t < threads; t++){
		//printf("RANK: %d THREAD: %d\n", myrank, t);
        //if(!pthread_join(t_id[i],NULL)) printf("%i: Thread %i exited with error.\n",myrank,i);
        
        void *thread_ret;
        int rv = pthread_join(t_id[t], &thread_ret);
        if (rv!=0){
            printf("%i: Thread %i exited with error %i .\n",myrank,t,rv);
        }
        else{
			res += *((float *)thread_ret);
			free(thread_ret);
            //printf("%i: %i:  %f: thread joined successfully\n",myrank,t,res);
            if (myrank < numprocs/2){
				printf("mutex reduce!\n");
			}
			else{
				printf("basic reduce!\n");
			}
        }
    }
    
    if (myrank < numprocs/2){
		MPI_Allreduce(&res, &mutex_runtime, 1, MPI_FLOAT, MPI_SUM, mu_comm);
	}
	else{
		MPI_Allreduce(&res, &basic_runtime, 1, MPI_FLOAT, MPI_SUM, ba_comm);
	}
   int b[SIZE];
 MPI_Request rrequest;
  MPI_Request srequest;
MPI_Status status;
   int flag=0;
  MPI_Irecv(b, SIZE, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &rrequest);  
  MPI_Isend(a, SIZE, MPI_INT,(myrank+(numprocs/2))%numprocs, 0, MPI_COMM_WORLD, &srequest);
 while(flag==0)
  {
     MPI_Test(&rrequest,&flag,&status);
  }
flag=0;
 while(flag==0)
  {
     MPI_Test(&srequest,&flag,&status);
   }
for(t=0;t<SIZE;t++)
  if(a[t]!=b[t])
     {printf("rank %i and %i are different\n",myrank,(myrank+(numprocs/2))%numprocs);break;}
    MPI_Finalize();
    if (myrank == 0){
		printf("Mutex runtime: %f\n", mutex_runtime);
	}
	else if (myrank == numprocs/2){
		printf("Basic runtime: %f\n", basic_runtime);
	}
	
    if (myrank==0)printf("Transactional memory test completed\n");
    return 0;
}
/*
void * tm_run (void * arg){

    int v, w, z;
    int b[SIZE];
    int * a = arg;
    float * runtime = (float *)malloc(sizeof(float));

    unsigned long long t0, t1;

    for (v=0; v<SIZE; v++)
    {
       // a[v] = v;
        b[v] = -v;
        printf("tm %i: b = %i\n",v,b[v]);
    }
    printf("For loop done\n");
    // Transactional Memory Run
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
    *runtime = (float)(t1-t0)/1600000000;
    printf("tm_atomic    ran in %f seconds \n", (float)(t1-t0)/1600000000);
    pthread_exit(runtime);
    return NULL;
}
* */
void * mutex_run (void * arg){

    int v, w, z;
    float res = 0.0;
    int b[SIZE];
    int * a = arg;
    float * runtime = (float *)malloc(sizeof(float));
    unsigned long long t0, t1;

        for (v=0; v<SIZE; v++)
        {
            a[v] = v;
            b[v] = -v;
            //printf("mu %i: b = %i\n",v,b[v]);
        }
        //Pthread Mutex Lock
        t0 = rdtsc();
        for (v=0; v<SIZE; v++)
          for (w=0; w<SIZE; w++)
            for (z=0; z<SIZE; z++)
            {
              pthread_mutex_lock(&lock);
              a[z] = a[w] + b[v];
              pthread_mutex_unlock(&lock);
            }
        t1 = rdtsc();
        *runtime = (float)(t1-t0)/1600000000;
        printf("mutex_lock   ran in %f seconds \n", *runtime);
        pthread_exit(runtime);
        return NULL;
}
void * basic_run (void * arg){

    int v, w, z;
    float test = 1.0;
    int b[SIZE];
    int * a = arg;
    float * runtime = (float *)malloc(sizeof(float));
    unsigned long long t0, t1;
        for (v=0; v<SIZE; v++)
        {
            a[v] = v;
            b[v] = -v;
            //printf("na %i: b = %i\n",v,b[v]);
        }
        //No lock or TM run
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
        *runtime = (float)(t1-t0)/1600000000;
        printf("non-atomic   ran in %f seconds \n", *runtime);
		pthread_exit(runtime);
		return NULL;
}	
