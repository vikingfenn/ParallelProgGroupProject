#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <pthread.h>

#define SIZE 128

pthread_mutex_t lock;

static __inline__ unsigned long long rdtsc(void){
    unsigned long long int result = 0;
    unsigned long int upper, lower, tmp;
    __asm__ volatile(
        "0: \n"
        "\tmftbu %0 \n"
        "\tmftb %1 \n"
        "\tmftbu %2 \n"
        "\tcmpw %2,%0 \n"
        "\tbne 0b \n"
        : "=r"(upper),"=r"(lower),"=r"(tmp)
    );
    result = upper;
    result = result << 32;
    result = result|lower;
    return (result);
}

void * tm_run(void * arg);
void * mutex_run(void * arg);
void * basic_run(void * arg);

int threads;

int main(int argc, char** argv) {

    int myrank, numprocs;

    if (argc < 2) {
        fprintf(stderr, "Useage: tm_test <threads>\n");
        return EXIT_FAILURE;
    }

    threads = atol(argv[1]);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int *ranks1 = malloc((numprocs / 3)*(sizeof (int)));
    int *ranks2 = malloc((numprocs / 3)*(sizeof (int)));
    int *ranks3 = malloc((numprocs / 3)*(sizeof (int)));

    int a[SIZE], i;
    float tm_runtime = 0.0;
    float mutex_runtime = 0.0;
    float basic_runtime = 0.0;

    pthread_t t_id[threads];
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    int t;
    int k;
    float resb = 0.0,rest=0.0,resm=0.0;

    // Create threads

    for (t = 0; t < threads; t++) {

        if (myrank < numprocs / 3){
            if ((pthread_create(&t_id[t], &attr, &tm_run, (void*) &a)) == 1){
                fprintf(stderr, "%d: Error during tm thread creation!\n", myrank);
                return EXIT_FAILURE;
            }
        } else if (myrank < 2 * numprocs / 3) {
            if ((pthread_create(&t_id[t], &attr, &mutex_run, (void*) &a)) == 1){
                fprintf(stderr, "%d: Error during mutex thread creation!\n", myrank);
                return EXIT_FAILURE;
            }
        } else {
            if ((pthread_create(&t_id[t], &attr, &basic_run, (void*) &a)) == 1){
                fprintf(stderr, "%d: Error during basic thread creation!\n", myrank);
                return EXIT_FAILURE;
            }
        }
    }

    pthread_attr_destroy(&attr);

    // Join threads

    for (t = 0; t < threads; t++) {

        void *thread_ret;
        int rv = pthread_join(t_id[t], &thread_ret);

        if (rv != 0) {
            fprintf(stderr, "%i: Thread %i exited with error %i!\n", myrank, t, rv);
        } else {
            if (myrank < numprocs/3) {
                rest += *((float *) thread_ret);
            } else if (myrank < 2*numprocs/3) {
                resm += *((float *) thread_ret);
            } else {
                resb += *((float *) thread_ret);
            }

            free(thread_ret);
        }
    }

    printf("%f,%f,%f\n", rest, resm, resb);
    MPI_Allreduce(&rest, &tm_runtime, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&resm, &mutex_runtime, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&resb, &basic_runtime, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

    int b[SIZE];
    MPI_Request rrequest;
    MPI_Request srequest;
    MPI_Status status;

    int flag = 0;
    MPI_Irecv(b, SIZE, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &rrequest);
    MPI_Isend(a, SIZE, MPI_INT, (myrank + (numprocs / 2)) % numprocs, 0, MPI_COMM_WORLD, &srequest);

    while (flag == 0) {
        MPI_Test(&rrequest, &flag, &status);
    }

    flag = 0;

    while (flag == 0) {
        MPI_Test(&srequest, &flag, &status);
    }

    /*
    for (t = 0; t < SIZE; t++) {
        if (a[t] != b[t]) {
            printf("rank %i and %i are different\n", myrank, (myrank + (numprocs / 2)) % numprocs);
            break;
        }
    }
    */

    MPI_Finalize();

    if (myrank != 0) {
        return EXIT_SUCCESS;
    }

    printf("Basic runtime: %f\n", basic_runtime);
    printf("Mutex runtime: %f\n", mutex_runtime);
    printf("Tm runtime: %f\n", tm_runtime);
    printf("Transactional memory test completed\n");

    tm_print_all_stats();

    return EXIT_SUCCESS;
}

void * tm_run(void * arg) {

    int v, w, z;
    int b[SIZE];
    int * a = arg;
    float * runtime = (float *) malloc(sizeof (float));

    unsigned long long t0, t1;

    for (v = 0; v < SIZE; v++) {
        a[v] = v;
        b[v] = -v;
    }

    t0 = rdtsc();

    for (v = 0; v < SIZE; v++) {
        for (w = 0; w < SIZE; w++) {
            for (z = 0; z < SIZE; z++) {

                #pragma tm_atomic
                {
                    a[z] = a[w] + b[v];
                }

            }
        }
    }

    t1 = rdtsc();
    *runtime = (float) (t1 - t0) / 1600000000;

    pthread_exit(runtime);

    return (void*)runtime;
}

void * mutex_run(void * arg) {

    int v, w, z;
    float res = 0.0;
    int b[SIZE];
    int * a = arg;
    float * runtime = (float *) malloc(sizeof (float));
    unsigned long long t0, t1;

    for (v = 0; v < SIZE; v++) {
        a[v] = v;
        b[v] = -v;
    }

    t0 = rdtsc();

    for (v = 0; v < SIZE; v++) {
        for (w = 0; w < SIZE; w++) {
            for (z = 0; z < SIZE; z++) {

                pthread_mutex_lock(&lock);
                a[z] = a[w] + b[v];
                pthread_mutex_unlock(&lock);

            }
        }
    }

    t1 = rdtsc();
    *runtime = (float) (t1 - t0) / 1600000000;

    pthread_exit(runtime);

    return (void*)runtime;
}

void * basic_run(void * arg) {

    int v, w, z;
    float test = 1.0;
    int b[SIZE];
    int * a = arg;
    float * runtime = (float *) malloc (sizeof(float));
    unsigned long long t0, t1;
    for (v = 0; v < SIZE; v++) {
        a[v] = v;
        b[v] = -v;
    }

    t0 = rdtsc();

    for (v = 0; v < SIZE; v++) {
        for (w = 0; w < SIZE; w++) {
            for (z = 0; z < SIZE; z++) {

                a[z] = a[w] + b[v];

            }
        }
    }

    t1 = rdtsc();
    *runtime = (float) (t1 - t0) / 1600000000;

    pthread_exit(runtime);

    return (void*)runtime;
} 
