#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <pthread.h>

#define CPU_FREQ 1600000000
#define SIZE 128

pthread_mutex_t lock;

//instruction counter used to time trials is returned
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

void *tm_run(void *arg);
void *mutex_run(void *arg);
void *basic_run(void *arg);

int threads;

int main(int argc, char **argv) {
    //Variable initialization
    
    int myrank, numprocs;
    //the first argument is the number of threads desired
    if (argc < 2) {
        fprintf(stderr, "Usage: tm_test <threads>\n");
        return EXIT_FAILURE;
    }

    threads = atol(argv[1]);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int *ranks1 = malloc((numprocs / 3)*(sizeof (int)));
    int *ranks2 = malloc((numprocs / 3)*(sizeof (int)));
    int *ranks3 = malloc((numprocs / 3)*(sizeof (int)));

    int a[SIZE];
    float tm_runtime = 0.0;
    float mutex_runtime = 0.0;
    float basic_runtime = 0.0;

    pthread_t t_id[threads];
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    float resb = 0.0, rest=0.0, resm=0.0;

    // Thread creation is attempted

    for (int i = 0; i < threads; i++) {

        if (myrank < numprocs / 3){
            if ((pthread_create(&t_id[i], &attr, &tm_run, (void*) &a)) == 1){
                fprintf(stderr, "%d: Error during tm thread creation!\n", myrank);
                return EXIT_FAILURE;
            }
        } else if (myrank < 2 * numprocs / 3) {
            if ((pthread_create(&t_id[i], &attr, &mutex_run, (void*) &a)) == 1){
                fprintf(stderr, "%d: Error during mutex thread creation!\n", myrank);
                return EXIT_FAILURE;
            }
        } else {
            if ((pthread_create(&t_id[i], &attr, &basic_run, (void*) &a)) == 1){
                fprintf(stderr, "%d: Error during basic thread creation!\n", myrank);
                return EXIT_FAILURE;
            }
        }
    }

    // Threads are joined

    for (int i = 0; i < threads; i++) {

        void *thread_ret;
        int rv = pthread_join(t_id[i], &thread_ret);

        if (rv != 0) {
            fprintf(stderr, "%d: Thread %d exited with error %d!\n", myrank, i, rv);
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

    pthread_attr_destroy(&attr);
    
    //Allreduce is run to sum the total time all ranks spent doing each computation
    
    printf("%f,%f,%f\n", rest, resm, resb);
    MPI_Allreduce(&rest, &tm_runtime, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&resm, &mutex_runtime, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&resb, &basic_runtime, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

    int b[SIZE];
    MPI_Request rrequest;
    MPI_Request srequest;
    MPI_Status status;
    //Messages are sent across the network to use bisectional bandwidth
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

    //Print a message if the arrays are different
    for (int i = 0; i < SIZE; i++) {
        if (a[i] != b[i]) {
            printf("Ranks %d and %d are different\n", myrank, (myrank + (numprocs / 2)) % numprocs);
            break;
        }
    }

    //All ranks but 0 are joined, zero prints run times and transactional memory statistics
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

void *tm_run(void *arg) {
    //variables are initialized, arg is renamed to a
    int *a = arg, b[SIZE];
    float *runtime = (float *) malloc(sizeof (float));

    unsigned long long start;
    //a and b matrices are set
    for (int i = 0; i < SIZE; i++) {
        a[i] = i;
        b[i] = -i;
    }
    //timer starts
    start = rdtsc();

    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            for (int k = 0; k < SIZE; k++) {
                //hardware transactional memory's critical region
                #pragma tm_atomic
                {
                    a[i] = a[j] + b[k];
                }

            }
        }
    }

    *runtime = (float) (rdtsc() - start) / CPU_FREQ;
    //return with total runtime
    pthread_exit(runtime);

    return (void *) runtime;
}

void *mutex_run(void *arg) {

    int *a = arg, b[SIZE];
    float *runtime = (float *) malloc(sizeof (float));
    unsigned long long start;

    for (int i = 0; i < SIZE; i++) {
        a[i] = i;
        b[i] = -i;
    }

    start = rdtsc();

    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            for (int k = 0; k < SIZE; k++) {
                //the lock is set, variables are manipulated, and lock is unset
                pthread_mutex_lock(&lock);
                a[i] = a[j] + b[k];
                pthread_mutex_unlock(&lock);

            }
        }
    }

    *runtime = (float) (rdtsc() - start) / CPU_FREQ;

    pthread_exit(runtime);

    return (void *) runtime;
}

void *basic_run(void *arg) {

    int *a = arg, b[SIZE];
    float *runtime = (float *) malloc(sizeof(float));
    unsigned long long start;

    for (int i = 0; i < SIZE; i++) {
        a[i] = i;
        b[i] = -i;
    }

    start = rdtsc();

    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            for (int k = 0; k < SIZE; k++) {
                //this is unsafe access of memory, is extremely chaotic
                a[i] = a[j] + b[k];

            }
        }
    }

    *runtime = (float) (rdtsc() - start) / CPU_FREQ;

    pthread_exit(runtime);

    return (void *) runtime;
}
