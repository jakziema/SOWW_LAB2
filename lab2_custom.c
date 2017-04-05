#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#define PRECISION 0.000001
#define RANGESIZE 0.02
#define DATA 0
#define RESULT 1
#define FINISH 2
#define MAX 5
#define DEBUG

//mpicc lab2_custom.c -lm
//tme mpirun -np 4 ./a.out

struct Pair {
    double a;
    double b;
};
struct Pair queue[MAX];
int front = 0;
int rear = -1;
int itemCount = 0;

struct Pair peek() {
    return queue[front];
}

int isEmpty() {
    return itemCount == 0;
}

int isFull() {
    return itemCount == MAX;
}

int size() {
    return itemCount;
}

void insert(struct Pair elem) {
    if (!isFull()) {
        if (rear == MAX-1) {
            rear = -1;
        }

        queue[++rear] = elem;
        itemCount++;
    }
}

struct Pair removeData() {
    struct Pair data = queue[front++];

    if (front == MAX) {
        front = 0;
    }

    itemCount--;
    return data;
}

double f(double x) {
    return sin(x)*sin(x)/x;
}

double SimpleIntegration(double a,double b) {
    double i;
    double sum=0;
    for (i=a;i<b;i+=PRECISION)
        sum+=f(i)*PRECISION;
    return sum;
}

int main(int argc, char **argv) {
    //requests sprawdza wykonanie
    MPI_Request *requests;
    int requestcount=0;
    int requestcompleted;
    int myrank,proccount;
    double a=0,b=4*M_PI;
    double *ranges;
    double range[2];
    double result=0;
    double *resulttemp;
    int sentcount=0;
    int recvcount=0;
    int i;
    int flag;
    struct Pair data;
    int finish = 0;
    MPI_Status status;
    // Initialize MPI
    MPI_Init(&argc, &argv);
    // find out my rank
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    // find out the number of processes in MPI_COMM_WORLD
    MPI_Comm_size(MPI_COMM_WORLD, &proccount);
    if (proccount<2) {
        printf("Run with at least 2 processes");
        MPI_Finalize();
        return -1;
    }
    if (((b-a)/RANGESIZE)<2*(proccount-1)) {
        printf("More subranges needed");
        MPI_Finalize();
        return -1;
    }

    // now the master will distribute the data and slave processes will perform computations
    if (myrank==0) {
        requests=(MPI_Request *)malloc(3*(proccount-1)*sizeof(MPI_Request));
        if (!requests) {
            printf("\nNot enough memory");
            MPI_Finalize();
            return -1;
        }
        ranges=(double *)malloc(4*(proccount-1)*sizeof(double));
        if (!ranges) {
            printf("\nNot enough memory");
            MPI_Finalize();
            return -1;
        }
        resulttemp=(double *)malloc((proccount-1)*sizeof(double));
        if (!resulttemp) {
            printf("\nNot enough memory");
            MPI_Finalize();
            return -1;
        }
        range[0]=a;
        // first distribute some ranges to all slaves
        for(i=1;i<proccount;i++) {
            //wysij range o roznicy 5s
            range[1]=range[0]+(MAX*RANGESIZE);
#ifdef DEBUG
            printf("\nMaster sending range %f,%f to process %d",range[0],range[1],i);
            fflush(stdout);
#endif
            // send it to process i
            MPI_Send(range,2,MPI_DOUBLE,i,DATA,MPI_COMM_WORLD);
            sentcount += MAX;
            range[0]=range[1];
        }
        // the first proccount requests will be for receiving, the latter ones for sending
        for(i=0;i<2*(proccount-1);i++)
            requests[i]=MPI_REQUEST_NULL; // none active at this point
        // start receiving for results from the slaves
        for(i=1;i<proccount;i++)
            MPI_Irecv(&(resulttemp[i- 1]),1,MPI_DOUBLE,i,RESULT,MPI_COMM_WORLD,&(requests[i-1]));


        while (range[1]<b) {
#ifdef DEBUG
            printf("\nMaster waiting for completion of requests");
            fflush(stdout);
#endif
            // wait for completion of any of the requests
            MPI_Waitany(2*proccount-2,requests,&requestcompleted,MPI_STATUS_IGNORE);
            // if it is a result then send new data to the process
            // and add the result
            if (requestcompleted<(proccount-1)) {
                result+=resulttemp[requestcompleted];
                recvcount++;
#ifdef DEBUG
                printf("\nMaster received %d result %f from process %d",recvcount,resulttemp[requestcompleted],requestcompleted+1);
                fflush(stdout);
#endif
// first check if the send has terminated
                MPI_Wait(&(requests[proccount-1+requestcompleted]),MPI_STATUS_IGNORE);
// now send some new data portion to this process
                range[1]=range[0]+RANGESIZE;
                if (range[1]>b) range[1]=b;
#ifdef DEBUG
                printf("\nMaster sending range %f,%f to process %d",range[0],range[1],requestcompleted+1);
                fflush(stdout);
#endif
                ranges[2*requestcompleted]=range[0];
                ranges[2*requestcompleted+1]=range[1];
                MPI_Isend(&(ranges[2*requestcompleted]),2,MPI_DOUBLE,requestcompleted+1,DATA,MPI_COMM_WORLD,&(requests[proccount-1+requestcompleted]));
                sentcount++;
                range[0]=range[1];
// now issue a corresponding recv
                MPI_Irecv(&(resulttemp[requestcompleted]),1,MPI_DOUBLE,requestcompleted+1,RESULT,MPI_COMM_WORLD,&(requests[requestcompleted]));
            }
        }
// now send the FINISHING ranges to the slaves
// shut down the slaves
        range[0]=range[1];
        for(i=1;i<proccount;i++) {
#ifdef DEBUG
            printf("\nMaster sending FINISHING range %f,%f to process %d",range[0],range[1],i);
            fflush(stdout);
#endif
            ranges[2*i-4+2*proccount]=range[0];
            ranges[2*i-3+2*proccount]=range[1];
            MPI_Isend(range,2,MPI_DOUBLE,i,DATA,MPI_COMM_WORLD,&(requests[2*proccount- 3+i]));
        }
#ifdef DEBUG
        printf("\nMaster before MPI_Waitall with total proccount=%d",proccount);
        fflush(stdout);
#endif
// now receive results from the processes - that is finalize the pending requests
        MPI_Waitall(3*proccount-3,requests,MPI_STATUSES_IGNORE);
#ifdef DEBUG
        printf("\nMaster after MPI_Waitall with total proccount=%d",proccount);
        fflush(stdout);
#endif
// now simply add the results
        for(i=0;i<(proccount-1);i++) {
            result+=resulttemp[i];
            recvcount++;
        }
// now receive results for the initial sends
        for(i=0;i<4*(proccount-1);i++) {
#ifdef DEBUG
            printf("\nMaster receiving result from process %d",i%(proccount-1+1));
            fflush(stdout);
#endif
            MPI_Recv(&(resulttemp[i%(proccount-1)]),1,MPI_DOUBLE,(i%(proccount-1))+1,RESULT,MPI_COMM_WORLD,&status);
            result+=resulttemp[i%(proccount-1)];
            recvcount++;
#ifdef DEBUG
            printf("\nMaster received %d result %f from process %d",recvcount,resulttemp[i%(proccount-1)],i%(proccount-1)+1);
            fflush(stdout);
#endif
        }
// now display the result
        printf("\nHi, I am process 0, the result is %f\n",result);
    } else { // slave
        requests=(MPI_Request *)malloc(2*sizeof(MPI_Request));
        if (!requests) {
            printf("\nNot enough memory");
            MPI_Finalize();
            return -1;
        }
        requests[0]=requests[1]=MPI_REQUEST_NULL;
        ranges=(double *)malloc(2*sizeof(double));
        if (!ranges) {
            printf("\nNot enough memory");
            MPI_Finalize();
            return -1;
        }
        resulttemp=(double *)malloc(2*sizeof(double));
        if (!resulttemp) {
            printf("\nNot enough memory");
            MPI_Finalize();
            return -1;
        }
// first receive the initial data
        MPI_Recv(range,2,MPI_DOUBLE,0,DATA,MPI_COMM_WORLD,&status);
        // dzielimy na 5 rangeow
        for (i = 0; i < MAX; i++) {
            data.a = range[0];
            range[0] += RANGESIZE;
            data.b = range[0];
            insert(data);
#ifdef DEBUG
            printf("\nSlave received range %f,%f  %d", data.a, data.b, myrank);
            fflush(stdout);
#endif
        }
        //asynchroniczne odebranie rangey
        MPI_Irecv(ranges, 2, MPI_DOUBLE, 0, DATA, MPI_COMM_WORLD, &(requests[0]));
        //dopoki nie jest skonczony lub pusty
        while (!(finish && isEmpty())) {
            //jezeli nie jest pusty to pobiera dane ktore beda wyslane
            if (!isEmpty()) {
                data = removeData();
                //wait until old send was completed
                
                MPI_Wait(&(requests[1]), MPI_STATUS_IGNORE);
                //oblicz range
                resulttemp[0] = SimpleIntegration(data.a, data.b);
#ifdef DEBUG
                printf("\nSlave just computed range %f,%f  %f", data.a, data.b, resulttemp[0]);
                fflush(stdout);
#endif          
                //wyslij obliczony range
                MPI_Isend(&resulttemp[0],1,MPI_DOUBLE,0,RESULT,MPI_COMM_WORLD,&(requests[1]));
            }
            //jesli nie jest pelny
            if (!isFull()) {
                //sprawdz czy poprzedni , update range jezeli poprzedni jest wykonany
                MPI_Test(&(requests[0]), &flag, MPI_STATUS_IGNORE);
                if (flag) {
                    data.a = ranges[0];
                    data.b = ranges[1];

#ifdef DEBUG
                    printf("\nSlave just received range %f,%f  %d", data.a, data.b, myrank);
                    fflush(stdout);
#endif
                    //poza zakresem
                    if (data.a == data.b) {
                        finish = 1;
                    } else {
                        //wloz na koniec kolejki
                        insert(data);
                    }
                    
                    if (!finish) {
                        //wyslij range jezeli nie skonczono
                        MPI_Irecv(ranges, 2, MPI_DOUBLE, 0, DATA, MPI_COMM_WORLD, &(requests[0]));
                    }
                }
            }
        }

// now finish sending the last results to the master
        MPI_Wait(&(requests[1]),MPI_STATUS_IGNORE);
    }
// Shut down MPI
    MPI_Finalize();
#ifdef DEBUG
    printf("\nProcess %d finished",myrank);
    fflush(stdout);
#endif
    return 0;
}
