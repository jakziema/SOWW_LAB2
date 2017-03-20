#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include <float.h>
//#define PRECISION 0.000001
#define RANGESIZE 1
#define DATA 0
#define RESULT 1
#define FINISH 2
#define EPSILON 0.0001
//#define DEBUG
double f(double x) {
        return sin(x)*sin(x)/x;
}
double SimpleIntegration(double a,double b,double precision) {
    double i;
    double sum=0;
    for (i=a;i<b;i+=precision)
        sum+=f(i)*precision;
    return sum;
}
int main(int argc, char **argv)
{
    int myrank,proccount;
    double a=1,b=100;
    double range[3], tempRange[3];
    double result=0,resulttemp[3];
    int sentcount=0;
    double precyzjaM=0.1;
    int i;
    MPI_Status status;
    // Initialize MPI
    MPI_Init(&argc, &argv);
    // find out my rank
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    // find out the number of processes in MPI_COMM_WORLD
    MPI_Comm_size(MPI_COMM_WORLD, &proccount);
    if (proccount<2)
    {
        printf("Run with at least 2 processes");
        MPI_Finalize();
        return -1;
    }
    if (((b-a)/RANGESIZE)<2*(proccount-1))
    {
        printf("More subranges needed");
        MPI_Finalize();
        return -1;
    }
    // now the master will distribute the data and slave processes will perform computations
    if (myrank==0)
    {
        double *tab=(double*) malloc ((proccount)*sizeof(double)); //tabela na poprzednie wartosci
        for(i=0;i<(proccount);i++)
            tab[i]=DBL_MAX; //infinity
        tab[0]=0.0;
        //for(i=0;i<(proccount);i++)
           // printf("%d:%1.16f\n",i,tab[i]);
        range[0]=a;
        range[2]=precyzjaM; //Precyzja początkowa
    // first distribute some ranges to all slaves
        for(i=1;i<proccount;i++)
        {
            range[1]=range[0]+RANGESIZE;

            #ifdef DEBUG
            printf("\nMaster sending range %f,%f precision %f to process %d",range[0],range[1],range[2],i);
            fflush(stdout);
            #endif
            // send it to process i
            MPI_Send(range,3,MPI_DOUBLE,i,DATA,MPI_COMM_WORLD);
            sentcount++;
            range[0]=range[1];
        }
        do
        {
            // distribute remaining subranges to the processes which have completed their parts
            MPI_Recv(resulttemp,3,MPI_DOUBLE,MPI_ANY_SOURCE,RESULT,MPI_COMM_WORLD,&status);
            

	    if(fabs((resulttemp[0]-tab[status.MPI_SOURCE])/resulttemp[0])>EPSILON)
            {
                #ifdef DEBUG
                printf("\nMaster NOOT result %f from process %d range %f - %f precision %f",resulttemp[0],status.MPI_SOURCE,resulttemp[2],resulttemp[2]+RANGESIZE,resulttemp[1]);
                fflush(stdout);
                #endif
                tab[status.MPI_SOURCE]=resulttemp[0];//wynik
                tempRange[0]=resulttemp[2];
                tempRange[1]=resulttemp[2]+RANGESIZE;

                tempRange[2]=resulttemp[1]*0.1;
                MPI_Send(tempRange,3,MPI_DOUBLE,status.MPI_SOURCE,DATA,MPI_COMM_WORLD);

            }
            else
            {
                //printf("Process: %d = %1.16f\n",status.MPI_SOURCE,resulttemp[0]);
                result+=resulttemp[0];
                tab[status.MPI_SOURCE]=DBL_MAX;
                #ifdef DEBUG
                printf("\nMaster approved result %f from process %d range %f - %f precision %f",resulttemp[0],status.MPI_SOURCE,resulttemp[2],resulttemp[2]+RANGESIZE,resulttemp[1]);
                fflush(stdout);
                #endif
                // check the sender and send some more data
                range[1]=range[0]+RANGESIZE;
                if (range[1]>b)
                    range[1]=b;
                #ifdef DEBUG
                printf("\nMaster sending range %f,%f to process %d",range[0],range[1],status.MPI_SOURCE);
                fflush(stdout);
                #endif
                range[2]=precyzjaM;
                MPI_Send(range,3,MPI_DOUBLE,status.MPI_SOURCE,DATA,MPI_COMM_WORLD);
                range[0]=range[1];
            }
        } while (range[1]<b);
        // now receive results from the processes
        int toContinue=1;
        double sss;
        for(;toContinue;)
        {
            MPI_Recv(&resulttemp,3,MPI_DOUBLE,MPI_ANY_SOURCE,RESULT,MPI_COMM_WORLD,&status);
            if(fabs((resulttemp[0]-tab[status.MPI_SOURCE])/resulttemp[0])>EPSILON)
            {
                #ifdef DEBUG
                printf("\nMaster NOOT result %f from process %d range %f - %f precision %f",resulttemp[0],status.MPI_SOURCE,resulttemp[2],resulttemp[2]+RANGESIZE,resulttemp[1]);
                fflush(stdout);
                #endif
                tab[status.MPI_SOURCE]=resulttemp[0];
                tempRange[0]=resulttemp[2];
                tempRange[1]=resulttemp[2]+RANGESIZE;
                tempRange[2]=resulttemp[1]*0.1;
                MPI_Send(tempRange,3,MPI_DOUBLE,status.MPI_SOURCE,DATA,MPI_COMM_WORLD);

            }
            else
            {
                #ifdef DEBUG
                printf("\nMaster END result %f from process %d range %f - %f precision %f",resulttemp[0],status.MPI_SOURCE,resulttemp[2],resulttemp[2]+RANGESIZE,resulttemp[1]);
                fflush(stdout);
                #endif
                sss=0;
                result+=resulttemp[0];
                tab[status.MPI_SOURCE]=0.0;
                for(i=0;i<(proccount);i++)
                        sss+=fabs(tab[i]);
                if(sss==0.0)
                    toContinue=0;
            }
        }
        // shut down the slaves
        for(i=1;i<proccount;i++) {
            MPI_Send(NULL,0,MPI_DOUBLE,i,FINISH,MPI_COMM_WORLD);
        }
        // now display the result
        printf("\nHi, I am process 0, the result is %f\n",result);
    }
    else
    { // slave
    // this is easy - just receive data and do the work
        do {
            MPI_Probe(0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
            if (status.MPI_TAG==DATA) {
                MPI_Recv(range,3,MPI_DOUBLE,0,DATA,MPI_COMM_WORLD,&status);
                // compute my part
                resulttemp[0]=SimpleIntegration(range[0],range[1],range[2]);
                resulttemp[1]=range[2];
                resulttemp[2]=range[0];
                // send the result back
                MPI_Send(&resulttemp,3,MPI_DOUBLE,0,RESULT,MPI_COMM_WORLD);
            }
        } while (status.MPI_TAG!=FINISH);
    }
    // Shut down MPI
    MPI_Finalize();
    return 0;
}

