#include <stddef.h>
#include <mpi.h>

int ogo_MPI_Init();
int ogo_MPI_Abort(int);
int ogo_MPI_Finalize();
int ogo_MPI_Comm_rank();
int ogo_MPI_Comm_size();
double ogo_MPI_Wtime();
int ogo_MPI_Probe_bytes(int*);
int ogo_MPI_Recv_bytes(void* buffer, int count, int source, int tag);
int ogo_MPI_Send_bytes(void* buffer, int count, int dest, int tag);
int ogo_MPI_Bcast_char(void* buffer, int count, int root);

int ogo_MPI_Init(){
    return MPI_Init(NULL, NULL);
}

int ogo_MPI_Abort(int err){
    return MPI_Abort(MPI_COMM_WORLD, err);
}

int ogo_MPI_Finalize(){
    return MPI_Finalize();
}

int ogo_MPI_Comm_rank(int* rank){
    return MPI_Comm_rank(MPI_COMM_WORLD, rank);
}

int ogo_MPI_Comm_size(int* size){
    return MPI_Comm_size(MPI_COMM_WORLD, size);
}

double ogo_MPI_Wtime(){
    return MPI_Wtime();
}

int ogo_MPI_Probe_bytes(int* data){

    const int mySource = (data[0] == -42) ? MPI_ANY_SOURCE : data[0];
    const int myTag = (data[1] == -42) ? MPI_ANY_TAG : data[1];

    MPI_Status stat;
    const int ret = MPI_Probe(mySource, myTag, MPI_COMM_WORLD, &stat);
    if (ret != MPI_SUCCESS) return ret;

    if (mySource == MPI_ANY_SOURCE) data[0] = stat.MPI_SOURCE;
    if (myTag == MPI_ANY_TAG) data[1] = stat.MPI_TAG;

    // assumes bytes
    int count = 0;
    const int ret2 = MPI_Get_count(&stat, MPI_BYTE, &count);
    data[2] = count;

    return ret2;
}

int ogo_MPI_Recv_bytes(void* buffer, int count, int source, int tag){

    MPI_Status stat;
    return MPI_Recv(buffer, count, MPI_BYTE, source, tag, MPI_COMM_WORLD, &stat);
}

int ogo_MPI_Bcast_char(void* buffer, int count, int root){
    return MPI_Bcast(buffer, count, MPI_CHAR, root, MPI_COMM_WORLD);
}

int ogo_MPI_Send_bytes(void* buffer, int count, int dest, int tag){
    return MPI_Send(buffer, count, MPI_BYTE, dest, tag, MPI_COMM_WORLD);
}
