
 /*  It is called after MPI_INIT
  with the lines:

   !$OMP PARALLEL
   !$OMP CRITICAL
     call thread_bind()
   !$OMP END CRITICAL
   !$OMP END PARALLEL



   The compile line for the C routine is

  mpcc_r -qsmp -c bind.c

  */


#include <stdio.h>
#include <sys/thread.h>
#include <sys/processor.h>
#include "mpi.h"
#include <math.h>
#include <strings.h>

void get_info(int *on_proc,int *proc_index, int *MYID);
void thread_bind()
{
#define ate 4
/*
ate       = number of processors / node
nproc     = number of MPI tasks on a node
proc_id   = to which one of those MPI tasks does this thread belong,1 to
nproc
myid      = mpid task id, not used just printed for debugging
thread_id = thread id returned by omp_get_thread_num
*/

     unsigned char my_name[MPI_MAX_PROCESSOR_NAME];
     int thread_id;
     cpu_t Where;
     short int a;
     int j,n;
     int nproc,proc_id,myid;
     int omp_get_thread_num();
     thread_id=omp_get_thread_num();
/*
     get_info(&nproc,&proc_id,&myid);
*/
     a=-1;
     j=ate/(nproc);
     n=j*((proc_id)-1)+(thread_id);
     Where=n;
     bindprocessor(BINDTHREAD, thread_self(),Where);
     a=mycpu();
/*
     MPI_Get_processor_name(my_name,&j);
     printf ("mpi task %d thread number %d is running on cpu %d of node %s\n",myid,thread_id,a,my_name);
*/
}

#if 0
void get_info(int *on_proc,int *proc_index, int *MYID) {

     unsigned char *the_names,my_name[MPI_MAX_PROCESSOR_NAME];
     int i,offset,myoffset;
     int numnodes,myid,ierr;

     MPI_Comm_size(MPI_COMM_WORLD,&numnodes);

     MPI_Comm_rank(MPI_COMM_WORLD,&myid);

     the_names=(unsigned char*)malloc(numnodes*MPI_MAX_PROCESSOR_NAME*sizeof(unsigned char));

     myoffset=myid*MPI_MAX_PROCESSOR_NAME*sizeof(unsigned char);

     for (i=0;i<=8*MPI_MAX_PROCESSOR_NAME;i++)
         the_names[i]=(unsigned char)0;

     MPI_Get_processor_name(my_name,&i);

     ierr=MPI_Allgather(&my_name[0],MPI_MAX_PROCESSOR_NAME,
		MPI_UNSIGNED_CHAR,&the_names[0],MPI_MAX_PROCESSOR_NAME,MPI_UNSIGNED_CHAR,MPI_COMM_WORLD);

     *on_proc=0;

     for( i=0;i<=numnodes-1;i++){
         offset=i*MPI_MAX_PROCESSOR_NAME;
         if(strncmp(&the_names[myoffset],&the_names[offset],MPI_MAX_PROCESSOR_NAME) == 0){
             *on_proc=*on_proc+1;
             if(i == myid)*proc_index=*on_proc;
         }
     }

     *MYID=myid;
     free(the_names);
}
#endif
