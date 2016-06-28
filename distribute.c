/*
 *  distribute.c
 *  
 *
 *  Created by Johannes Langguth on 26.08.14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */
#include <mpi.h>
#include "commonDefs.h"
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h> 
#include <assert.h> 
#include <string.h>
#include "utility.h" 

//utility
int sddMPI_inlist(int element,int * list,int listsize)
{
	for (int i= 0; i < listsize; i++)	
	{
		if(element == list[i])
		{
			return(i);
		}
	}
	return(-1);
}

//utility
int sddMPI_getpart(int i,localparttype* L)
{
	int k=0;
	while(i >= L->limits[k+1])
		k++;
	assert(k <= L->parts);
	return(k);
}


//all MPI ranks
void sddMPI_localRenumbering(ELLmatrixVector* LZ, localparttype* L)
{	
		 printf("sepcount %d ",L->sepcount);
	if(L->limits[L->myrank]>0)						//renumbering of local tetrahedra
	for (int i= L->sepcount; i < L->mysize; i++)	
	{
		for (int k=0;  k<RNZ;  k++ ) 
		{
			//LZ->I[i*16+k]=LZ->I[i*16+k]-L->limits[L->myrank];
			if(LZ->I[i*16+k]<0)
				printf(" at %d step vertex %d k %d limits %d idx %d A %e \n ",L->myrank,i,k,L->limits[L->myrank],LZ->I[i*16+k],LZ->A[i*16+k]);
		}
	}
	 
	//setup receive requests
	int* remotelist = (int*)calloc(RNZ*L->sepcount,sizeof(int));
	int* tmpreadlength = (int*)calloc(L->parts+1,sizeof(int));
	
	int** tmpremotelists = (int**)calloc(L->parts,sizeof(int*));
	for (int i= 0; i < L->parts; i++)	
	{
		tmpremotelists[i] = (int*)calloc(RNZ*L->sepcount,sizeof(int));		//potentially too large
	}
	
	for (int i= 0; i < L->sepcount; i++)	
	{
		for (int k=0;  k<RNZ;  k++ )											
		{
			int elem = LZ->I[i*16+k]+L->limits[L->myrank];
			//int part = L->getpart[i*16+k];
			int part = sddMPI_getpart(elem,L);

		 printf("\n  %d %d ",elem,part);

			if(part!=L->myrank)							//off proc tetrahedron neighbour.
			{	
				int l = sddMPI_inlist(elem,tmpremotelists[part],tmpreadlength[part]);
				if(l == -1)	    //check if already encountered 
				{
					tmpremotelists[part][tmpreadlength[part]]=elem;
					if(elem >= L->limits[part+1])
						printf(" on rank %d elem %d part %d limits %d \n ",L->myrank,elem , part, L->limits[part+1]);
					if(elem < L->limits[part])
						printf(" on rank %d elem %d part %d limits %d \n ",L->myrank,elem , part, L->limits[part]);
					LZ->I[i*16+k]=tmpreadlength[part]+L->mysize;
					tmpreadlength[part]++;
				}
				else {
					LZ->I[i*16+k]=l+L->mysize;
				}
			}
			else {
			//	LZ->I[i*16+k]=elem-L->limits[L->myrank];				//on proc tetrahedron. Normal renumbering.
			}
		}
	}
	
	L->remoteVcount=0;
	for (int i= 0; i < L->parts; i++)	
	{
		L->remoteVcount+=tmpreadlength[i];  
	}
	
	
	L->recvneighbourcount=0;													//count the number of nodes to receive from, reallocate remote lists
	for (int i= 0; i < L->parts; i++)	
	{
		if(tmpreadlength[i]>0)
		{
			L->recvneighbourcount++;
		}
	}	
	
	L->recvecounts = (int*)calloc(L->recvneighbourcount,sizeof(int));
	L->recvneighbours = (int*)calloc(L->recvneighbourcount,sizeof(int));
	
	int tmp=0;																	//for each node to receive from, store number of elements to receive
	for (int i= 0; i < L->parts; i++)	
	{
		if(tmpreadlength[i]>0)
		{
			L->recvecounts[tmp]=tmpreadlength[i];
			L->recvneighbours[tmp]=i;											//for each node to receive from, store number of nodes in recvneighbours contigous list
			tmp++;
		}
	}	
	
	//setup send structure
	
	
	int* tmpsendlist = (int*)calloc(L->parts,sizeof(int));		
	MPI_Alltoall(tmpreadlength, 1, MPI_INT,tmpsendlist, 1, MPI_INT,MPI_COMM_WORLD);		//communicate expected receive numbers to sender
	
	L->sendneighbourcount=0;													//count the number of nodes to send to
	for (int i= 0; i < L->parts; i++)	
	{
		if(tmpsendlist[i]>0)
			L->sendneighbourcount++;
	}	
	
	L->sendcounts = (int*)calloc(L->sendneighbourcount,sizeof(int));
	L->sendneighbours = (int*)calloc(L->sendneighbourcount,sizeof(int));	
	
	tmp=0;
	for (int i= 0; i < L->parts; i++)											//for each node to send to, store number of elements to send
	{	
		if(tmpsendlist[i]>0)
		{
			L->sendcounts[tmp]=tmpsendlist[i];
			L->sendneighbours[tmp]=i;											//for each node to send to, store number of nodes in sendneighbours contigous list
			tmp++;
		}
	}

	//build send lists
	
	MPI_Request* ar_send_req=(MPI_Request*)malloc(L->recvneighbourcount*sizeof(MPI_Request));
	MPI_Status* ar_status=(MPI_Status*)malloc(L->sendneighbourcount*sizeof(MPI_Status));
	MPI_Request* ar_recv_req=(MPI_Request*)malloc(L->sendneighbourcount*sizeof(MPI_Request));
	
	L->sendlists = (int**)calloc(L->sendneighbourcount,sizeof(int*));
	L->sendbuffer = (double**)calloc(L->sendneighbourcount,sizeof(double*));
	
	//int sum=0;
	for (int i = 0; i < L->recvneighbourcount; i ++) {							//send list of requested numbers to future sender
		MPI_Isend(tmpremotelists[L->recvneighbours[i]], L->recvecounts[i], MPI_INT, L->recvneighbours[i], 0, MPI_COMM_WORLD, &ar_send_req[i]);
	}
	
	for (int i = 0; i < L->sendneighbourcount; i ++) {							//receive sendlists
		L->sendlists[i] = (int*)calloc(L->sendcounts[i],sizeof(int));
		L->sendbuffer[i] = (double*)calloc(L->sendcounts[i],sizeof(double));
		MPI_Irecv(L->sendlists[i], L->sendcounts[i], MPI_INT, L->sendneighbours[i], 0, MPI_COMM_WORLD, &ar_recv_req[i]);
	}
	
	if(L->sendneighbourcount>0)
		MPI_Waitall( L->sendneighbourcount, ar_recv_req, ar_status);

	//last step: renumber entries in seperators and sendlists

	for (int i= 0; i < L->sendneighbourcount; i++)								//change entries in sendlists from global to local renumbering
	{
		for (int j=0;  j<L->sendcounts[i];  j++ ) 
		{	
			//printf("  on %d sendlists[%d][%d] = %d  limit %d \n ",L->myrank,i,j, L->sendlists[i][j], L->limits[L->myrank]);
			L->sendlists[i][j]=L->sendlists[i][j]-L->limits[L->myrank];			
		}
	}
	
	L->recvstart = (int*)calloc(L->parts,sizeof(int));
	for (int j=1;  j<L->parts ;  j++ ) {
		L->recvstart[j]=L->recvstart[j-1]+L->recvecounts[j-1];
	}
	
	count2sum(tmpreadlength,L->parts);
	L->remoteVcount=tmpreadlength[L->parts];
	for (int i= 0; i < L->sepcount; i++)	
	{
		for (int k=0;  k<RNZ;  k++ )											
		{
			int elem = LZ->I[i*16+k];
			//int part = L->getpart[i*16+k];
			int part = sddMPI_getpart(elem,L);
			if(part!=L->myrank)							//off proc tetrahedron neighbour.
			{	
				LZ->I[i*16+k]=elem+tmpreadlength[part];
			}
		}
	}
	

	free(remotelist);
	free(tmpreadlength);
	free(tmpsendlist);
}




void sdd_Retreive(ELLmatrixVector* Z,ELLmatrixVector* LZ,localparttype* L)
{
	
	int* scattercount = (int*)calloc(L->parts,sizeof(int));

	for (int j = 0; j < L->parts; j++)					{
		scattercount[j]=L->limits[j+1]-L->limits[j];
	}

	MPI_CHECK(MPI_Gatherv(LZ->V, L->mysize, MPI_DOUBLE, Z->V, scattercount, L->limits, MPI_DOUBLE, 0, MPI_COMM_WORLD));	
}


	
