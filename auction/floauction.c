/**************************************************************
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "Madan Sathe, Ye Zhao, University of Basel, Switzerland" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

For any kind of questions or comments please send an email to Madan Sathe: madan.sathe@unibas.ch

****************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <omp.h>
#include <assert.h>
#include "floauction.h"

#define INFINITY_ 1000000000
#define log2(x) (log(x)/log(2))

int MAX_ITER = 1000000000;
int ITER_AUCTION = 100;
int ITER_ALL = 0;
int INFO = 0;
double ERR = 0.0;
double GMIN = 0.000001;
double GMAX = 2;
double SCALE = 1.0;
double BEGIN_EPS = 1.0;//1 - 1000
double END_EPS = 0.00001;// 4 
//double LIMIT = 0.000000000001; //--> dense problems
double LIMIT = 1e-5; //--> dense problems
//double LIMIT = 0.0; 
//double LIMIT = 0.00000001;
int FACTOR = 2; // 4
int INCR_FACTOR = 2; // 4
int THRESH_MAX = 100;//100
int THETA = 12; //40
int OMEGA = 16; //200
double APPROX_FACTOR = 1.0;

int isMatrixDistributed = 0;
int isBoundaryInfoDistributed = 0;
int isMatrixDense = 1;
int returnDualVariables = 0;
int useEpsilonScaling = 1;
int areValuesScaled = 0;

double TIME_DISTRIBUTION = 0;
double TIME_WEIGHT_SCALING = 0;
double TIME_GET_BEST_OBJECT = 0;
double TIME_BID = 0;
double TIME_UPDATE = 0;
double TIME_AUCTION = 0;
double TIME_MERGE_MATCHING = 0;
double TIME_ADD_MATCHING = 0;
double TIME_ALL = 0.0;
double TIME_COMMUNICATION_MERGE = 0;
double TIME_COMMUNICATION_UPDATE = 0;
double TIME_PARTITION = 0.0;
double TIME_PERMUTE = 0.0;
double TIME_SCALING = 0.0;
long int MSGS = 0;
			

int mergeAndUpdatePrice(double * changed_idx_price_local, double ** changed_idx_price_all, int nChanged, int * changed_displs, MPIType mpiVar){

	double starttime1 = 0.0, stoptime1 = 0.0;
	int i = 0, j = 0;
	int nChanged_all = 0;
	int * recvcounts = NULL;
	MPI_Request  * request1 = NULL,  * request2 = NULL;
	MPI_Status * status = NULL;	
	
	starttime1 = MPI_Wtime();
			
	nChanged = 2 * nChanged;

	//MPI_Gather(&nChanged, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	if(mpiVar.mype != 0){
		request1 = (MPI_Request *) malloc(sizeof(MPI_Request));
		request2 = (MPI_Request *) malloc(sizeof(MPI_Request));
		MPI_Isend(&nChanged, 1, MPI_INT, 0, mpiVar.mype, MPI_COMM_WORLD, request1);
		MPI_Isend(&changed_idx_price_local[0], nChanged, MPI_DOUBLE, 0, mpiVar.npes+mpiVar.mype, MPI_COMM_WORLD, request2);
	}
	else{		
		recvcounts = (int *) malloc(mpiVar.npes * sizeof(int));
		status = (MPI_Status *) malloc( (mpiVar.npes-1) * sizeof(MPI_Status));
		request1 = (MPI_Request *) malloc( (mpiVar.npes-1) * sizeof(MPI_Request));

		for(i = 1; i < mpiVar.npes; i++){
			MPI_Irecv(&recvcounts[i], 1, MPI_INT, i, i, MPI_COMM_WORLD, &request1[i-1]);
		}
		
		recvcounts[0] = nChanged;
		memset(status, 0, (mpiVar.npes-1) * sizeof(MPI_Status));
		MPI_Waitall((mpiVar.npes-1), request1, status);

		for(i=0;i<mpiVar.npes;i++){
			changed_displs[i] = nChanged_all;
			nChanged_all += recvcounts[i];
		}
	
		*changed_idx_price_all = (double *) malloc(nChanged_all * sizeof(double));

		for(i = 1; i < mpiVar.npes; i++){
			MPI_Irecv(&(*changed_idx_price_all)[changed_displs[i]], recvcounts[i], MPI_DOUBLE, i, mpiVar.npes+i, MPI_COMM_WORLD, &request1[i-1]);
		}

		for(j = 0; j < nChanged; j++) (*changed_idx_price_all)[j] = changed_idx_price_local[j];
		
		memset(status, 0, (mpiVar.npes-1) * sizeof(MPI_Status));
		MPI_Waitall((mpiVar.npes-1), request1, status);
		
		free(status);
		free(request1);
		free(recvcounts);
	}
		
	MPI_Bcast(&nChanged_all, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(mpiVar.mype != 0) {
		*changed_idx_price_all = (double *) malloc(nChanged_all * sizeof(double));	
		free(request1);
		free(request2);
	}
	MPI_Bcast(&(*changed_idx_price_all)[0], nChanged_all, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	

	
	stoptime1 = MPI_Wtime();
	TIME_COMMUNICATION_UPDATE =  TIME_COMMUNICATION_UPDATE + (stoptime1 - starttime1);

	MSGS = MSGS + nChanged_all;

	return nChanged_all;
}

int mergeMatching(int * row, int * col, int * perm, double * scaling, int neqns_total, int * ia, int * ja, double * a, int neqns, int * displs, MPIType mpiVar, double * amax, double * price){
	double starttime = 0.0, stoptime = 0.0, starttime1 = 0.0, stoptime1 = 0.0;
	starttime = MPI_Wtime();
	int i = 0, j = 0, l = 0;
	int nMatch = 0;	
	int idx_object = -1;
	
	int * recvcounts = NULL;	
	double scaled_max = 0.0, val = 0.0, value = 0.0;

	double * u = NULL;	
	if(returnDualVariables == 1){
		u = malloc(neqns * sizeof(double));
	}
	
	MPI_Request * request = NULL;
	MPI_Status * status = NULL;
	
	if(mpiVar.mype == 0){
		recvcounts = malloc(mpiVar.npes * sizeof(int));		
		
		for(i=0;i<mpiVar.npes;i++){
			recvcounts[i] = displs[i+1] - displs[i];
		}	
	}

	//MPI_Gatherv(&row[displs[mpiVar.mype]], neqns, MPI_INT, perm, recvcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);

	value = GMIN/GMAX;
	
	starttime1 = MPI_Wtime();
	
	if(mpiVar.mype != 0){
		request = (MPI_Request *) malloc(2*sizeof(MPI_Request));
		MPI_Isend(&row[displs[mpiVar.mype]], neqns, MPI_INT, 0, mpiVar.mype, MPI_COMM_WORLD, &request[0]);
		
		if(returnDualVariables == 1){	
			#pragma omp parallel for private(i, j, val, scaled_max)
			for(i=0;i<neqns;i++){
				scaled_max = 0.0;
				for(j=ia[i];j<ia[i+1];j++){
					val = a[j] - price[ja[j]];
					if(val > scaled_max){
						scaled_max = val;
					}
				}
				if(scaled_max != 0.0){
					u[i] = scaled_max;
					u[i] = exp(u[i] / SCALE);
					u[i] = 1.0 / (u[i] * amax[i] * value);
				}
				else{
					u[i] = 1.0;
				}
			}
		}
	}else{
		request = (MPI_Request *) malloc((mpiVar.npes-1) * sizeof(MPI_Request));
		for(i = 1; i < mpiVar.npes; i++){
			MPI_Irecv(&perm[displs[i]], recvcounts[i], MPI_INT, i, i, MPI_COMM_WORLD, &request[i-1]);
		}
		#pragma omp parallel for private(l)		
		for(l = 0; l < neqns; l++) {
			perm[l] = row[l];
		}
		if(returnDualVariables == 1){
			#pragma omp parallel for private(i, j, val, scaled_max)
			for(i = 0; i < neqns; i++){
				scaled_max = 0.0;
				for(j = ia[i]; j < ia[i+1]; j++){
					val = a[j] - price[ja[j]];
					if(val > scaled_max){
						scaled_max = val;
					}
				}
				if(scaled_max != 0.0){
					u[i] = scaled_max;
					u[i] = exp(u[i] / SCALE);
					u[i] = 1.0 / (u[i] * amax[i] * value);
				}else{
					u[i] = 1.0;
				}
			}
		}
	}
	if(mpiVar.mype == 0){
		status = (MPI_Status *) malloc((mpiVar.npes-1) * sizeof(MPI_Status));
		MPI_Waitall((mpiVar.npes-1), request, status);
	}
	
	if(returnDualVariables == 1){
		if(mpiVar.mype != 0){
			MPI_Isend(u, neqns, MPI_DOUBLE, 0, mpiVar.mype, MPI_COMM_WORLD, &request[1]);
		}else{
			for(i = 1; i < mpiVar.npes; i++){		
				MPI_Irecv(&scaling[displs[i]], recvcounts[i], MPI_DOUBLE, i, i, MPI_COMM_WORLD, &request[i-1]);
			}
			#pragma omp parallel for private(j)
			for(j = 0; j < neqns; j++){
		 		scaling[j] = u[j];
			}

			#pragma omp parallel for private(j)
			for(j = 0; j < neqns_total; j++){
				scaling[neqns_total+j] = 1.0 / (exp(price[j] / SCALE));
			}
		
			MPI_Waitall((mpiVar.npes-1), request, status);		
		}	
	}

	
	stoptime1 = MPI_Wtime();		
	TIME_COMMUNICATION_MERGE = TIME_COMMUNICATION_MERGE + (stoptime1 - starttime1);

	if(mpiVar.mype == 0){		
		//#pragma omp parallel for private (i, idx_object) reduction(+:nMatch)
		for(i = 0; i < neqns_total; i++){
			idx_object = perm[i];

			if(idx_object > -1 && row[i] == -1 && col[idx_object] == -1){
				row[i] = idx_object;
				col[idx_object] = i;
			}
			perm[i] = row[i];
			
			if(perm[i] > -1){
				nMatch++;
			}
		}
				
		free(status);
		free(recvcounts);
	}

	MPI_Bcast(&nMatch, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	free(request);
	
	if(returnDualVariables == 1){
		free(u);
	}
	
	stoptime = MPI_Wtime();
	TIME_MERGE_MATCHING = TIME_MERGE_MATCHING + (stoptime - starttime);
	
	return nMatch;
}

int AuctionAlgorithmParallel( int * ia, int * ja, double * a, int neqns, int * displs, int * row, int * col, int * perm, double * scaling, double * amax, int neqns_total, MPIType mpiVar){
	double starttime = 0.0, stoptime = 0.0, starttime0 = 0.0, stoptime0 = 0.0, starttime2 = 0.0, stoptime2 = 0.0;
	starttime0 = MPI_Wtime();				
	int i = 0, j = 0;
	int nMatch = 0, nChanged = 0, count = 0, thresh = 0, nolist = 0, nolist_all = 0;	
	int idx_person = -1, idx_object = -1, idx_best_object = -1;
	double epsilon = 0.0, start_incr = 0.0, incr = 0.0, profit = 0.0;
	double price_object = 0.0, inv_scale = 0.0;	
	int index = -1;
	int * nonewlist = (int *) malloc(sizeof(int));
	*nonewlist = 0;
	int * list = (int *) malloc(neqns * sizeof(int));
	int * listBestObject = (int *) malloc(neqns * sizeof(int));
	double * listMaxProfit = (double *) malloc(neqns * sizeof(double));
	double * listSecMaxProfit = (double *) malloc(neqns * sizeof(double));
	double * listCost = (double *) malloc(neqns * sizeof(double));
	double * listBid = (double *) malloc(neqns * sizeof(double));
	
	double * price = (double *) malloc(neqns_total * sizeof(double));
        double * tmp_bid = (double *) malloc(neqns_total * sizeof(double));

	double * changed_idx_price_local = (double *) malloc(2*neqns * sizeof(double));
	double * changed_idx_price_all;
	
	int * changed_displs = (int *) malloc(mpiVar.npes * sizeof(int));

	int cnt;

	double scal_val = 1.0 / SCALE;
	inv_scale = scal_val < LIMIT ? scal_val : LIMIT;

	BEGIN_EPS = SCALE / THETA;
	if(BEGIN_EPS < END_EPS) BEGIN_EPS = END_EPS;
	epsilon = BEGIN_EPS;
	
	if(isMatrixDense == 1){
                BEGIN_EPS = 1.0 / BEGIN_EPS; 
		start_incr = END_EPS; //> BEGIN_EPS ? BEGIN_EPS : END_EPS;//for dense matrix
	}else{
		start_incr = BEGIN_EPS;//for sparse matrix
	}
			
	#pragma omp parallel for private (i)
	for(i = 0; i < neqns_total; i++){
		price[i] = 0.0;
                tmp_bid[i] = 0.0;
	}
	for(i = 0; i < neqns; i++){
		idx_person = i + displs[mpiVar.mype];
		if(idx_person < 0) {
			printf("Error: Index of %d is smaller than 0 on proc %d.\n", idx_person, mpiVar.mype);
			MPI_Abort(MPI_COMM_WORLD, i);
		}
		if(row[idx_person] == -1){
			list[*nonewlist] = idx_person;
			nonewlist[0]++;
		}
	}
	nolist = *nonewlist;

	if(useEpsilonScaling == 1){	
		MPI_Allreduce (&nolist, &nolist_all, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		thresh = nolist_all / FACTOR < neqns_total / OMEGA ? nolist_all / FACTOR : neqns_total / OMEGA;
	}
	if(INFO > 2) printf("epsilon start value %f\n", start_incr);
	//START AUCTION
	incr = start_incr;
	do{
		if(epsilon <= END_EPS){
			thresh = 0;
		}
		if(useEpsilonScaling == 1) incr = start_incr;
		if(incr > epsilon){
			incr = epsilon;
		}
		do{
			starttime2 = 0.0; starttime = 0.0; stoptime = 0.0; stoptime2 = 0.0;
			starttime = MPI_Wtime();	
			starttime2 = MPI_Wtime();
			cnt = 0;

			nChanged = 0;
			*nonewlist = 0;
		
			#pragma omp parallel for schedule(static) private (i, j, idx_object, profit, idx_person, index)
			for(i = 0 ; i < nolist; i++){
			
				listBestObject[i] = -1;
				listMaxProfit[i] = 0.0;
				listSecMaxProfit[i] = 0.0;
				listCost[i] = 0.0;
				
				idx_person = list[i];
				idx_object = -1;
				index = idx_person - displs[mpiVar.mype];
				if(ia[index+1] - ia[index] == 1){
					if(a[ia[index]] > 0.0){
						listBestObject[i] = ja[ia[index]];
						listCost[i] = a[ia[index]];
					}
				}else{
					// get best_object, max_profit and sec_max_profit
					for(j = ia[index]; j < ia[index+1]; j++){
						if(a[j] > 0.0){
							idx_object = ja[j];
							profit = a[j] - price[idx_object];
							if(profit > listSecMaxProfit[i]){
								if(profit > listMaxProfit[i]){
									listSecMaxProfit[i] = listMaxProfit[i];
									listMaxProfit[i] = profit;
									listBestObject[i] = idx_object;
									listCost[i] = a[j];
								}else{
									listSecMaxProfit[i] = profit;
								}
							}
						}
					}
				}
				// calculate the bid for given person: bid = price_j + max_profit_i - sec_max_profit + incr
				if(listBestObject[i] > -1) {
				    //incr = 0.01;
                        	        tmp_bid[listBestObject[i]] = 0.0;
					listBid[i] = price[listBestObject[i]] + listMaxProfit[i] - listSecMaxProfit[i] + incr;
					if(INFO > 3) printf("bid of person %d for object %d with bid %lf (1st %lf, 2nd %lf, incr %lf) \n", list[i], listBestObject[i], listBid[i], listMaxProfit[i], listSecMaxProfit[i], incr);
				}
			}
			stoptime2 = MPI_Wtime();		
				//MPI_Abort(MPI_COMM_WORLD, i);
			TIME_GET_BEST_OBJECT = TIME_GET_BEST_OBJECT + (stoptime2 - starttime2);
			////////////////////////////////////////////////////////////////////////////////////
			// assign the object to the person and set the price with bid, the old owner left this object
			for(i = 0; i < nolist; i++){
  		           if(listBestObject[i] > -1){
                       	       if(listBid[i] > tmp_bid[listBestObject[i]])
                        	  tmp_bid[listBestObject[i]] = listBid[i];
		           }
            		}



			for(i = 0; i < nolist; i++){
			
				if(listBestObject[i] > -1){
					idx_best_object = listBestObject[i];
					idx_person = list[i];	
					
					//if(listBid[i] > price[idx_best_object]){			
					if(listBid[i] == tmp_bid[idx_best_object]){			
						price[idx_best_object] = listBid[i];
						row[idx_person] = idx_best_object;

						if(col[idx_best_object] > -1){
							row[col[idx_best_object]] = -1;
							if(col[idx_best_object] >= displs[mpiVar.mype] && col[idx_best_object] < displs[mpiVar.mype+1]){
								list[*nonewlist] = col[idx_best_object];
								nonewlist[0]++;
							}
						}
					
						col[idx_best_object] = idx_person;

						changed_idx_price_local[2*nChanged] = idx_best_object;
						changed_idx_price_local[2*nChanged+1] = price[idx_best_object];
						
						nChanged++;
						
						if(INFO > 3) printf("assigned is person %d to object %d with bid %lf\n", list[i], listBestObject[i], listBid[i]);
					}
					else{
						list[*nonewlist] = idx_person;
						nonewlist[0]++;
					
					}
				}
			}		
			stoptime = MPI_Wtime();
			TIME_BID = TIME_BID  + (stoptime - starttime);
			count++;

			if(INFO > 2) printf("proc(%d) PAA (bid) %d time %f: nChanged %d Matching %d/%d, eps: %.15e\n", mpiVar.mype, count, stoptime - starttime, nChanged, neqns - *nonewlist, neqns, incr);	
			
			if(useEpsilonScaling == 0){	
					if(incr > 0.0) incr = incr - inv_scale;
				else{
					incr = incr + inv_scale;
				}
			}

			/////////////////////////////////////////////////////////////////////////////////////////////////////////
			// merge changed price and update local price
			if(mpiVar.npes > 1){
				starttime = MPI_Wtime();	
				nChanged = mergeAndUpdatePrice(changed_idx_price_local, &changed_idx_price_all, nChanged, changed_displs, mpiVar);
				
				for(i = 0; i < nChanged * 0.5; i++){			
					idx_object = changed_idx_price_all[2*i];
					price_object = changed_idx_price_all[2*i+1];						

					if(price_object > price[idx_object]){
						price[idx_object] = price_object;
						idx_person = col[idx_object];

						if(idx_person > -1){
							if(idx_person >= displs[mpiVar.mype] && idx_person < displs[mpiVar.mype+1]){
								list[*nonewlist] = idx_person;
						
								nonewlist[0]++;
							}
							row[idx_person] = -1;
							col[idx_object] = -1;
						}
					}
				}
				free(changed_idx_price_all);	
				
				stoptime = MPI_Wtime();	
				TIME_UPDATE = TIME_UPDATE + (stoptime - starttime);
				
				if(INFO > 2) printf("proc(%d) PAA (update) %d time %f: Matching %d/%d\n", mpiVar.mype, count, stoptime - starttime, neqns - *nonewlist, neqns);	
			}
					
			nolist = *nonewlist;
			
			if(useEpsilonScaling == 1){	
				MPI_Allreduce (&nolist, &nolist_all, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
				incr = incr * INCR_FACTOR;
				if(incr > epsilon){
					incr = epsilon;
				}
				// REDUCE EPSILON
				if(nolist_all <= thresh){			  	
					thresh = thresh / FACTOR ;
					if( epsilon > END_EPS){
						epsilon = epsilon / FACTOR;
						if(epsilon > incr){
							epsilon = epsilon / FACTOR;
						}
						if(epsilon < END_EPS){
							epsilon = END_EPS;
						}
					}
					if(start_incr < epsilon){
						start_incr =  start_incr * FACTOR;
						break;
					}
				}
			}
			
			cnt++;

		}while(nolist_all > thresh && count < MAX_ITER);

	}while(nChanged > 0 && count < MAX_ITER);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// merge matching
	nMatch = neqns - nolist;		
	
	nMatch = mergeMatching(row, col, perm, scaling, neqns_total, ia, ja, a, neqns, displs, mpiVar, amax, price);
		
	free(list);
	free(nonewlist);	
	free(listBestObject);
	free(listMaxProfit);
	free(listSecMaxProfit);
	free(listCost);
	free(listBid);
	free(price);
        free(tmp_bid);
	free(changed_idx_price_local);
	free(changed_displs);

	
	stoptime0 = MPI_Wtime();
	TIME_AUCTION = TIME_AUCTION + (stoptime0 - starttime0);
	ITER_ALL = ITER_ALL + count;
	
	if(mpiVar.mype == 0 && INFO > 1) printf("PAA: Matching %d/%d, # iter %d\n", nMatch, neqns_total, count);	
	
	return nMatch;
}

AdjMatrixType distributeVtxs(int neqns, int * ia, int * ja, double * a, int * displs, MPIType mpiVar){
	double starttime = 0.0, stoptime = 0.0;
	starttime = MPI_Wtime();	

	MPI_Status status;
	MPI_Status * status_send = NULL, * status_recv;
	MPI_Request * request_send = NULL, * request_recv;
	
	int * rBlockSize;
	int * procs_nnz,  * procs_nnz_bndry = NULL;
	int * p_ia, * p_ja;
	double * p_a;
	
	int i;
	int nVtxs, nEdges;
	
	if(mpiVar.mype == 0){
		procs_nnz = malloc((mpiVar.npes+1) * sizeof(int));
		procs_nnz_bndry = malloc((mpiVar.npes+1) * sizeof(int));		
		rBlockSize = malloc(mpiVar.npes * sizeof(int));
		request_send = malloc(4*(mpiVar.npes-1)*sizeof(MPI_Request));		
		status_send = malloc(4*(mpiVar.npes-1)*sizeof(MPI_Status));		
		
		procs_nnz_bndry[0] = 0;
		#pragma omp parallel for private(i)
		for(i = 0; i < mpiVar.npes; i++){
			rBlockSize[i] = displs[i+1] - displs[i] + 1;
			procs_nnz[i] = ia[displs[i+1]] - ia[displs[i]];
		}
		for(i = 0; i < mpiVar.npes; i++){
			procs_nnz_bndry[i+1] = procs_nnz_bndry[i] + procs_nnz[i];
		}
	}
	
	nVtxs = displs[mpiVar.mype + 1] - displs[mpiVar.mype];
	p_ia = malloc((nVtxs+1) * sizeof(int));
	
	if(mpiVar.mype == 0){
		memset(request_send, 0, 4*(mpiVar.npes-1)*sizeof(MPI_Request));
		for(i = 1; i < mpiVar.npes; i++){
			MPI_Isend (&procs_nnz[i], 1, MPI_INT, i, 4*(i-1), MPI_COMM_WORLD, &request_send[4*(i-1)]);
			MPI_Isend (&ia[displs[i]], rBlockSize[i], MPI_INT, i, 4*(i-1)+1, MPI_COMM_WORLD, &request_send[4*(i-1)+1]);
			MPI_Isend (&ja[procs_nnz_bndry[i]], procs_nnz[i], MPI_INT, i, 4*(i-1)+2, MPI_COMM_WORLD, &request_send[4*(i-1)+2]);
			MPI_Isend (&a[procs_nnz_bndry[i]], procs_nnz[i], MPI_DOUBLE, i, 4*(i-1)+3, MPI_COMM_WORLD, &request_send[4*(i-1)+3]);
		}
		nEdges = procs_nnz[0];
		
		p_ja = malloc(nEdges * sizeof(int));
		p_a = malloc(nEdges * sizeof(double));
		
		memcpy(p_ia, ia, (nVtxs+1) * sizeof(int));
		memcpy(p_ja, ja, nEdges * sizeof(int));
		memcpy(p_a, a, nEdges * sizeof(double));
		
		int begin = p_ia[0];
		if(begin == 1){
			#pragma omp parallel for private (i) 
			for(i=0;i<=nVtxs;i++){
				p_ia[i] = p_ia[i] - begin;
			}
		}
		
		memset(status_send, 0, 4*(mpiVar.npes-1)*sizeof(MPI_Status));
		MPI_Waitall(4*(mpiVar.npes-1), request_send, status_send);
		
		free(procs_nnz);
		free(procs_nnz_bndry);
		free(rBlockSize);
		free(status_send);
		free(request_send);
		
	}else{
		memset(&status, 0, sizeof(MPI_Status));
		MPI_Recv(&nEdges, 1, MPI_INT, 0, 4*(mpiVar.mype-1), MPI_COMM_WORLD, &status);
		p_ja = malloc(nEdges * sizeof(int));
		p_a = malloc(nEdges * sizeof(double));
		request_recv = malloc(3 * sizeof(MPI_Request));
		status_recv = malloc(3 * sizeof(MPI_Status));
		
		MPI_Irecv (p_ia, nVtxs+1, MPI_INT, 0, 4*(mpiVar.mype-1)+1, MPI_COMM_WORLD, &request_recv[0]);
		MPI_Irecv (p_ja, nEdges, MPI_INT, 0, 4*(mpiVar.mype-1)+2, MPI_COMM_WORLD, &request_recv[1]);
		MPI_Irecv (p_a, nEdges, MPI_DOUBLE, 0, 4*(mpiVar.mype-1)+3, MPI_COMM_WORLD, &request_recv[2]);
		
		memset(status_recv, 0, 3*sizeof(MPI_Status));
		MPI_Waitall(3, request_recv, status_recv);
		
		int begin = p_ia[0];
		#pragma omp parallel for private (i) 
		for(i=0;i<=nVtxs;i++){
			p_ia[i] = p_ia[i] - begin;
		}
		free(request_recv);
		free(status_recv);
	}
	// printf("proc %d TIME_DISTRIBUTION %f\n", mpiVar.mype, TIME_DISTRIBUTION);
	if(INFO > 2) printf("proc %d neqns %d nnz %d\n", mpiVar.mype, nVtxs, nEdges);	
	
	AdjMatrixType adjMatrix;
	adjMatrix.ia = p_ia;
	adjMatrix.ja = p_ja;
	adjMatrix.a = p_a;
	adjMatrix.neqns = nVtxs;
	adjMatrix.nnz = nEdges;
	
	stoptime = MPI_Wtime();	
	TIME_DISTRIBUTION =  stoptime - starttime;

	return adjMatrix;
}

int scaleWeight(int * ia, int * ja, double * a, double * amax, int neqns, int nnz){
	double starttime = 0.0, stoptime = 0.0;
	starttime = MPI_Wtime();	
	int i = 0, j = 0;
	double min = INFINITY_, minShared = INFINITY_;
	double max = 0.0, maxShared = 0.0;
	double scaled_max = 0.0;

	if(areValuesScaled == 0){
		#pragma omp parallel private (i,j) shared(minShared, maxShared) firstprivate(min, max)
		{
			#pragma omp for nowait
			for(i=0;i<neqns;i++){
				amax[i] = 0.0;
				for(j=ia[i];j<ia[i+1];j++){
					if(a[j]!=0.0){
						amax[i] = fabs(a[j]) > amax[i] ? fabs(a[j]) : amax[i];
						if(fabs(a[j]) < min){
							min = fabs(a[j]) < min ? fabs(a[j]) : min;
						}
						if(fabs(a[j]) > max){
							max = fabs(a[j]) > max ? fabs(a[j]) : max;
						}
					}
				}
			}
			#pragma omp critical
			{
				if(min < minShared) minShared = min;
				if(max > maxShared) maxShared = max;
			}
		}
		MPI_Allreduce (&minShared, &GMIN, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		MPI_Allreduce (&maxShared, &GMAX, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		
		ERR = log(GMIN/GMAX);

    
    	#pragma omp parallel for private (i,j, scaled_max) 
    	for(i=0;i<neqns;i++){
    		scaled_max = 0.0;
    		
    		if(areValuesScaled == 0){
    			scaled_max = log(amax[i]);
    		}
    		
    		for(j=ia[i];j<ia[i+1];j++){
    		    if(INFO > 4) printf("before %d %d, %lf \n", i, ja[j], a[j]);
    			if(a[j] != 0.0){
    				a[j] = SCALE * (log(fabs(a[j])) - scaled_max - ERR);
    			}else{
    				a[j] = -SCALE * fabs(log(GMAX));
    			}
    		    if(INFO > 4) printf("after %d %d, %lf \n", i, ja[j], a[j]);
    		}
    	}
	}	
	stoptime = MPI_Wtime();	
	TIME_WEIGHT_SCALING = stoptime - starttime;
	
	return 0;
}
/**
	Parameters
	- *ia: length of neqns + 1
	- *ja: legnth of nnz
	- *a: weight, legnth of nnz
	- *neqns: number of vertices
	- *nnz: number of edges
	- *perm: permutation vector with length of neqns
	- *scaling: for scaling variables with length of 2 * neqns
	- *option: 
	    -option[0]: 0 no output, 1 with info of summary, 2 with number of iterations in each round, 3 for debug, return -1 if it did not find perfect matching
		-option[1]: 
		-option[2]: 
		-option[3]: maximal iterations in auction algorithm, default 1000000000, or use small number of iterations as an approximation algorithm
		-option[4]: maximal iterations of forward auction algorithm, default 100
	return: number of matching
*/
int auction_(int * ia, int * ja, double * a, int * local_neqns, int * local_nnz, int * perm, double * scaling, int * rBlockBndry, double * option, int * neqns){
	double start_time0 = 0.0, stop_time0 = 0.0, start_time2 = 0.0, stop_time2 = 0.0,start_time3 = 0.0, stop_time3 = 0.0;	
	
	start_time0 = MPI_Wtime();
	
	int mype, npes;
	MPI_Comm_size(MPI_COMM_WORLD, &npes);
	MPI_Comm_rank(MPI_COMM_WORLD, &mype);

	MPIType mpiVar;
	mpiVar.npes = npes;
	mpiVar.mype = mype;
	
	MPI_Bcast(neqns, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// every processsor reads the available options	
	INFO = option[0];
	if(option[1] > 0) FACTOR = option[1];
	if(option[2] > 0) THETA = option[2];
	if(option[3] > 0) OMEGA = option[3];
	if(option[4] > 0) MAX_ITER = option[4];
	if(option[5] > 0) ITER_AUCTION = option[5];
	if(option[6] > 0) APPROX_FACTOR = option[6];
	if(option[7] >= 0) isMatrixDistributed = option[7];
	if(option[8] >= 0) isBoundaryInfoDistributed = option[8];
	if(option[9] >= 0) isMatrixDense = option[9];
	if(option[10] >= 0) returnDualVariables = option[10];
	if(option[11] >= 0) useEpsilonScaling = option[11];
	if(option[12] >= 0) areValuesScaled = option[12];
	
	if(INFO > 2){	
		if(mpiVar.mype == 0){
			int i;
			for(i = 0; i < 13; i++){
				printf("Option %d is set to %e\n", i, option[i]);
			}
		}
	}

	ITER_ALL = 0;
	TIME_DISTRIBUTION = 0.0;
	TIME_WEIGHT_SCALING = 0.0;
	TIME_GET_BEST_OBJECT = 0.0;
	TIME_BID = 0.0;
	TIME_UPDATE = 0.0;
	TIME_AUCTION = 0.0;
	TIME_MERGE_MATCHING = 0.0;
	TIME_ADD_MATCHING = 0.0;
	TIME_ALL = 0.0;
	TIME_COMMUNICATION_MERGE = 0.0;
	TIME_COMMUNICATION_UPDATE = 0.0;
	TIME_PARTITION = 0.0;
	TIME_SCALING = 0.0;
		
	MSGS = 0;
	
	AdjMatrixType matrix;
	int nMatch = 0;
	int i = 0, j = 0, tmp = 0; 
	
	int * row = (int *) malloc(*neqns * sizeof(int));	
	int * col = (int *) malloc(*neqns * sizeof(int));
	int * displs = (int *) malloc((mpiVar.npes+1) * sizeof(int));
	
	#pragma omp parallel for private (i) 
	for(i = 0; i < *neqns; i++){
		row[i] = -1;
		col[i] = -1;
	}
	
	//memset(row, -1, *neqns*sizeof(int));
	//memset(col, -1, *neqns*sizeof(int));
		
	if(isBoundaryInfoDistributed == 0){	
		if(mpiVar.mype == 0){
			displs[0] = 0;
	
			int avg_neqns = *neqns / mpiVar.npes;
			for(i = 0; i < mpiVar.npes-1; i++){
				displs[i+1] = displs[i] + avg_neqns;
			}
	    	displs[mpiVar.npes] = *neqns;
		}
	}else{
		if(mpiVar.mype == 0){
			#pragma omp parallel for private (i) 
			for(i = 0; i < mpiVar.npes; i++){
				displs[i] = rBlockBndry[i];
			}
			displs[mpiVar.npes] = *neqns;
		}	
	}

	MPI_Bcast(displs, mpiVar.npes+1, MPI_INT, 0, MPI_COMM_WORLD);
	
	TIME_PARTITION = TIME_PARTITION + (stop_time2 - start_time2);
	
	if(isMatrixDistributed == 0){	
		matrix = distributeVtxs(*neqns, ia, ja, a, displs, mpiVar);
	}else{
		matrix.ia = ia;
		matrix.ja = ja;
		matrix.a = a;
		matrix.neqns = *local_neqns; 
		matrix.nnz = *local_nnz; 
	}

#if 0
	if(mpiVar.mype == 0){
		FILE *file;
		file = fopen("rank0.txt", "a+");
		fprintf(file, "%d\n", matrix.neqns);
		fprintf(file, "%d\n", matrix.nnz);
		for(i=0; i<=matrix.neqns;i++){
			fprintf(file, "%d\n", matrix.ia[i]);
		}
		for(i=0; i<matrix.nnz;i++){
			fprintf(file, "%d\n", matrix.ja[i]);
		}
		for(i=0; i<matrix.nnz;i++){
			fprintf(file, "%f\n", matrix.a[i]);
		}
	
	}else{
		
		FILE *file;
		file = fopen("rank1.txt", "a+");
		fprintf(file, "%d\n", matrix.neqns);
		fprintf(file, "%d\n", matrix.nnz);
		for(i=0; i<=matrix.neqns;i++){
			fprintf(file, "%d\n", matrix.ia[i]);
		}
		for(i=0; i<matrix.nnz;i++){
			fprintf(file, "%d\n", matrix.ja[i]);
		}
		for(i=0; i<matrix.nnz;i++){
			fprintf(file, "%f\n", matrix.a[i]);
		}
	}
#endif
	//printf("rank %d neqns %d, nnz %d\n", mpiVar.mype, matrix.neqns, matrix.nnz);
	
	double * amax = malloc(matrix.neqns * sizeof(double));
	
	SCALE = *neqns+1;
	double scal_val = 1.0 / SCALE;
	END_EPS = scal_val < LIMIT ? scal_val : LIMIT;
	
	scaleWeight(matrix.ia, matrix.ja, matrix.a, amax, matrix.neqns, matrix.nnz);
	
	/////////////////////////////////////////////
	for(i = 0; i < ITER_AUCTION; i++){

		nMatch = AuctionAlgorithmParallel(matrix.ia, matrix.ja, matrix.a, matrix.neqns, displs, row, col, perm, scaling, amax, *neqns, mpiVar);
		if(nMatch >= (*neqns*APPROX_FACTOR) || ITER_ALL >= MAX_ITER) break;
	
		MPI_Bcast(row, *neqns, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(col, *neqns, MPI_INT, 0, MPI_COMM_WORLD);
	}
	/////////////////////////////////////////////
	if(mpiVar.mype == 0){
		if(nMatch < *neqns){
			option[1] = -1;
			printf("MATCHED VERTICES: nMatch %d/%d\n", nMatch, *neqns);	

			#pragma omp parallel for private (i) 
			for(i = 0; i < *neqns; i++) col[i] = -1;
			
			#pragma omp parallel for private (i) 
			for(i = 0; i < *neqns; i++) {
				if(perm[i] > -1) {
					col[perm[i]] = 0;
				}
			}
			nMatch = 0; tmp = 0;
			
			for(i = 0; i < *neqns; i++){	
				if(perm[i] == -1){
					for (j = tmp; j < *neqns ; j++) {	
						if(col[j] == -1){
							perm[i] = j;
							col[j] = i;
							tmp = j + 1;
							break;
						}
					}
					if(perm[i] == -1) printf("ERROR: no col found for row %d\n",i);
				}
				if(perm[i] > -1){
					nMatch++;
				}
			}
			TIME_ADD_MATCHING = TIME_ADD_MATCHING + (stop_time3 - start_time3);
		}
	}

	if(isMatrixDistributed == 0){
		free(matrix.ia);
		free(matrix.ja);
		free(matrix.a);
	}
	
	free(amax);
	free(displs);
	free(col);
	free(row);	
	
	stop_time0 = MPI_Wtime();
	TIME_ALL = stop_time0 - start_time0;

    if(mpiVar.mype == 0){ 
         	
        if(INFO >= 2){
        printf("\n--------------------------\n");
        printf("PAA: npes %d, Matching %d/%d, # iter %d\n", mpiVar.npes, nMatch, *neqns, ITER_ALL);
        printf("--------------------------\n");
        printf("TIME COMPLETE AUCTION: %f\n", TIME_ALL);
        printf("--------------------------\n");
        printf("--TIME PARTITIONING: %f\n", TIME_PARTITION);
        printf("--------------------------\n");
        printf("--TIME DISTRIBUTION: %f\n", TIME_DISTRIBUTION);
        printf("--------------------------\n");
        printf("--TIME ACTUAL AUCTION: %f\n", TIME_ALL - TIME_DISTRIBUTION - TIME_PARTITION);
        printf("----TIME COMPUTATION TOTAL %f\n", TIME_BID + TIME_WEIGHT_SCALING);
        printf("----TIME COMMUNICATION TOTAL %f\n", TIME_COMMUNICATION_MERGE + TIME_COMMUNICATION_UPDATE);
        printf("--------------------------\n");
        printf("-------TIME WEIGHT SCALING: %f\n", TIME_WEIGHT_SCALING);
        printf("--------------------------\n");
        printf("-------TIME PARALLEL AUCTION ALGORITHM: %f\n", TIME_AUCTION);
        printf("----------TIME BIDDING %f (GETTING BEST OBJECT %f )\n", TIME_BID, TIME_GET_BEST_OBJECT);
        printf("----------TIME UPDATING %f (COMMUNICATION %f, # MSG %ld)\n", TIME_UPDATE, TIME_COMMUNICATION_UPDATE, MSGS*2*mpiVar.npes);
        printf("--------------------------\n");
        printf("----------TIME MERGE_MATCHING: %f (COMMUNICATION %f)\n", TIME_MERGE_MATCHING, TIME_COMMUNICATION_MERGE);
        printf("--------------------------\n");
        printf("--TIME ADD_MATCHING: %f\n", TIME_ADD_MATCHING);
        printf("--TIME SCALING: %f\n", TIME_SCALING);
        printf("--------------------------\n");
        printf("--------------------------\n\n\n");
        printf("PAA: npes %d/%d, Matching %d/%d, # iter %d\n", mpiVar.npes, omp_get_max_threads(), nMatch, *neqns, ITER_ALL);
        printf("TIMINGS: AUCTION %f, TOTAL %f (COMMUNICATION %f, COMPUTATION %f, MSG %ld)\n", TIME_ALL - TIME_DISTRIBUTION - TIME_PARTITION, TIME_ALL, TIME_COMMUNICATION_MERGE + TIME_COMMUNICATION_UPDATE, TIME_WEIGHT_SCALING + TIME_BID + TIME_MERGE_MATCHING-TIME_COMMUNICATION_MERGE+TIME_UPDATE-TIME_COMMUNICATION_UPDATE, MSGS*2*mpiVar.npes);
        }else{
            if(INFO == 1){     
                printf("PAA: npes %d/%d, Matching %d/%d, # iter %d\n", mpiVar.npes, omp_get_max_threads(), nMatch, *neqns, ITER_ALL);
                printf("TIMINGS: AUCTION %f, TOTAL %f (COMMUNICATION %f, COMPUTATION %f, MSG %ld)\n", TIME_ALL - TIME_DISTRIBUTION - TIME_PARTITION, TIME_ALL, TIME_COMMUNICATION_MERGE + TIME_COMMUNICATION_UPDATE,  TIME_WEIGHT_SCALING + TIME_BID + TIME_MERGE_MATCHING-TIME_COMMUNICATION_MERGE+TIME_UPDATE-TIME_COMMUNICATION_UPDATE, MSGS*2*mpiVar.npes);
            }
        }
    }
	return nMatch;
}
