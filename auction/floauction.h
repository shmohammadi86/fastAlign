/**************************************************************
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "Madan Sathe, Ye Zhao, University of Basel, Switzerland" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

For any kind of questions or comments please send an email to Madan Sathe: madan.sathe@unibas.ch

****************************************************************/


#ifndef FLOAUCTION_H
#define FLOAUCTION_H
struct pairdef {
	int vtx1;
	int vtx2;
	double wPair;
};
typedef struct pairdef PairType;

struct objectdef {
	int idx_object;
	double cost;
	double bid;
};
typedef struct objectdef ObjectType;

struct adjMatrixdef {
	int * ia;
	int * ja;
	double * a;
	double * amin;
	int neqns;
	int nnz;
};
typedef struct adjMatrixdef AdjMatrixType;

struct mpidef {
	int npes;
	int mype;
	MPI_Comm comm;
};
typedef struct mpidef MPIType;


/**
	Parameters
	- *ia: length of local_neqns + 1
	- *ja: legnth of local_nnz
	- *a: weight, legnth of local_nnz
	- *local_neqns: number of local neqns for each proc different
	- *nnz: number of local number of nonzeros
	- *perm: permutation vector of size "global" neqns
	- *scaling: optional: scaling vector, the dual variables, length of 2 * "global" neqns
	- *rBlockBndry: row indizes for distribution the matrix, the start index for each row on the procs
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

	
	return: number of matching
*/
int auction_(int * ia, int * ja, double * a, int * local_neqns, int * local_nnz, int * perm, double * scaling, int * rBlockBndry, double * option, int * neqns);

#endif /*AUCTION_H*/
