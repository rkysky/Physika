//atomicAdd(XX, 1) -> atomicInc !!!!!!!
#define OUTPUT_TXT

// CUDA Runtime
#include <cuda_runtime.h>

#include <cuda_profiler_api.h>
#include <assert.h>

// Utilities and system includes
#include <helper_functions.h>  // helper for shared functions common to CUDA SDK samples
#include <helper_cuda.h>       // helper for CUDA error checking

#include "Physika_Collision/data_struct/vec3.cuh"
#include "Physika_Collision/data_struct/tools.cuh"
#include "Physika_Collision/data_struct/box.cuh"
#include "Physika_Collision/data_struct/tri3f.cuh"
#include "Physika_Collision/bvh/bvh.cuh"
#include "Physika_Collision/data_struct/pair.cuh"
#include "Physika_Collision/tri-contact/tri-contact.cuh"

#include <math.h>
#include <stdarg.h>

#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <thrust/sort.h>

#include <string>

using namespace std;

typedef unsigned int uint;


typedef struct {
	uint  numFace, numVert;
	REAL3 *_dx, *_dx0;
	tri3f *_df;
	g_box *_dfBx;

	// init function
	void init()
	{
		numFace = 0;
		numVert = 0;
		_dx0 = _dx = NULL;
		_df = NULL;
		_dfBx = NULL;
	}

	void destroy()
	{
		if (_dx == NULL) return;

		checkCudaErrors(cudaFree(_dx));
		checkCudaErrors(cudaFree(_dx0));
		checkCudaErrors(cudaFree(_df));
		checkCudaErrors(cudaFree(_dfBx));
	}

	void computeWSdata(REAL thickness, bool ccd);
} g_mesh;

//=======================================================

cudaDeviceProp deviceProp;
extern void initPairsGPU();

void initGPU()
{
	int devID = 0;
	cudaGetDevice(&devID);
	checkCudaErrors(cudaGetDeviceProperties(&deviceProp, devID));

	initPairsGPU();
}

//=======================================================

g_mesh theCloth;
g_bvh* theBVH;
g_front theFront;
g_pair thePairs[2]; // potentially colliding pairs
g_pair retPairs; //results

//rky
int BVH_NUM;
//=======================================================
//rky
void init(int num){
	BVH_NUM = num;
	theBVH=new g_bvh[num];
	//cudaMalloc((void **)&theFront, num*sizeof(g_front));
	theFront.init();
}
//rky
void front_clear() {
	theFront.clear();
}

void initPairsGPU()
{
	//pairs[0].init(MAX_PAIR_NUM); // MAX_PAIR_NUM);
	thePairs[1].init(MAX_PAIR_NUM);
	retPairs.init(MAX_PAIR_NUM / 10);
}

void pushMesh2GPU(int  numFace, int numVert, void *faces, void *nodes)
{
	theCloth.init();

	theCloth.numFace = numFace;
	theCloth.numVert = numVert;

	cudaMalloc((void **)&theCloth._df, numFace*sizeof(tri3f));
	cudaMalloc((void **)&theCloth._dfBx, numFace*sizeof(g_box));
	cudaMalloc((void **)&theCloth._dx, numVert*sizeof(REAL3));
	cudaMalloc((void **)&theCloth._dx0, numVert*sizeof(REAL3));

	cudaMemcpy(theCloth._df, faces, sizeof(tri3f)*numFace, cudaMemcpyHostToDevice);
	cudaMemcpy(theCloth._dx, nodes, sizeof(REAL3)*numVert, cudaMemcpyHostToDevice);
	cudaMemcpy(theCloth._dx0, theCloth._dx, sizeof(REAL3)*numVert, cudaMemcpyDeviceToDevice);

	theCloth.computeWSdata(0, false);
}

void updateMesh2GPU(void *nodes)
{
	cudaMemcpy(theCloth._dx0, theCloth._dx, sizeof(REAL3)*theCloth.numVert, cudaMemcpyDeviceToDevice);
	cudaMemcpy(theCloth._dx, nodes, sizeof(REAL3)*theCloth.numVert, cudaMemcpyHostToDevice);
	theCloth.computeWSdata(0, false);
}

//=======================================================

void pushBVHIdx(int max_level, unsigned int *level_idx, int i)
{
	theBVH[i]._max_level = max_level;
	theBVH[i]._level_idx = new uint[max_level];
	memcpy(theBVH[i]._level_idx, level_idx, sizeof(uint)*max_level);
}

void pushBVH(unsigned int length, int *ids, int i)
{
	theBVH[i]._num = length;
	checkCudaErrors(cudaMalloc((void**)&theBVH[i]._bvh, length*sizeof(int) * 2));
	checkCudaErrors(cudaMemcpy(theBVH[i]._bvh, ids, length*sizeof(int) * 2, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMalloc((void**)&theBVH[i]._bxs, length*sizeof(g_box)));
	checkCudaErrors(cudaMemset(theBVH[i]._bxs, 0, length*sizeof(g_box)));
	theBVH[i].hBxs = NULL;

	//rky
	theBVH[i]._triBxs = theCloth._dfBx;
	//theBVH[i]._triBxs =NULL;
	theBVH[i]._triCones = NULL;
}

void pushBVHLeaf(unsigned int length, int *idf, int i)
{
	checkCudaErrors(cudaMalloc((void**)&theBVH[i]._bvh_leaf, length*sizeof(int)));
	checkCudaErrors(cudaMemcpy(theBVH[i]._bvh_leaf, idf, length*sizeof(int), cudaMemcpyHostToDevice));
}

//======================================================


void refitBVH_Serial(int bvh_id, int length)
{

	refit_serial_kernel << <1, 1, 0 >> >
		(theBVH[bvh_id]._bvh, theBVH[bvh_id]._bxs, theBVH[bvh_id]._triBxs,
		theBVH[bvh_id]._cones, theBVH[bvh_id]._triCones,
		length == 0 ? theBVH[bvh_id]._num : length);

	getLastCudaError("refit_serial_kernel");
	cudaThreadSynchronize();
}

void refitBVH_Parallel(int bvh_id, int st, int length)
{
	BLK_PAR(length);

	refit_kernel << < B, T >> >
		(theBVH[bvh_id]._bvh, theBVH[bvh_id]._bxs, theBVH[bvh_id]._triBxs,
		theBVH[bvh_id]._cones, theBVH[bvh_id]._triCones,
		st, length);

	getLastCudaError("refit_kernel");
	cudaThreadSynchronize();
}

//rky
void refitBVH(int bvh_id)
{
	// before refit, need to get _tri_boxes !!!!
	// copying !!!
	for (int i = theBVH[bvh_id]._max_level - 1; i >= 0; i--) {
		int st = theBVH[bvh_id]._level_idx[i];
		int ed = (i != theBVH[bvh_id]._max_level - 1) ?
			theBVH[bvh_id]._level_idx[i + 1] - 1 : theBVH[bvh_id]._num - 1;

		int length = ed - st + 1;
		if (i < 5) {
			refitBVH_Serial(bvh_id, length + st);
			break;
		}
		else
		{
			refitBVH_Parallel(bvh_id, st, length);
		}
	}
}


void refitBVH()
{
	// before refit, need to get _tri_boxes !!!!
	// copying !!!
	for (int k = 0; k < BVH_NUM; k++){
		int bvh_id = k;
		for (int i = theBVH[bvh_id]._max_level - 1; i >= 0; i--) {
			int st = theBVH[bvh_id]._level_idx[i];
			int ed = (i != theBVH[bvh_id]._max_level - 1) ?
				theBVH[bvh_id]._level_idx[i + 1] - 1 : theBVH[bvh_id]._num - 1;

			int length = ed - st + 1;
			if (i < 5) {
				refitBVH_Serial(bvh_id, length + st);
				break;
			}
			else
			{
				refitBVH_Parallel(bvh_id, st, length);
			}
		}
	}
}

//===============================================

void pushFront(int num, unsigned int *data,unsigned int *data_id)
{
	g_front *f = &theFront;

	//rky
	//f->init();
	f->push(num, (uint4 *)data,data_id);
}

//===============================================
// show memory usage of GPU
void  reportMemory(char *tag)
{
	//return;

#ifdef OUTPUT_TXT
	size_t free_byte;
	size_t total_byte;
	cudaError_t cuda_status = cudaMemGetInfo(&free_byte, &total_byte);

	if (cudaSuccess != cuda_status) {
		printf("Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(cuda_status));
		exit(1);
	}

	REAL free_db = (REAL)free_byte;
	REAL total_db = (REAL)total_byte;
	REAL used_db = total_db - free_db;
	printf("%s: GPU memory usage: used = %f, free = %f MB, total = %f MB\n",
		tag, used_db / 1024.0 / 1024.0, free_db / 1024.0 / 1024.0, total_db / 1024.0 / 1024.0);
#endif
}

//===============================================

#define STACK_SIZE 50
#define EMPTY (nIdx == 0)

#define PUSH_PAIR(nd1, nd2)  {\
	nStack[nIdx].x = nd1;\
	nStack[nIdx].y = nd2;\
	nIdx++;\
}

#define POP_PAIR(nd1, nd2) {\
	nIdx--;\
	nd1 = nStack[nIdx].x;\
	nd2 = nStack[nIdx].y;\
}

#define NEXT(n1, n2) 	POP_PAIR(n1, n2)


inline __device__ void pushToFront(int a, int b, uint4 *front, uint *idx, uint ptr)
{
	//	(*idx)++;
	if (*idx < MAX_FRONT_NUM)
	{
		uint offset = atomicAdd(idx, 1);
		front[offset] = make_uint4(a, b, 0, ptr);
	}
}

inline __device__ void sproutingAdaptive(int left, int right,
	int *bvhA, g_box *bxsA, int *bvhB, g_box *bxsB,
	uint4 *front, uint *frontIdx,
	uint2 *pairs, uint *pairIdx, bool update, uint ptr)
{
	uint2 nStack[STACK_SIZE];
	uint nIdx = 0;

	for (int i = 0; i<4; i++)
	{
		if (isLeaf(left, bvhA) && isLeaf(right, bvhB)) {
			pushToFront(left, right, front, frontIdx, ptr);
		}
		else {
			if (!overlaps(left, right, bxsA, bxsB)) {
				pushToFront(left, right, front, frontIdx, ptr);
			}
			else {
				if (isLeaf(left, bvhA)) {
					PUSH_PAIR(left, getLeftChild(right, bvhB));
					PUSH_PAIR(left, getRightChild(right, bvhB));
				}
				else {
					PUSH_PAIR(getLeftChild(left, bvhA), right);
					PUSH_PAIR(getRightChild(left, bvhA), right);
				}
			}
		}

		if (EMPTY)
			return;

		NEXT(left, right);
	}

	while (!EMPTY) {
		NEXT(left, right);
		pushToFront(left, right, front, frontIdx, ptr);
	}
}

inline __device__ void sprouting(int left, int right,
	int *bvhA, g_box *bxsA, int *bvhB, g_box *bxsB,
	uint4 *front, uint *frontIdx,
	int2 *pairs, uint *pairIdx, bool update, uint ptr)
{
	uint2 nStack[STACK_SIZE];
	uint nIdx = 0;

	while (1)
	{
		if (isLeaf(left, bvhA) && isLeaf(right, bvhB)) {
			if (update)
				pushToFront(left, right, front, frontIdx, ptr);

			if (overlaps(left, right, bxsA, bxsB))
				addPair(getTriID(left, bvhA), getTriID(right, bvhB), pairs, pairIdx);
		}
		else {
			if (!overlaps(left, right, bxsA, bxsB)) {
				if (update)
					pushToFront(left, right, front, frontIdx, ptr);

			}
			else {
				if (isLeaf(left, bvhA)) {
					PUSH_PAIR(left, getLeftChild(right, bvhB));
					PUSH_PAIR(left, getRightChild(right, bvhB));
				}
				else {
					PUSH_PAIR(getLeftChild(left, bvhA), right);
					PUSH_PAIR(getRightChild(left, bvhA), right);
				}
			}
		}

		if (EMPTY)
			return;

		NEXT(left, right);
	}
}

__device__ void doPropogate(
	uint4 *front,uint *bvh_id,g_bvh *bvh, uint *frontIdx, int num,
	int2 *pairs, uint *pairIdx, bool update, tri3f *Atris, int idx, bool *flags)
{
	uint4 node = front[idx];

	uint _bvhid[2];
	_bvhid[0] = bvh_id[idx*2];
	_bvhid[1] = bvh_id[idx*2 + 1];

	if (node.z != 0) {
#if defined(_DEBUG) || defined(OUTPUT_TXT)
		atomicAdd(frontIdx + 1, 1);
#endif
		return;
	}

#ifdef USE_NC
	if (flags != NULL && flags[node.w] == 0) {
#if defined(_DEBUG) || defined(OUTPUT_TXT)
		atomicAdd(frontIdx + 2, 1);
#endif
		return;
	}
#endif


	uint left = node.x;
	uint right = node.y;
	int *bvhA = bvh[_bvhid[0]]._bvh;
	g_box *bxsA = bvh[_bvhid[0]]._bxs;
	int *bvhB = bvh[_bvhid[1]]._bvh;
	g_box *bxsB = bvh[_bvhid[1]]._bxs;

	if (isLeaf(left, bvhA) && isLeaf(right, bvhB)) {
		if (overlaps(left, right, bxsA, bxsB))
			if (_bvhid[0] != _bvhid[1])
				addPair(getTriID(left, bvhA), getTriID(right, bvhB), pairs, pairIdx);
			else { // for self ccd, we need to remove adjacent triangles, they will be processed seperatedly with orphan set
				if (!covertex(getTriID(left, bvhA), getTriID(right, bvhB), Atris))
					addPair(getTriID(left, bvhA), getTriID(right, bvhB), pairs, pairIdx);
			}
			return;
	}

	if (!overlaps(left, right, bxsA, bxsB))
		return;

	if (update)
		front[idx].z = 1;

	int ptr = node.w;
	if (isLeaf(left, bvhA)) {
		sprouting(left, getLeftChild(right, bvhB), bvhA, bxsA, bvhB, bxsB, front, frontIdx, pairs, pairIdx, update, ptr);
		sprouting(left, getRightChild(right, bvhB), bvhA, bxsA, bvhB, bxsB, front, frontIdx, pairs, pairIdx, update, ptr);
	}
	else {
		sprouting(getLeftChild(left, bvhA), right, bvhA, bxsA, bvhB, bxsB, front, frontIdx, pairs, pairIdx, update, ptr);
		sprouting(getRightChild(left, bvhA), right, bvhA, bxsA, bvhB, bxsB, front, frontIdx, pairs, pairIdx, update, ptr);
	}
}


__global__ void kernelPropogate(uint4 *front, uint *bvh_id,g_bvh *bvh, uint *frontIdx, int num,
	int2 *pairs, uint *pairIdx, bool update, tri3f *Atris, int stride, bool *flags)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	for (int i = 0; i<stride; i++) {
		int j = idx*stride + i;
		if (j >= num)
			return;

		doPropogate(front,bvh_id, bvh,frontIdx, num,
		 pairs, pairIdx, update, Atris, j, flags);
	}
}

int g_front::propogate(bool &update, bool ccd)
{
	uint dummy[1];
	cutilSafeCall(cudaMemcpy(dummy, _dIdx, 1 * sizeof(uint), cudaMemcpyDeviceToHost));
#ifdef OUTPUT_TXT
	//printf("Before propogate, length = %d\n", dummy[0]);
#endif

#if defined(_DEBUG) || defined(OUTPUT_TXT)
	uint dummy2[5] = { 0, 0, 0, 0, 0 };
	cutilSafeCall(cudaMemcpy(_dIdx + 1, dummy2, 5 * sizeof(int), cudaMemcpyHostToDevice));
#endif

	if (dummy[0] != 0) {
		//rky
		//g_bvh *pb1 = &theBVH[id1];
		//g_bvh *pb2 = &theBVH[id2];
		tri3f *faces = theCloth._df;

		int stride = 4;
#ifdef FIX_BT_NUM
		BLK_PAR2(dummy[0], stride);
#else
		BLK_PAR3(dummy[0], stride, getBlkSize((void *)kernelPropogate));
#endif
		g_bvh* _thebvh;
		cutilSafeCall(cudaMalloc((void**)&_thebvh, BVH_NUM*sizeof(g_bvh)));
		cutilSafeCall(cudaMemcpy(_thebvh, theBVH, BVH_NUM*sizeof(g_bvh),cudaMemcpyHostToDevice));
			
		//rky
		kernelPropogate << < B, T >> >
			(_dFront,bvh_id,_thebvh, _dIdx, dummy[0],
			thePairs[1]._dPairs, thePairs[1]._dIdx, update, faces, stride,  theBVH[0]._ctFlags );
		//thePairs[self]._dPairs, thePairs[self]._dIdx, update, faces, stride, (self && !ccd) ? theBVH[1]._ctFlags : NULL);

		cudaThreadSynchronize();
		getLastCudaError("kernelPropogate");
	}

	cutilSafeCall(cudaMemcpy(dummy, _dIdx, 1 * sizeof(uint), cudaMemcpyDeviceToHost));
#ifdef OUTPUT_TXT
	//printf("After propogate, length = %d\n", dummy[0]);
#endif

#if defined(_DEBUG) || defined(OUTPUT_TXT)
	cutilSafeCall(cudaMemcpy(dummy2, _dIdx + 1, 5 * sizeof(int), cudaMemcpyDeviceToHost));
	//printf("Invalid = %d, NC culled = %d\n", dummy2[0], dummy2[1]);
#endif

	if (update && dummy[0] > SAFE_FRONT_NUM) {
		printf("Too long front, stop updating ...\n");
		update = false;
	}

	if (dummy[0] > MAX_FRONT_NUM) {
		printf("Too long front, exiting ...\n");
		exit(0);
	}
	return dummy[0];
}

//===============================================

__global__ void
kernel_face_ws(tri3f *face, REAL3 *x, REAL3 *ox, g_box *bxs, bool ccd, REAL thickness, int num)
{
	LEN_CHK(num);

	int id0 = face[idx].id0();
	int id1 = face[idx].id1();
	int id2 = face[idx].id2();

	REAL3 ox0 = ox[id0];
	REAL3 ox1 = ox[id1];
	REAL3 ox2 = ox[id2];
	REAL3 x0 = x[id0];
	REAL3 x1 = x[id1];
	REAL3 x2 = x[id2];

	bxs[idx].set(ox0, ox1);
	bxs[idx].add(ox2);

	if (ccd) {
		bxs[idx].add(x0);
		bxs[idx].add(x1);
		bxs[idx].add(x2);
	}
	//else
	bxs[idx].enlarge(thickness);
}

void g_mesh::computeWSdata(REAL thickness, bool ccd)
{
	if (numFace == 0)
		return;

	{
		int num = numFace;
		BLK_PAR(num);
		kernel_face_ws << <B, T >> > (
			_df, _dx, _dx, _dfBx, ccd, thickness, num);
		getLastCudaError("kernel_face_ws");
	}
}

//===============================================

__global__ void kernelGetCollisions(
	int2 *pairs, int num, 
	REAL3 *cx, tri3f *ctris, int2 *pairRets, uint *pairIdx,
	int stride)
{
	int idxx = blockDim.x * blockIdx.x + threadIdx.x;

	for (int i = 0; i<stride; i++) {

		int j = idxx*stride + i;
		if (j >= num)
			return;

		int idx = j;

		int2 pair = pairs[idx];
		int fid1 = pair.x;
		int fid2 = pair.y;

		tri3f t1 = ctris[fid1];
		tri3f t2 = ctris[fid2];

#ifdef FOR_DEBUG
		bool find = false;
		if (fid1 == 369 && fid2 == 3564)
			find = true;
		if (fid2 == 369 && fid1 == 3564)
			find = true;
#endif

		REAL3 p0 = cx[t1.id0()];
		REAL3 p1 = cx[t1.id1()];
		REAL3 p2 = cx[t1.id2()];
		REAL3 q0 = cx[t2.id0()];
		REAL3 q1 = cx[t2.id1()];
		REAL3 q2 = cx[t2.id2()];

#ifdef FOR_DEBUG
		if (find) {
			printf("%d: %lf, %lf, %lf\n", t1.id0(), p0.x, p0.y, p0.z);
			printf("%d: %lf, %lf, %lf\n", t1.id1(), p1.x, p1.y, p1.z);
			printf("%d: %lf, %lf, %lf\n", t1.id2(), p2.x, p2.y, p2.z);
			printf("%d: %lf, %lf, %lf\n", t2.id0(), q0.x, q0.y, q0.z);
			printf("%d: %lf, %lf, %lf\n", t2.id1(), q1.x, q1.y, q1.z);
			printf("%d: %lf, %lf, %lf\n", t2.id2(), q2.x, q2.y, q2.z);
		}
#endif
		

		if (tricontact::tri_contact(p0, p1, p2, q0, q1, q2))
			//if (fid1 > fid2)
				//addPair(fid2, fid1, pairRets, pairIdx);
			//else
				addPair(fid1, fid2, pairRets, pairIdx);
	}
}

//===============================================

int g_pair::getCollisions(bool self, g_pair &ret)
{
	int num = length();

#ifdef OUTPUT_TXT
	if (self)
		//printf("self pair = %d\n", num);
	//else
		//printf("inter-obj pair = %d\n", num);
#endif

	if (num == 0)
		return 0;

	ret.clear();
	
	int stride = 4;
#ifdef FIX_BT_NUM
	BLK_PAR3(num, stride, 32);
#else
	BLK_PAR3(num, stride, getBlkSize((void *)kernelGetCollisions));
#endif


	kernelGetCollisions << < B, T >> >(_dPairs, num,
		theCloth._dx, theCloth._df, ret._dPairs, ret._dIdx,stride);

	getLastCudaError("kernelGetCollisions");

	int len = ret.length();
#ifdef OUTPUT_TXT
	//printf("collision num = %d\n", len);
#endif

	return len;
}
//===============================================

int getCollisionsGPU(int *rets)
{
	bool update = false;
	int len = 0;

	TIMING_BEGIN
	thePairs[1].clear();

	refitBVH();

	theFront.propogate(update,false);
	cudaThreadSynchronize();
	

	len = thePairs[1].getCollisions(true, retPairs);
	cudaThreadSynchronize();

	TIMING_END("$$$get_collisions_gpu")

	if (len > 0) {
		cutilSafeCall(cudaMemcpy(rets, retPairs._dPairs, sizeof(uint)*2*len, cudaMemcpyDeviceToHost));
	}

	return len;
}