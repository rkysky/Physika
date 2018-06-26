#include "Physika_Collision/mesh/cmesh.h"
#include <set>
#include <iostream>
#include <stdio.h>

using namespace std;
#include "Physika_Collision/data_struct/mat3f.h"
#include "Physika_Collision/data_struct/box.h"
#include "Physika_Collision/bvh/tmbvh.hpp"

//rky
#include "Physika_Collision/collision/collid.h"
//rky

bool findd;

#include <omp.h>
//rky
double tmp_timing_start;
# define	TIMING_BEGIN \
	tmp_timing_start = omp_get_wtime();\

# define	TIMING_END(message) \
	{double tmp_timing_finish = omp_get_wtime();\
	double  tmp_timing_duration = tmp_timing_finish - tmp_timing_start;\
	printf("%s: %2.5f seconds\n", (message), tmp_timing_duration);}

//#define POVRAY_EXPORT
#define OBJ_DIR "e:\\temp\\output-objs\\"

//#define VEC_CLOTH

#pragma warning(disable: 4996)

extern void mesh_id(int id, vector<mesh *> &m, int &mid, int &fid);

extern void pushMesh2GPU(int  numFace, int numVert, void *faces, void *nodes);
extern void updateMesh2GPU(void *nodes);
static tri3f *s_faces;
static vec3f *s_nodes;
static int s_numFace = 0, s_numVert = 0;

void updateMesh2GPU(vector <mesh *> &ms)
{
	vec3f *curVert = s_nodes;
	for (int i = 0; i < ms.size(); i++) {
		mesh *m = ms[i];
		memcpy(curVert, m->_vtxs, sizeof(vec3f)*m->_num_vtx);
		curVert += m->_num_vtx;
	}

	updateMesh2GPU(s_nodes);
}

void pushMesh2GPU(vector<mesh *> &ms)
{
	for (int i = 0; i < ms.size(); i++) {
		s_numFace += ms[i]->_num_tri;
		s_numVert += ms[i]->_num_vtx;
	}

	s_faces = new tri3f[s_numFace];
	s_nodes = new vec3f[s_numVert];

	int curFace = 0;
	int vertCount = 0;
	vec3f *curVert = s_nodes;
	for (int i = 0; i < ms.size(); i++) {
		mesh *m = ms[i];
		for (int j = 0; j < m->_num_tri; j++) {
			tri3f &t = m->_tris[j];
			s_faces[curFace++] = tri3f(t.id0() + vertCount, t.id1() + vertCount, t.id2() + vertCount);
		}
		vertCount += m->_num_vtx;

		memcpy(curVert, m->_vtxs, sizeof(vec3f)*m->_num_vtx);
		curVert += m->_num_vtx;
	}

	pushMesh2GPU(s_numFace, s_numVert, s_faces, s_nodes);
}

extern int getCollisionsGPU(int *);
extern void initGPU();

//rky
extern void init(int num);
extern void front_clear();

bool cmp(vector<tri_pair> a, vector<tri_pair> b){
	unsigned int ta[4],tb[4];
	for (int i = 0; i < 2; i++){
		a[i].get(ta[i * 2],ta[i * 2 + 1]);
		b[i].get(tb[i * 2], tb[i * 2 + 1]);
	}
	if (ta[0] != tb[0])
		return ta[0] < tb[0];
	else if (ta[2] != tb[2])
		return ta[2] < tb[2];
	else if (ta[1] != tb[1])
		return ta[1] < tb[1];
	else
		return ta[3] < tb[3];
}

//static vector<bvh*> bvhC;
//static std::vector<mesh *> meshes;


void body_collide_gpu(vector<mesh_pair> mpair, vector<CollisionDate> bodys, vector<vector<tri_pair>> &contacts){
	static vector<bvh*> bvhC;
	//front_list fIntra;
	static std::vector<mesh *> meshes;

	vector<tri_pair> fret;
	static vector<int> _tri_offset;
	if (bvhC.size() == 0) {
		for (int i = 0; i < bodys.size(); i++) {
			meshes.push_back(bodys[i].ms);
			_tri_offset.push_back(i == 0 ? bodys[0].ms->_num_tri : (_tri_offset[i - 1] + bodys[i].ms->_num_tri));
			bvh* tem = new bvh(meshes, i);
			bvhC.push_back(tem);

		}

		initGPU();
		//rky
		init(bvhC.size());
		pushMesh2GPU(meshes);
		for (int i = 0; i < bvhC.size(); i++) {
			bvhC[i]->push2GPU(i, i == 0 ? 0 : _tri_offset[i - 1]);
		}

		int mesh_1, mesh_2;
		for (int i = 0; i < mpair.size(); i++)
		{
			mesh_1 = mpair[i].first;
			mesh_2 = mpair[i].second;
			front_list tem_fl;

			if (mesh_2 != -1)
			{
				//fret.clear();
				bvhC[mesh_1]->collide(bvhC[mesh_2], tem_fl);
				tem_fl.push2GPU(bvhC[mesh_1]->root(), bvhC[mesh_2]->root(), mesh_1, mesh_2);
			}

		}
	}//if

	updateMesh2GPU(meshes);

#ifdef FOR_DEBUG
	vec3f *pts = meshes[0]->getVtxs() + 3126;
	printf("XXXXXXXXXXXX3126: %lf, %lf, %lf\n", pts->x, pts->y, pts->z);
#endif

	int *buffer = new int[10240 * 10];
	int count = getCollisionsGPU(buffer);
	printf("Collision num %d\n",count);

	tri_pair *pairs = (tri_pair *)buffer;
	vector<tri_pair> ret(pairs, pairs + count);

	//Find mesh id and face id
	for (int i = 0; i < count; i++){
		vector<tri_pair> tem;
		int mid1, mid2;
		unsigned int fid1, fid2;
		ret[i].get(fid1, fid2);

		for (int j = 0; j < _tri_offset.size(); j++){
			if (fid1 <= _tri_offset[j]){
				mid1 = j == 0 ? 0 : j;
				break;
			}			
		}
		tem.push_back(tri_pair(mid1, fid1 -1- (mid1 == 0 ? 0 : _tri_offset[mid1 - 1])));

		for (int j = 0; j < _tri_offset.size(); j++){
			if (fid2 <= _tri_offset[j]){
				mid2 = j == 0 ? 0 : j;
				break;
			}
		}

		tem.push_back(tri_pair(mid2, fid2 -1- (mid2 == 0 ? 0 : _tri_offset[mid2 - 1])));

		contacts.push_back(tem);
	}

	std::sort(contacts.begin(), contacts.end(), cmp);
	
	delete[] buffer;	
}
