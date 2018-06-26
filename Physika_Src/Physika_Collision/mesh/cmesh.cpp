#if defined(WIN32)
#define WIN32_LEAN_AND_MEAN
#  include <windows.h>
#endif

#include <stdio.h>
#include <string.h>
#include "cmesh.h"
#include "physika_collision/data_struct/box.h"
#include "physika_collision/data_struct/vec3f.h"

#include <set>
using namespace std;

// for fopen
#pragma warning(disable: 4996)

inline vec3f update(vec3f &v1, vec3f &v2, vec3f &v3)
{
	vec3f s = (v2-v1);
	return s.cross(v3-v1);
}

inline vec3f
update(tri3f &tri, vec3f *vtxs)
{
	vec3f &v1 = vtxs[tri.id0()];
	vec3f &v2 = vtxs[tri.id1()];
	vec3f &v3 = vtxs[tri.id2()];

	return update(v1, v2, v3);
}

void mesh::update(matrix3f &trf, vec3f &shift)
{
	if (_first) {
		_first = false;
		memcpy(_ivtxs, _vtxs, sizeof(vec3f)*_num_vtx);
	}

	for (unsigned int i=0; i<_num_vtx; i++)
		_vtxs[i] = _ivtxs[i]*trf + shift;

}

void mesh::updateNrms()
{
	for (unsigned int i=0; i<_num_vtx; i++)
		_nrms[i] = vec3f::zero();

	for (unsigned int i=0; i<_num_tri; i++) {
		vec3f n = ::update(_tris[i], _vtxs);
		n.normalize();

		_nrms[_tris[i].id0()] += n;
		_nrms[_tris[i].id1()] += n;
		_nrms[_tris[i].id2()] += n;
	}

	for (unsigned int i=0; i<_num_vtx; i++)
		_nrms[i].normalize();
}


mesh::mesh(unsigned int numVtx, unsigned int numTri, tri3f *tris, vec3f *vtxs)
{
	_first = true;

	_num_vtx = numVtx;
	_num_tri = numTri;

	_tris = tris;
	_vtxs = vtxs;
	_ivtxs = new vec3f[numVtx];
	_ovtxs = new vec3f[numVtx];

	_nrms = new vec3f[numVtx];
}

mesh::~mesh()
{
	delete [] _tris;
	delete [] _vtxs;
	delete [] _ivtxs;
	delete [] _ovtxs;
	delete [] _nrms;
}









