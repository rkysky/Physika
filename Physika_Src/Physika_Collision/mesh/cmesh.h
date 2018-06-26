#pragma once

#include "Physika_Collision/data_struct/vec3f.h"
#include "mesh_basic.h"
#include "Physika_Collision/data_struct/tri.h"
#include "Physika_Collision/data_struct/box.h"
#include "Physika_Collision/data_struct/mat3f.h"
#include <set>
using namespace std;

class mesh:public mesh_basic{
public:
	unsigned int _num_vtx;
	unsigned int _num_tri;
	
	tri3f *_tris;

	// used by time integration
	vec3f *_vtxs;
	vec3f *_ivtxs; // initial positions
	vec3f *_ovtxs; // previous positions
	vec3f *_nrms;
	bool _first;

	mesh(unsigned int numVtx, unsigned int numTri, tri3f *tris, vec3f *vtxs);
	virtual ~mesh();

	virtual unsigned int getNbVertices() const { return _num_vtx; }
	virtual unsigned int getNbFaces() const { return _num_tri; }
	virtual vec3f *getVtxs() const { return _vtxs; }
	virtual vec3f *getOVtxs() const { return _ovtxs;}
	virtual vec3f *getIVtxs() const {return _ivtxs;}

	void update(matrix3f &, vec3f &);

	// calc norms, and prepare for display ...
	void updateNrms();

	BOX bound() {
		BOX a;

		for (int i=0; i<_num_vtx; i++)
			a += _vtxs[i];

		return a;
	}
};
