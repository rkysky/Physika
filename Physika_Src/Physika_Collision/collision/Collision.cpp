#include "Collision.h"
#include "collid.h"

Collision* Collision::instance = NULL;

void Collision::Transform_Pair(unsigned int a,unsigned int b) {
	if(_is_first ==true)
		mesh_pairs.push_back(mesh_pair(a, b));
}

void Collision::Transform_Mesh(unsigned int numVtx, unsigned int numTri, vector<unsigned int> tris, vector<double> vtxs,int m_id,
	bool able_selfcollision) {

	tri3f* _tris=new tri3f[numTri];
	vec3f* _vtxs = new vec3f[numVtx];

	for (int i = 0; i < numVtx; i++) {
		_vtxs[i] = vec3f(vtxs[i * 3], vtxs[i * 3 + 1], vtxs[i * 3 + 2]);
	}

	for (int i = 0; i < numTri; i++) {
		_tris[i] = tri3f(tris[i * 3], tris[i * 3 + 1], tris[i * 3 + 2]);
	}
	
	//mesh* m =new mesh(numVtx, numTri, _tris, _vtxs);
	//dl_mesh.push_back(m);
	if (_is_first == true) {
		mesh* m =new mesh(numVtx, numTri, _tris, _vtxs);
		bodys.push_back(CollisionDate(m, able_selfcollision));
	}
		
	else {
		memcpy(bodys[m_id].ms->_vtxs,_vtxs, numVtx*sizeof(vec3f));
		delete _tris;
		delete _vtxs;
	}
		
}

 void Collision::Collid(){
	contact_pairs.clear();
	body_collide_gpu(mesh_pairs, bodys, contact_pairs);

	_is_first = false;
}