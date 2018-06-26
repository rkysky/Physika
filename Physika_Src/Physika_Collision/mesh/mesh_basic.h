#pragma once
#include "Physika_Collision/data_struct/vec3f.h"
#include "Physika_Collision/data_struct/tri.h"
#include "Physika_Framework\Framework\ModuleCollision.h"

namespace Physika {

	class mesh_basic:public CollisionModule{
		unsigned int _num_vtx;
		unsigned int _num_tri;

		tri3f *_tris;

		// used by time integration
		vec3f *_vtxs;
		vec3f *_ivtxs; // initial positions
		vec3f *_ovtxs; // previous positions
		vec3f *_nrms;
		bool _first;

		virtual unsigned int getNbVertices() const = 0;
		virtual unsigned int getNbFaces()const = 0;
		virtual vec3f *getVtxs() const = 0;
		virtual vec3f *getOVtxs() const = 0;
		virtual vec3f *getIVtxs() const = 0;

	};
}