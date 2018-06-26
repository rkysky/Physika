#pragma once
#include "physika_collision/data_struct/def.cuh"
#include "Physika_Framework\Framework\ModuleCollision.h"
// very robust triangle intersection test
// uses no divisions
// works on coplanar triangles


namespace Physika {
	class tricontact:public CollisionModule 
	{
	public:
		tricontact();
		virtual ~tricontact();

		static __device__ int project3(const REAL3 &ax,
			const REAL3 &p1, const REAL3 &p2, const REAL3 &p3);

		static __device__ int project6(REAL3 &ax,
			REAL3 &p1, REAL3 &p2, REAL3 &p3,
			REAL3 &q1, REAL3 &q2, REAL3 &q3);

		static __device__ bool tri_contact(REAL3 &P1, REAL3 &P2, REAL3 &P3,
			REAL3 &Q1, REAL3 &Q2, REAL3 &Q3);

	};
}
