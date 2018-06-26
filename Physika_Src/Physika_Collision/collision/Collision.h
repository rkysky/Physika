#ifndef SELFCOLLISION
#define SELFCOLLISION

#include "Physika_Collision/mesh/cmesh.h"
#include "Physika_Collision/bvh/tmbvh.hpp"
#include "Physika_Collision/collision/CollisionDate.h"
#include "Physika_Framework/Framework/ModuleCollision.h"

#include "iostream"
#include<vector>

namespace Physika
{

	typedef pair<int, int> mesh_pair;

	class Collision :public CollisionModule{
	public:

		virtual ~Collision() {}

		//Invoke internal function for collision detect.
		void Collid();

		//碰撞体对
		void Transform_Pair(unsigned int a, unsigned int b);

		//碰撞体网格
		//m_id：碰撞体的序号
		void Transform_Mesh(unsigned int numVtx, unsigned int numTri, vector<unsigned int> tris, vector<double> vtxs, int m_id,
			bool able_selfcollision = false
		);

		vector<vector<tri_pair>> getContactPairs() { return contact_pairs; }

		//int Find_Mesh(CollisionDate m);

		int getNumContacts() { return contact_pairs.size(); }

		static Collision* getInstance() {
			if (instance == NULL) {
				instance = new Collision();
				return instance;
			}
			else
				return instance;
		}

		static Collision* instance;

	private:
		vector<CollisionDate> bodys;
		vector<mesh_pair> mesh_pairs;
		vector<vector<tri_pair>> contact_pairs;
		//vector<mesh*> dl_mesh;//delete mesh points

		Collision() :_is_first(true) {}

		bool _is_first;//是否第一次传入网格

	};

}
#endif