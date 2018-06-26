#ifndef COLLID
#define COLLID

#include "physika_collision/collision/CollisionDate.h"
#include "physika_collision/collision/Collision.h"

//void body_collide(vector<mesh_pair> mpair, vector<CollisionDate> bodys, vector<vector<tri_pair>> &contacts);

void body_collide_gpu(vector<mesh_pair> mpair, vector<CollisionDate> bodys, vector<vector<tri_pair>> &contacts);

#endif