# 使用基于前线的BVH碰撞检测库

*编译环境
win7+VS2015+CUDA8.0.

*目录
库代码在CollisionDetect\目录下.

*demo
集成了一个基于PBD中RigidBodyCollisionDemo的demo，整个工程文件在build文件夹下.

*使用
-单例类Collision的成员函数
	-Transform_Pair()：
		-导入需要碰撞检测的一对物体的序号mesh_id；
	-Transform_Mesh():
		-导入物体的网格信息；
	-Collid():
		-调用内部函数，完成碰撞检测；
	-getContactPairs():
		-返回检测到碰撞的碰撞对的序号mesh_id和面片的序号face_id；
		-例子((mesh_id1,face_id1),(mesh_id2,face_id2));
	-getNumContacts()：
		-返回发生碰撞的对数；

-调用顺序
	-Transform_Pair(unsigned int a,unsigned int b)
	参数a，b是一对需要碰撞检测的物体的序号，所以在这之前需要按一定顺序给物体编号。这个函数会在单例类内保存所有的碰撞对，并且只在第一次导入的时候更改，之后默认不会有新的碰撞对。
	
	-Transform_Mesh(unsigned int numVtx, unsigned int numTri，
	vector<unsigned int> tris, vector<double> vtxs,int m_id,bool able_selfcollision=false)
	这个函数是将物体网格用单例类内的mesh数据结构进行保存，numVtx和numTri分别是点和面的数目，tris是一个保存了所有三角面片的顶点序号的数组，vtxs是一个保存所有顶点的数组，m_id是对应物体的序号。

	-Collid()
	这个函数会调用内部的碰断检测代码进行碰撞检测，并保存结果。

-e.g
	//导入网格信息和碰撞对
	-for(int i=0;i<nums of object;i++){
		for(int j=0;j<nums of object;j++){
			if(i!=j){
			Transform_Pair(i,j);
			}
		}

		Transform_Mesh(numVtx, numTri，tris, vtxs,m_id);
	}

	Colid();
	
	//如果需要，返回检测到的碰撞对，返回值是一个二维数组，
	保存形式如下：
	(mesh_id1,face_id1,mesh_id2,face_id2)
	(mesh_id3,face_id3,mesh_id4,face_id4)
	：
	：
	：
	
	getContactPairs();
	
	//如果需要，返回检测到的碰撞对的数目
	getNumContacts();




