#include <body_force.h>

void body_force(Eigen::VectorXd &q, Eigen::VectorXd &qdot, Eigen::Vector3d &g, double dt){
	Eigen::VectorXd g_gen;
	g_gen.resize(q.size());
    
    for(int i = 0; i < q.size()/3; i ++){
    	g_gen.segment<3>(i*3) = g;
    }

	qdot += dt*g_gen;
}