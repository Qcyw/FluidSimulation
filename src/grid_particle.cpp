#include <grid_particle.h>

void grid_particle(Eigen::VectorXd &q, Eigen::VectorXd &qdot, Eigen::VectorXd &u, Eigen::VectorXd &v, Eigen::VectorXd &u_prev, Eigen::VectorXd &v_prev, const Eigen::Vector2i corner, int mode, double alpha, double dx){

	//for flip:
	Eigen::VectorXd du, dv;
	int nx, ny; 
	nx = corner(0); 
	ny = corner(1);
	du.resize(u.size());
	dv.resize(v.size());
	du.setZero();
	dv.setZero();
	du = u-u_prev;
	dv = v-u_prev;

	switch(mode){
		case 0: 
			bilinearInterpolation(q, qdot, du, dv, corner, dx);
			break;
		case 1:
			bilinearInterpolation(q, qdot, u, v, corner, dx);
			break;
		case 2:
			Eigen::VectorXd qdotflip;
			qdotflip.resize(qdot.size());
			bilinearInterpolation(q, qdot, u, v, corner, dx);
			bilinearInterpolation(q, qdotflip, du, dv, corner, dx);
			qdot = (1-alpha)*qdot + alpha*qdotflip;
			break;
	}
}