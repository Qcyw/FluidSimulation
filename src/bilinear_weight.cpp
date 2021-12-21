#include <bilinear_weight.h>

void bilinear_weight(Eigen::VectorXd &q, Eigen::VectorXd &qdot, Eigen::Vector4d &vertex, Eigen::Vector4d &weight, const Eigen::Vector2i corner, int qi){

	// assumes vertex are stored like bottom left, bottom right, top left, top right
	weight.setZero();
	double x0, x1, y0, y1;
	int nx, ny;
	nx = corner(0);
	ny = corner(1);
	weight.setZero();	
	x0  = vertex(0); //[0](0)
	x1  = vertex(1); //[1](0)
	y0  = vertex(2); //[0](1)
	y1  = vertex(3); //[2](1)

	// we can always shift the cell to unit cell with bottom left vertex at (0,0)
	double norm_qx = q(qi*3)*nx/2 - x0;
	double norm_qy = q(qi*3+1)*ny/2 - y0;

	weight << (1-norm_qx)*(1-norm_qy), (1- norm_qy)*norm_qx, (1-norm_qx)*norm_qy, norm_qx*norm_qy;
}	