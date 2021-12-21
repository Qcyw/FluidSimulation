#include <update_cell_type.h>

void update_cell_type(Eigen::VectorXd &cell_type, Eigen::VectorXd &q, Eigen::VectorXd &qdot, const Eigen::Vector2i corner, double dx){
	Eigen::VectorXd Xbar;
	Eigen::VectorXd xi; 
	Eigen::VectorXd x;
	double barnorm, rbar, h, ri, k, phi;
	int nx, ny; nx = corner(0), ny = corner(1);

	Xbar.resize(3);
	xi.resize(3);
	x.resize(3);
	cell_type.resize(nx*ny);
	x.setZero();
	cell_type.setZero();
	Xbar.setZero();
	xi.setZero();

	ri = dx/2;
	h = 3*ri;
	barnorm = 0;
	phi = 0;
	

	for(int i = 0; i < nx; i ++){
		for(int j = 0; j < ny; j ++){
			phi = 0;
			x(0) = i*dx + dx/2;
			x(1) = j*dx + dx/2; 
			x(3) = 0;
			// blobbies(phi, q, qdot, x, dx);
			for(int l = 0; l < q.size()/3; l ++){
				xi = q.segment<3>(l*3);
				if((x-xi).norm() < h){
					k = pow((1-pow((x-xi).norm()/(h),2)),3);
					Xbar += k*xi;
					rbar += k*ri;
					barnorm += k;
				}
			}
			if(barnorm != 0){
				Xbar /= barnorm;
				rbar /= barnorm;		
			}
			phi = (x-Xbar).norm() - rbar;	
			// std::cout << i << "|" << j << ": " << phi << std::endl;
			if(phi < 0){
				cell_type(i + j*nx) = 1.;
			}
		}
	}
}