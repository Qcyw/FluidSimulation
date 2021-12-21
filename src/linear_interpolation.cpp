#include <linear_interpolation.h>

void linear_interpolation(double &alpha, Eigen::VectorXd &q, Eigen::VectorXd &qdot, Eigen::Vector3d &xf, Eigen::Vector3d &xa, double dx){
	
	double barnorm, rbar, h, ri, k, phi;
	double eta, tol;

	Eigen::VectorXd Xbar;
	Eigen::VectorXd xi; 
	Eigen::VectorXd x; 
	Xbar.resize(3);
	xi.resize(3);
	x.resize(3);
	x.setZero();
	Xbar.setZero();
	xi.setZero();

	eta = 1E-3;
	tol = 1E-3;
	ri = dx/2;
	h = 3*ri;
	barnorm = 0;
	phi = 0;
	
	for(int i = 0; i < 1000; i ++)
	{
		x = xf + (xa-xf)*eta;
		// blobbies
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

		eta += 1E-3;
		if(phi > 0)
		{
			x = xf - 2*(xa-xf)*eta;
			eta *= 0.5;
		}	
		if(abs(phi) < tol)
		{
			break;
		}
		alpha = (x-xf).norm()*eta;
	}	
	alpha = 1 - alpha;
}
