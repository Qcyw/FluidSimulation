#include <bilinearInterpolation.h>

void bilinearInterpolation(Eigen::VectorXd &q, Eigen::VectorXd &qdot, Eigen::VectorXd &u, Eigen::VectorXd &v, const Eigen::Vector2i corner, double dx){
	int nx, ny; nx = corner(0), ny =  corner(1);
	int n = q.size()/3;
	double x, y, qdotx, qdoty;
	int xr, yr;
	bool ul, vb;
	Eigen::Vector4d weight;
	Eigen::Vector4d vertex;

	for(int i = 0; i < n; i ++){
		Eigen::Vector3d nqdot;
		nqdot.setZero();
		
		x = q(i*3)*nx/2.;
		y = q(i*3+1)*ny/2.;
		qdotx = qdot(i*3);
		qdoty = qdot(i*3 + 1);
		xr = round(x);
		yr = round(y);
		vertex.setZero();
		// std::cout << "iteration: " << i << "---" << xr << "|" << yr <<"|" << x <<"|" << y << std::endl;
		
		if(xr != 0 && yr != 0 && xr != nx && yr != ny){

			//==================================== linear interpolation of u ====================================//
			if(xr <= x){
				// take right cell
				// vertex << xr, yr - 0.5, 0,
				// 	xr + 1, yr - 0.5, 0, 
				// 	xr, yr + 0.5, 0, 
				// 	xr + 1, yr + 0.5, 0;
				vertex << xr, xr+1, yr-0.5, yr+0.5; 
				ul = false;
			}
			else if(xr > x){
				// left cell
				// vertex << xr - 1, yr - 0.5, 0, 
				// 	xr, yr - 0.5, 0, 
				// 	xr - 1, yr + 0.5, 0, 
				// 	xr, yr + 0.5, 0;
				vertex << xr-1, xr, yr-0.5, yr+0.5; 
				ul = true;
			}
			bilinear_weight(q, qdot, vertex, weight, corner, i);
			if(ul == false){
				nqdot(0) = u(xr + (yr-1)*(nx+1))*weight(0) + 
					u(xr + 1 + (yr-1)*(nx+1))*weight(1) + 
					u(xr + (yr)*(nx+1))*weight(2) + 
					u(xr + 1 + (yr)*(nx+1))*weight(3);
			} else{
				nqdot(0) = u(xr - 1 + (yr-1)*nx)*weight(0) + 
					u(xr + (yr-1)*nx)*weight(1) + 
					u(xr - 1 + yr*nx)*weight(2) + 
					u(xr + yr*nx)*weight(3);
			}

			//==================================== linear interpolation of v ====================================//
			if(yr <= y){
				// take top cell
				// vertex << xr - 0.5, yr, 0,
				// 		xr + 0.5, yr, 0,
				// 		xr - 0.5, yr + 1, 0,
				// 		xr + 0.5, yr + 1, 0;
				vertex << xr-0.5, xr+0.5, yr, yr+1; 
				vb = false;
			}
			else if(yr > y){
				// take down cell
				// vertex << xr - 0.5, yr - 1, 0,
				// 		xr + 0.5, yr - 1, 0,
				// 		xr - 0.5, yr, 0,
				// 		xr + 0.5, yr, 0;
				vertex << xr-0.5, xr+0.5, yr-1, yr; 
				vb = true;
			}
			bilinear_weight(q, qdot, vertex, weight, corner, i);

			if(vb == false){
				nqdot(1) = 
					v(xr - 1 + (yr)*(nx))*weight(0)+
					v(xr + (yr)*(nx))*weight(1)+
					v(xr - 1 + (yr + 1)*(nx))*weight(2)+
					v(xr + (yr + 1)*(nx))*weight(3);
			} else{
				nqdot(1) = 
					v(xr - 1 + (yr-1)*(nx))*weight(0)+
					v(xr + (yr-1)*(nx))*weight(1)+
					v(xr - 1 + yr*(nx))*weight(2)+
					v(xr + yr*(nx))*weight(3);
			}
			qdot.segment<3>(i*3) = nqdot;
		}
		// std::cout << "iteration: " << i <<", bilinear interpolation: \n " << nqdot  << std::endl;
	}
}