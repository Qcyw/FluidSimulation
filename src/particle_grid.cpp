#include <particle_grid.h>

void particle_grid(Eigen::VectorXd &q, Eigen::VectorXd &qdot, Eigen::VectorXd &u, Eigen::VectorXd &v, const Eigen::Vector2i corner, double dx){
	
	Eigen::Vector4d weight;
	Eigen::Vector4d vertex;
	Eigen::VectorXd u_weight, v_weight;

	int n = q.size()/3;
	int nx, ny, xr, yr; 
	double x, y, qdotx, qdoty, aqdotx, aqdoty;
	bool ul, vb;

	nx = corner(0);
	ny = corner(1);

	u.resize((nx+1)*ny);
	v.resize(nx*(ny+1));
	u_weight.resize((nx+1)*ny);
	v_weight.resize(nx*(ny+1));

	u.setZero();
	v.setZero();
	u_weight.setZero();
	v_weight.setZero();
	// std::cout << "here1" << std::endl;
	for(int i = 0; i < n; i ++){
		// iterate over all the particles in the grid cells
		x = q(i*3)*nx/2.;
		y = q(i*3+1)*ny/2.;
		qdotx = qdot(i*3);
		qdoty = qdot(i*3+1);
		aqdotx = abs(qdot(i*3));
		aqdoty = abs(qdot(i*3+1));
		xr = round(x);
		yr = round(y);
		
		// std::cout << "xr: " <<xr << " | "<<"xr: " <<yr << std::endl;
		if((xr != 0 && yr != 0) && (xr != nx && yr != ny) && (xr != 0 && yr != ny) && (yr != 0 && xr != nx)){
			//==================================== linear interpolation of u ====================================//
			if((xr == 0 || xr == nx)){

				if(xr <= x){
					vertex << xr, xr+1, yr-0.5, yr+0.5; 
					ul = false;
				}
				else{
					vertex << xr-1, xr, yr-0.5, yr+0.5; 
					ul = true;
				}
				bilinear_weight(q, qdot, vertex, weight, corner, i);
				if(ul == false){
					u(xr + (yr-1)*(nx+1)) += weight(0)*qdotx;
					u(xr + 1 + (yr-1)*(nx+1)) += weight(1)*qdotx;
					u(xr + (yr)*(nx+1)) += weight(2)*qdotx;
					u(xr + 1 + (yr)*(nx+1)) += weight(3)*qdotx;
				} 
				else{
					u(xr - 1 + (yr-1)*(nx+1)) += weight(0)*qdotx;
					u(xr + (yr-1)*(nx+1)) += weight(1)*qdotx;
					u(xr - 1 + yr*(nx+1)) += weight(2)*qdotx;
					u(xr + yr*(nx+1)) += weight(3)*qdotx;
				}
			}else{
				if(xr <= x){
					vertex << xr, xr+1, yr-0.5, yr+0.5; 
					ul = false;
				}
				else if(xr > x){
					vertex << xr-1, xr, yr-0.5, yr+0.5; 
					ul = true;
				}
				bilinear_weight(q, qdot, vertex, weight, corner, i);
				if(ul == false){
					u(xr + (yr-1)*(nx+1)) += weight(0)*qdotx;
					u(xr + 1 + (yr-1)*(nx+1)) += weight(1)*qdotx;
					u(xr + (yr)*(nx+1)) += weight(2)*qdotx;
					u(xr + 1 + (yr)*(nx+1)) += weight(3)*qdotx;
				} else{
					u(xr - 1 + (yr-1)*(nx+1)) += weight(0)*qdotx;
					u(xr + (yr-1)*(nx+1)) += weight(1)*qdotx;
					u(xr - 1 + yr*(nx+1)) += weight(2)*qdotx;
					u(xr + yr*(nx+1)) += weight(3)*qdotx;
				}
			}

			//==================================== linear interpolation of v ====================================//
			if((yr == 0 || yr == ny)){
				if(yr <= y){
					vertex << xr-0.5, xr+0.5, yr, yr+1; 
					vb = false;
				}
				else{
					vertex << xr-0.5, xr+0.5, yr-1, yr; 
					vb = true;
				}
				bilinear_weight(q, qdot, vertex, weight, corner, i);
				if(vb == false){
					v(xr - 1 + (yr)*(nx)) += weight(0)*qdoty;
					v(xr + (yr)*(nx)) += weight(1)*qdoty;
					v(xr - 1 + (yr+1)*(nx)) += weight(2)*qdoty;
					v(xr + (yr+1)*(nx)) += weight(3)*qdoty;
				} else{ 
					v(xr - 1 + (yr-1)*(nx)) += weight(0)*qdoty;
					v(xr + (yr-1)*(nx)) += weight(1)*qdoty;
					v(xr - 1 + yr*(nx)) += weight(2)*qdoty;
					v(xr + yr*(nx)) += weight(3)*qdoty;
				}	
			}else{
	 			if(yr <= y){
					vertex << xr-0.5, xr+0.5, yr, yr+1; 
					vb = false;
				}
				else if(yr > y){
					// take down cell
					vertex << xr-0.5, xr+0.5, yr-1, yr; 
					vb = true;
				}
				bilinear_weight(q, qdot, vertex, weight, corner, i);
				if(vb == false){
					v(xr - 1 + (yr)*(nx)) += weight(0)*qdoty;
					v(xr + (yr)*(nx)) += weight(1)*qdoty;
					v(xr - 1 + (yr+1)*(nx)) += weight(2)*qdoty;
					v(xr + (yr+1)*(nx)) += weight(3)*qdoty;
				} else{ 
					v(xr - 1 + (yr-1)*(nx)) += weight(0)*qdoty;
					v(xr + (yr-1)*(nx)) += weight(1)*qdoty;
					v(xr - 1 + yr*(nx)) += weight(2)*qdoty;
					v(xr + yr*(nx)) += weight(3)*qdoty;
				}	
			}
		}	
	}
}

