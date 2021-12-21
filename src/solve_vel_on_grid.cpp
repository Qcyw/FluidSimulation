#include <solve_vel_on_grid.h>

void solve_vel_on_grid( Eigen::VectorXd &P, Eigen::VectorXd &u, Eigen::VectorXd &v, const Eigen::RowVector2i corner, double rho, double dx, double dt){

	int adj_cell = 4;
	int nx, ny; 
	Eigen::MatrixXd D;
	Eigen::VectorXd pj, qj;
	Eigen::VectorXd Pu, Pv;

	Pu.resize(u.size());
	Pv.resize(v.size());
	D.resize(adj_cell, adj_cell+1);
	pj.resize(adj_cell+1);
	qj.resize(adj_cell); 
	nx = corner(0); 
	ny = corner(1);

	D.setZero();
	Pu.setZero();
	Pv.setZero();

	D << -1, 0, 1, 0, 0,
		0, 1, -1, 0, 0,
		0, 0, 1, -1, 0,
		0, 0, -1, 0, 1;
	D /= dx;
	
	for(int i = 1; i < nx-1; i ++){
		for(int j = 0; j < ny; j ++){
			// don't need the P j+1/j-1
			pj << P(i-1+j*nx), P(i+1+j*nx), P(i+j*nx), 0, 0;
			qj = D*pj; // divergence of p
			Pu(i+j*(nx+1)) += qj(0);
			Pu(i+1+j*(nx+1)) += qj(1);
		}
	}
	for(int j = 1; j < ny-1; j ++){
		for(int i = 0; i < nx; i ++){
			// Don't need P i-1/i+1
			// pj << P(i-1+j*nx), P(i+1+j*nx), P(i+j*nx), P(i+(j-1)*nx), P(i+(j+1)*nx); 
			pj << 0,0, P(i+j*nx), P(i+(j-1)*nx), P(i+(j+1)*nx); 
			qj = D*pj; // divergence of p
			Pv(i+j*(nx)) += qj(2);
			Pv(i+(j+1)*(nx)) += qj(3);
		}
	}
	u -= dt*Pu/rho;
	v -= dt*Pv/rho;
}