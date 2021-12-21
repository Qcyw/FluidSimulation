#include <pressure.h>
#include <solid_boundary.h>

void pressure(Eigen::MatrixXd &A, Eigen::VectorXd &f, Eigen::VectorXd &u, Eigen::VectorXd &v, const Eigen::Vector2i &corner, Eigen::RowVector3i &index, double rho, double dx, double dt){

	int adj_cell = 4;
	int nx, ny, i, j; 
	Eigen::MatrixXd B, P, D;
	Eigen::VectorXd qj;

	nx = corner(0);
	ny = corner(1);
	i = index(0);
	j = index(1);

	B.resize(1, adj_cell);
	D.resize(adj_cell, adj_cell+1);
	P.resize(adj_cell, adj_cell);
	qj.resize(adj_cell);
	A.resize(1, adj_cell+1);
	f.resize(1);

	D.setZero();
	B.setZero();
	B<< -1., 1., -1., 1.;
	qj << u(i + j*(nx+1)), u(i+1 + j*(nx+1)), v(i + j*(nx)), v(i + (j+1)*(nx));
	solid_boundary(P, index, corner);

	D << -1, 0, 1, 0, 0,
		0, 1, -1, 0, 0,
		0, 0, 1, -1, 0,
		0, 0, -1, 0, 1;
	D /= dx;
 	B /= dx;

	A = B*P*D;
	f = rho*B*P*qj/dt;
}