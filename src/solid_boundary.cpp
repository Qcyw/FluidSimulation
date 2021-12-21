#include<solid_boundary.h>

void solid_boundary(Eigen::MatrixXd &P, Eigen::RowVector3i &index, const Eigen::Vector2i &corner){
	int adj_cell = 4;
	int nx, ny, i, j;; 
	nx = corner(0);
	ny = corner(1);
	i = index(0);
	j = index(1);

	P.resize(adj_cell, adj_cell);
	P.setIdentity();
	
	if(i == 0){ P(0,0) = 0;}
	if(j == 0) {P(2,2) = 0;}
	if(i == nx-1) {P(1,1) = 0;}
	if(j == ny-1) {P(3,3) = 0;}
}