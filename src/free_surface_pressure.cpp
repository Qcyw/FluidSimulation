#include<free_surface_pressure.h>

typedef Eigen::Triplet<double> T;
void free_surface_pressure(Eigen::VectorXd &P, Eigen::VectorXd &q, Eigen::VectorXd &qdot, Eigen::VectorXd &u, Eigen::VectorXd &v, Eigen::VectorXd &cell_type, const Eigen::Vector2i corner, double rho, double dx, double dt){
	
	int adj_cell = 4;
	int nx, ny;  nx = corner(0), ny = corner(1);
	int n = nx*ny;
	double alpha;
	std::vector<T> tripletList;
	Eigen::SparseMatrixd A;
	Eigen::VectorXd b; 
	Eigen::Vector3d xf, xa;
	// Eigen::MatrixXd  B, bdry, Aj;
	// Eigen::RowVector3i index;
			
	// qj.resize(adj_cell);
	// B.resize(1, adj_cell);

	A.resize(n,n);
	P.resize(n);
	b.resize(n);
	// Aj.resize(1, adj_cell+1);
	// bdry.resize(adj_cell, adj_cell);
	// D.resize(adj_cell, adj_cell+1);

	A.setZero();
	P.setZero();
	b.setZero();
	// Aj.setZero();
	// D.setZero();
	// B.setZero();

	// B<< -1., 1., -1., 1.;
	// D << -1, 0, 1, 0, 0,
	// 	0, 1, -1, 0, 0,
	// 	0, 0, 1, -1, 0,
	// 	0, 0, -1, 0, 1;
	// D /= dx;
	
	int track = 0;
	int row = 0;
	int col = 0; 

	for(int j = 0; j < ny; j++){
		for(int i = 0; i < nx; i ++){
			Eigen::RowVector3i index;
			index << i, j, 0;			
			Eigen::MatrixXd Aj;
			Eigen::VectorXd fj;
			pressure(Aj, fj, u, v, corner, index, rho, dx, dt);

			// if(cell_type(i + j*nx) == 1)
			// {
				xf << (i)*dx + dx/2, j*dx + dx/2, 0; 
				if(col == 0 && row == 0){
					// bottom left
					tripletList.push_back(T(track, i+1 + j*nx, Aj(1)))	;
					tripletList.push_back(T(track, i+0 + j*nx, Aj(2)))	;
					tripletList.push_back(T(track, i + (j+1)*nx, Aj(4)));	
					
					if(cell_type(col + 1 + row*(nx)) == 0)
					{
						xa << (i+1)*dx + dx/2, (j)*dx + dx/2, 0; 
						linear_interpolation(alpha, q, qdot, xf, xa, dx);
						tripletList.push_back(T(track, i + j*nx, (alpha/(1-alpha))*Aj(2)));
						// tripletList.push_back(T(track, i + 1 + (j)*nx, -Aj(1)));
					} 
					if(cell_type(col + (row+1)*(nx)) == 0)
					{
						xa << (i)*dx + dx/2, (j+1)*dx + dx/2, 0; 
						linear_interpolation(alpha, q, qdot, xf, xa, dx);
						// tripletList.push_back(T(track, i + (j+1)*nx, -Aj(4)));
						tripletList.push_back(T(track, i + j*nx, (alpha/(1-alpha))*Aj(2)));
					} 
				}
				else if( row == 0 && col > 0 && col < nx - 1){
					// bottom edge
					tripletList.push_back(T(track, i-1 + j*nx, Aj(0)))	;
					tripletList.push_back(T(track, i+1 + j*nx, Aj(1)))	;
					tripletList.push_back(T(track, i+0 + j*nx, Aj(2)))	;
					tripletList.push_back(T(track, i + (j+1)*nx, Aj(4)));	
				
					if(cell_type(col - 1 + row*nx) == 0){
						xa << (i-1)*dx + dx/2, j*dx + dx/2, 0;
						linear_interpolation(alpha, q, qdot, xf, xa, dx);
						tripletList.push_back(T(track, i + j*nx, (alpha/(1-alpha))*Aj(2)));
						// tripletList.push_back(T(track, i-1 + j*nx, -Aj(0)));
					}
					if(cell_type(col + 1 + row*nx) == 0){
						xa << (i+1)*dx + dx/2, j*dx + dx/2, 0;
						linear_interpolation(alpha, q, qdot, xf, xa, dx);
						tripletList.push_back(T(track, i + j*nx, (alpha/(1-alpha))*Aj(2)));
						// tripletList.push_back(T(track, i+1 + j*nx, -Aj(1)));
					}
					if(cell_type(col + (row + 1)*nx) == 0){
						xa << (i)*dx + dx/2, (j+1)*dx + dx/2, 0;
						linear_interpolation(alpha, q, qdot, xf, xa, dx);
						tripletList.push_back(T(track, i + j*nx, (alpha/(1-alpha))*Aj(2)));
						// tripletList.push_back(T(track, i + (j+1)*nx, -Aj(4)));
					}
				}
				else if( row == 0 && col == nx - 1){
					// bottom right
					tripletList.push_back(T(track, i-1 + j*nx, Aj(0)))	;
					tripletList.push_back(T(track, i+0 + j*nx, Aj(2)))	;
					tripletList.push_back(T(track, i + (j+1)*nx, Aj(4)));	
					if(cell_type(col - 1 + row*nx) == 0){
						xa << (i-1)*dx + dx/2, j*dx + dx/2, 0;
						linear_interpolation(alpha, q, qdot, xf, xa, dx);
						tripletList.push_back(T(track, i + j*nx, (alpha/(1-alpha))*Aj(2)));
						// tripletList.push_back(T(track, i-1 + j*nx, -Aj(0)));
					}
					if(cell_type(col + (row+1)*nx) == 0){
						xa << (i)*dx + dx/2, (j+1)*dx + dx/2, 0;
						linear_interpolation(alpha, q, qdot, xf, xa, dx);
						tripletList.push_back(T(track, i + j*nx, (alpha/(1-alpha))*Aj(2)));
						// tripletList.push_back(T(track, i + (j+1)*nx, -Aj(4)));
					}
				}
				else if( col == 0 && row > 0 && row < ny - 1){
					// left edge
					tripletList.push_back(T(track, i+1 + j*nx, Aj(1)))	;
					tripletList.push_back(T(track, i+0 + j*nx, Aj(2)))	;
					tripletList.push_back(T(track, i + (j-1)*nx, Aj(3)));	
					tripletList.push_back(T(track, i + (j+1)*nx, Aj(4)));

					if(cell_type(col + 1 + row*nx) == 0){
						xa << (i+1)*dx + dx/2, j*dx + dx/2, 0;
						linear_interpolation(alpha, q, qdot, xf, xa, dx);
						tripletList.push_back(T(track, i + j*nx, (alpha/(1-alpha))*Aj(2)));
						// tripletList.push_back(T(track, i+1 + j*nx, -Aj(1)));
					}
					if(cell_type(col + (row-1)*nx) == 0){
						xa << (i)*dx + dx/2, (j-1)*dx + dx/2, 0;
						linear_interpolation(alpha, q, qdot, xf, xa, dx);
						tripletList.push_back(T(track, i + j*nx, (alpha/(1-alpha))*Aj(2)));
						// tripletList.push_back(T(track, i + (j-1)*nx, -Aj(3)));
					}
					if(cell_type(col + (row + 1)*nx) == 0){
						xa << (i)*dx + dx/2, (j+1)*dx + dx/2, 0;
						linear_interpolation(alpha, q, qdot, xf, xa, dx);
						tripletList.push_back(T(track, i + j*nx, (alpha/(1-alpha))*Aj(2)));
						// tripletList.push_back(T(track, i + (j+1)*nx, -Aj(4)));
					}	
				}
				else if( col == nx - 1 && row > 0 && row < ny - 1){
					// right edge 
					tripletList.push_back(T(track, i-1 + j*nx, Aj(0)))	;
					tripletList.push_back(T(track, i+0 + j*nx, Aj(2)))	;
					tripletList.push_back(T(track, i + (j-1)*nx, Aj(3)));	
					tripletList.push_back(T(track, i + (j+1)*nx, Aj(4)));

					if(cell_type(col - 1 + row*nx) == 0){
						xa << (i-1)*dx + dx/2, j*dx + dx/2, 0;
						linear_interpolation(alpha, q, qdot, xf, xa, dx);
						tripletList.push_back(T(track, i + j*nx, (alpha/(1-alpha))*Aj(2)));
						// tripletList.push_back(T(track, i-1 + j*nx, -Aj(0)));
					}
					if(cell_type(col + (row-1)*nx) == 0){
						xa << (i)*dx + dx/2, (j-1)*dx + dx/2, 0;
						linear_interpolation(alpha, q, qdot, xf, xa, dx);
						tripletList.push_back(T(track, i + j*nx, (alpha/(1-alpha))*Aj(2)));
						// tripletList.push_back(T(track, i + (j-1)*nx, -Aj(3)));
					}
					if(cell_type(col + (row + 1)*nx) == 0){
						xa << (i)*dx + dx/2, (j+1)*dx + dx/2, 0;
						linear_interpolation(alpha, q, qdot, xf, xa, dx);
						tripletList.push_back(T(track, i + j*nx, (alpha/(1-alpha))*Aj(2)));
						// tripletList.push_back(T(track, i + (j+1)*nx, -Aj(4)));
					}
				}
				else if( col == 0 && row == ny - 1 ){
					// top left
					tripletList.push_back(T(track, i+1 + j*nx, Aj(1)))	;
					tripletList.push_back(T(track, i+0 + j*nx, Aj(2)))	;
					tripletList.push_back(T(track, i + (j-1)*nx, Aj(3)));	

					if(cell_type(col + 1 + row*(nx)) == 0)
					{
						xa << (i+1)*dx + dx/2, (j)*dx + dx/2, 0; 
						linear_interpolation(alpha, q, qdot, xf, xa, dx);
						tripletList.push_back(T(track, i + j*nx, (alpha/(1-alpha))*Aj(2)));
						// tripletList.push_back(T(track, i + 1 + (j)*nx, -Aj(1)));
					} 
					if(cell_type(col + (row-1)*(nx)) == 0)
					{
						xa << (i)*dx + dx/2, (j-1)*dx + dx/2, 0; 
						linear_interpolation(alpha, q, qdot, xf, xa, dx);
						// tripletList.push_back(T(track, i + (j-1)*nx, -Aj(3)));
						tripletList.push_back(T(track, i + j*nx, (alpha/(1-alpha))*Aj(2)));
					} 
				}
				else if( col > 0 && col < nx - 1 && row == ny -1){
					// top edge
					tripletList.push_back(T(track, i-1 + j*nx, Aj(0)))	;
					tripletList.push_back(T(track, i+1 + j*nx, Aj(1)))	;
					tripletList.push_back(T(track, i+0 + j*nx, Aj(2)))	;
					tripletList.push_back(T(track, i + (j-1)*nx, Aj(3)));	

					if(cell_type(col - 1 + row*nx) == 0){
						xa << (i-1)*dx + dx/2, j*dx + dx/2, 0;
						linear_interpolation(alpha, q, qdot, xf, xa, dx);
						tripletList.push_back(T(track, i + j*nx, (alpha/(1-alpha))*Aj(2)));
						// tripletList.push_back(T(track, i-1 + j*nx, -Aj(0)));
					}
					if(cell_type(col + (row-1)*nx) == 0){
						xa << (i)*dx + dx/2, (j-1)*dx + dx/2, 0;
						linear_interpolation(alpha, q, qdot, xf, xa, dx);
						tripletList.push_back(T(track, i + j*nx, (alpha/(1-alpha))*Aj(2)));
						// tripletList.push_back(T(track, i + (j-1)*nx, -Aj(3)));
					}
					if(cell_type(col + 1+ (row)*nx) == 0){
						xa << (i+1)*dx + dx/2, (j)*dx + dx/2, 0;
						linear_interpolation(alpha, q, qdot, xf, xa, dx);
						tripletList.push_back(T(track, i + j*nx, (alpha/(1-alpha))*Aj(2)));
						// tripletList.push_back(T(track, i + 1 + (j)*nx, -Aj(1)));
					}
				}
				else if( col == nx - 1 && row == ny - 1){
					// top right
					tripletList.push_back(T(track, i-1 + j*nx, Aj(0)))	;
					tripletList.push_back(T(track, i+0 + j*nx, Aj(2)))	;
					tripletList.push_back(T(track, i + (j-1)*nx, Aj(3)));

					if(cell_type(col - 1 + row*nx) == 0){
						xa << (i-1)*dx + dx/2, j*dx + dx/2, 0;
						linear_interpolation(alpha, q, qdot, xf, xa, dx);
						tripletList.push_back(T(track, i + j*nx, (alpha/(1-alpha))*Aj(2)));
						// tripletList.push_back(T(track, i-1 + j*nx, -Aj(0)));
					}
					if(cell_type(col + (row-1)*nx) == 0){
						xa << (i)*dx + dx/2, (j-1)*dx + dx/2, 0;
						linear_interpolation(alpha, q, qdot, xf, xa, dx);
						tripletList.push_back(T(track, i + j*nx, (alpha/(1-alpha))*Aj(2)));
						// tripletList.push_back(T(track, i + (j-1)*nx, -Aj(3)));
					}
				
				}
				else{
					tripletList.push_back(T(track, i-1 + j*nx, Aj(0)))	;
					tripletList.push_back(T(track, i+1 + j*nx, Aj(1)))	;
					tripletList.push_back(T(track, i+0 + j*nx, Aj(2)))	;
					tripletList.push_back(T(track, i + (j-1)*nx, Aj(3)));	
					tripletList.push_back(T(track, i + (j+1)*nx, Aj(4)));	

					if(cell_type(col - 1 + row*nx) == 0){
						xa << (i-1)*dx + dx/2, j*dx + dx/2, 0;
						linear_interpolation(alpha, q, qdot, xf, xa, dx);
						tripletList.push_back(T(track, i + j*nx, (alpha/(1-alpha))*Aj(2)));
						// tripletList.push_back(T(track, i-1 + j*nx, -Aj(0)));
					}
					if(cell_type(col + 1+ (row)*nx) == 0){
						xa << (i+1)*dx + dx/2, (j)*dx + dx/2, 0;
						linear_interpolation(alpha, q, qdot, xf, xa, dx);
						tripletList.push_back(T(track, i + j*nx, (alpha/(1-alpha))*Aj(2)));
						// tripletList.push_back(T(track, i + 1 + (j)*nx, -Aj(1)));
					}
					if(cell_type(col + (row-1)*nx) == 0){
						xa << (i)*dx + dx/2, (j-1)*dx + dx/2, 0;
						linear_interpolation(alpha, q, qdot, xf, xa, dx);
						tripletList.push_back(T(track, i + j*nx, (alpha/(1-alpha))*Aj(2)));
						// tripletList.push_back(T(track, i + (j-1)*nx, -Aj(3)));
					}
					if(cell_type(col + (row + 1)*nx) == 0){
						xa << (i)*dx + dx/2, (j+1)*dx + dx/2, 0;
						linear_interpolation(alpha, q, qdot, xf, xa, dx);
						tripletList.push_back(T(track, i + j*nx, (alpha/(1-alpha))*Aj(2)));
						// tripletList.push_back(T(track, i + (j+1)*nx, -Aj(4)));
					}
				}
			// }
			b(track) = fj(0);
			// std::cout << track << std::endl;
			track += 1;
			col += 1;	
		}
		col = 0;
		row += 1;	
	}

	A.setFromTriplets(tripletList.begin(), tripletList.end());
	// std::cout << A <<std::endl;
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> cg;
	cg.compute(A);
	P = cg.solve(-b);
	for(int i = 0; i < cell_type.size(); i ++){
		if(cell_type(i) == 0){
			P(i) = 0;
		}
	}
	
}