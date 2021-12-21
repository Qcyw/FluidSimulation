#include <assemble_pressure.h>

typedef Eigen::Triplet<double> T;
void assemble_pressure(Eigen::VectorXd &P, Eigen::VectorXd &u, Eigen::VectorXd &v, const Eigen::Vector2i corner, double rho, double dx, double dt){

	int nx, ny;  nx = corner(0), ny = corner(1);
	int n = nx*ny;
	std::vector<T> tripletList;
	Eigen::SparseMatrixd A;
	Eigen::VectorXd b, Ptmp; 

	A.resize(n,n);
	P.resize(n);
	Ptmp.resize(n);
	b.resize(n);

	A.setZero();
	P.setZero();
	Ptmp.setZero();
	b.setZero();


	int track = 0;
	int row = 0;
	int col = 0; 

	for(int j = 0; j < ny; j ++){
		for(int i = 0; i < nx; i++){	
			Eigen::RowVector3i index;
			index << i, j, 0;
			Eigen::MatrixXd Aj;
			Eigen::VectorXd fj;
			pressure(Aj, fj, u, v, corner, index, rho, dx, dt);

			if(col == 0 && row == 0){
				// bottom left
				tripletList.push_back(T(track, i+1 + j*nx, Aj(1)))	;
				tripletList.push_back(T(track, i+0 + j*nx, Aj(2)))	;
				tripletList.push_back(T(track, i + (j+1)*nx, Aj(4)));	
			}
			else if( row == 0 && col > 0 && col < nx - 1){
				// bottom edge
				tripletList.push_back(T(track, i-1 + j*nx, Aj(0)))	;
				tripletList.push_back(T(track, i+1 + j*nx, Aj(1)))	;
				tripletList.push_back(T(track, i+0 + j*nx, Aj(2)))	;
				tripletList.push_back(T(track, i + (j+1)*nx, Aj(4)));	
			}
			else if( row == 0 && col == nx - 1){
				// bottom right
				tripletList.push_back(T(track, i-1 + j*nx, Aj(0)))	;
				tripletList.push_back(T(track, i+0 + j*nx, Aj(2)))	;
				tripletList.push_back(T(track, i + (j+1)*nx, Aj(4)));	
			}
			else if( col == 0 && row > 0 && row < ny - 1){
				// left edge
				tripletList.push_back(T(track, i+1 + j*nx, Aj(1)))	;
				tripletList.push_back(T(track, i+0 + j*nx, Aj(2)))	;
				tripletList.push_back(T(track, i + (j-1)*nx, Aj(3)));	
				tripletList.push_back(T(track, i + (j+1)*nx, Aj(4)));	
			}
			else if( col == nx - 1 && row > 0 && row < ny - 1){
				// right edge 
				tripletList.push_back(T(track, i-1 + j*nx, Aj(0)))	;
				tripletList.push_back(T(track, i+0 + j*nx, Aj(2)))	;
				tripletList.push_back(T(track, i + (j-1)*nx, Aj(3)));	
				tripletList.push_back(T(track, i + (j+1)*nx, Aj(4)));	
			}
			else if( col == 0 && row == ny - 1 ){
				// top left
				tripletList.push_back(T(track, i+1 + j*nx, Aj(1)))	;
				tripletList.push_back(T(track, i+0 + j*nx, Aj(2)))	;
				tripletList.push_back(T(track, i + (j-1)*nx, Aj(3)));	
			}
			else if( col > 0 && col < nx - 1 && row == ny -1){
				// top edge
				tripletList.push_back(T(track, i-1 + j*nx, Aj(0)))	;
				tripletList.push_back(T(track, i+1 + j*nx, Aj(1)))	;
				tripletList.push_back(T(track, i+0 + j*nx, Aj(2)))	;
				tripletList.push_back(T(track, i + (j-1)*nx, Aj(3)));	
			}
			else if( col == nx - 1 && row == ny - 1){
				// top right
				tripletList.push_back(T(track, i-1 + j*nx, Aj(0)))	;
				tripletList.push_back(T(track, i+0 + j*nx, Aj(2)))	;
				tripletList.push_back(T(track, i + (j-1)*nx, Aj(3)));	
			}
			else{
				tripletList.push_back(T(track, i-1 + j*nx, Aj(0)))	;
				tripletList.push_back(T(track, i+1 + j*nx, Aj(1)))	;
				tripletList.push_back(T(track, i+0 + j*nx, Aj(2)))	;
				tripletList.push_back(T(track, i + (j-1)*nx, Aj(3)));	
				tripletList.push_back(T(track, i + (j+1)*nx, Aj(4)));	
			}
			b(track) = fj(0);
			track += 1;
			col += 1;	
		}
		col = 0;
		row += 1;
	}
	A.setFromTriplets(tripletList.begin(), tripletList.end());
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> cg;
	cg.compute(A);
	P = cg.solve(-b);
	std::cout << "#iterations:     " << cg.iterations() << std::endl;
	std::cout << "estimated error: " << cg.error()      << std::endl;

}