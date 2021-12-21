#include <Eigen/Dense>
#include <EigenTypes.h>

#include <solid_boundary.h>


void pressure(Eigen::MatrixXd &A, Eigen::VectorXd &f, Eigen::VectorXd &u, Eigen::VectorXd &v, const Eigen::Vector2i &corner, Eigen::RowVector3i &index, double rho, double dx, double dt);