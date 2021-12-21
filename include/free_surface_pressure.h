#include <Eigen/Dense>
#include <EigenTypes.h>
#include <linear_interpolation.h>
#include <pressure.h>
void free_surface_pressure(Eigen::VectorXd &P, Eigen::VectorXd &q, Eigen::VectorXd &qdot, Eigen::VectorXd &u, Eigen::VectorXd &v, Eigen::VectorXd &cell_type, const Eigen::Vector2i corner, double rho, double dx, double dt);