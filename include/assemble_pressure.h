#include <Eigen/Dense>
#include <EigenTypes.h>
#include <pressure.h>

void assemble_pressure(Eigen::VectorXd &P, Eigen::VectorXd &u, Eigen::VectorXd &v, const Eigen::Vector2i corner, double rho, double dx, double dt);