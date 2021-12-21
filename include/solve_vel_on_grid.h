#include <Eigen/Dense>
#include <EigenTypes.h>

void solve_vel_on_grid(Eigen::VectorXd &P, Eigen::VectorXd &u, Eigen::VectorXd &v, const Eigen::RowVector2i corner, double rho, double dx, double dt);