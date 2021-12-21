#include <Eigen/Dense>
#include <EigenTypes.h>

// Input:
// u - the velocity field
// dt - time step
// q - current field quantity
void advect(Eigen::VectorXd &q, Eigen::VectorXd &qdot, const Eigen::Vector2i corner, double dt);

