#include <Eigen/Dense>
#include <EigenTypes.h>

//Input: 
// q position of the particles
// g gravity
// dt time step
void body_force(Eigen::VectorXd &q, Eigen::VectorXd &qdot, Eigen::Vector3d &g, double dt);