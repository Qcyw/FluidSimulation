#include <Eigen/Dense>
#include <EigenTypes.h>

// Input:

// q the position
// qdot velocity
// u velocity component in x
// v velocity component in y
// vertex the value at each vertex bottom left, bottom right, top left, top right
// pi the particle interested in

// output: weight value for each particle
void bilinear_weight(Eigen::VectorXd &q, Eigen::VectorXd &qdot, Eigen::Vector4d &vertex, Eigen::Vector4d &weight, const Eigen::Vector2i corner, int qi);