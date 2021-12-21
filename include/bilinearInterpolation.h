#include <Eigen/Dense>
#include <EigenTypes.h>
#include <bilinear_weight.h>
void bilinearInterpolation(Eigen::VectorXd &q, Eigen::VectorXd &qdot, Eigen::VectorXd &u, Eigen::VectorXd &v, const Eigen::Vector2i corner, double dx);