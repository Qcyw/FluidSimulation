#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

void init_state(Eigen::VectorXd &q, Eigen::VectorXd &qdot, const Eigen::Vector2i corner, Eigen::MatrixXd &P1, Eigen::MatrixXd &P2, Eigen::VectorXd &offset);