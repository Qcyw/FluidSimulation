#include <Eigen/Dense>
#include <EigenTypes.h>

// identify cell type naively, if there is a particle then it's fluid with entry 1, air cell with entry 0. 
void update_cell_type(Eigen::VectorXd &cell_type, Eigen::VectorXd &q, Eigen::VectorXd &qdot, const Eigen::Vector2i corner, double dx);