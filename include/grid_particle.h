#include <Eigen/Dense>
#include <EigenTypes.h>
#include <bilinearInterpolation.h>


// Input
// u, v:  the velocity field after pressure projection 
// u_prev, v_prev: velocity field before pressure projection saved after the particle_grid
void grid_particle(Eigen::VectorXd &q, Eigen::VectorXd &qdot, Eigen::VectorXd &u, Eigen::VectorXd &v, Eigen::VectorXd &u_prev, Eigen::VectorXd &v_prev, const Eigen::Vector2i corner, int mode, double alpha, double dx);