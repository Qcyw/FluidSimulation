
#include <advect.h>

void advect(Eigen::VectorXd &q, Eigen::VectorXd &qdot, const Eigen::Vector2i corner, double dt)
{     
    double x, y, eps, x_prev, y_prev;
    int nx, ny;
    nx = corner(0);
    ny = corner(1);
    double dx = 2/nx;
    eps = dx/4;

    for(int i = 0; i < q.size()/3; i ++){
        q.segment<3>(i*3) += qdot.segment<3>(i*3)*dt;
        x_prev = q(i*3);
        y_prev = q(i*3 + 1);
        x = q(i*3);
        y = q(i*3 + 1);

        if( x < eps ){
            q(i*3) = eps*2;
        }
        if( y < eps ){
            q(i*3+1) = eps*2;
        }
        if( x > 2.0 - eps ){
            q(i*3) = 2 - eps*2;
        }
        if( y > 2.0 -eps ){
            q(i*3+1) = 2 - eps*2;
        }
    }
}