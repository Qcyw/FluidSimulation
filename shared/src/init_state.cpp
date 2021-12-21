#include <init_state.h>

void init_state(Eigen::VectorXd &q, Eigen::VectorXd &qdot, const Eigen::Vector2i corner, Eigen::MatrixXd &P1, Eigen::MatrixXd &P2, Eigen::VectorXd &offset) {
    // ball of radius 2
    const double eps = 0.04;
    Eigen::RowVector3d pos;
    int nx, ny;
    nx = corner(0);
    ny = corner(1);
    int k = 100;
    q.resize(k/10*k/10*3);
    q.setZero();
    qdot.resize(k/10*k/10*3);
    qdot.setZero();
    offset.resize(q.size());
    offset.setZero();

    double pos_x, pos_y, pos_z;
    double init_pos_x, init_pos_y, init_pos_z; 
    pos_z = 0;
    init_pos_x = 1.1;
    init_pos_y = 1.1;
    pos_x = init_pos_x;
    pos_y = init_pos_y;
    int track;
    // 2d fluid
    for(int y = 0; y < k/10; y ++){
        for(int x = 0; x < k/10; x ++){
            pos << pos_x, pos_y, pos_z;
            q.segment<3>(x*3+3*k/10*y) = pos;
            pos_x += eps;   
            track += 1;
        }
        pos_x = init_pos_x;
        pos_y += eps;
    }

    P1 << -2, -2, 0,
          -2, -2, 0,
          -2, 2, 0, 
          2, -2, 0;

    P2 << -2, 2, 0,
          2, -2, 0,
          2, 2, 0,
          2, 2, 0;

    P1/=2;
    P2/=2;

    for(int i = 0; i < offset.size()/3; i++){
        offset(i*3) = -1;
        offset(i*3+1) = -1;
    }
}