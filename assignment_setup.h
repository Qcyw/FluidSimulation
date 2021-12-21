#ifndef ASSIGNMENT_SETUP_H
#define ASSIGNMENT_SETUP_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>
// #include <visualization.h>
#include <init_state.h>
#include <advect.h>
#include <body_force.h>
#include <update_cell_type.h>
#include <particle_grid.h>
#include <assemble_pressure.h>
#include <free_surface_pressure.h>
#include <solve_vel_on_grid.h>
#include <grid_particle.h>
#include <stdlib.h>

Eigen::Vector3d gravity; 
Eigen::VectorXd u, v, u_prev, v_prev, P, cell_type;
Eigen::Vector2i corner;
int track;
int rest;
int mode;
int iter;
double rho;
bool debug, free_surf;
double dx;

inline void simulate(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double t){
   
        double max = 0;
        for(int i = 0;i < q.size(); i ++){
            if(max < abs(qdot(i))) max = abs(qdot(i));
        }
        if(max != 0) dt = dx/max;
        advect(q, qdot, corner, dt);
        body_force(q, qdot, gravity, dt);
        update_cell_type(cell_type, q, qdot, corner, dx);
        particle_grid(q, qdot, u, v, corner, dx);
        particle_grid(q, qdot, u_prev, v_prev, corner, dx);
        free_surface_pressure(P, q, qdot, u, v, cell_type, corner, rho, dx, dt);
        solve_vel_on_grid(P, u, v, corner, rho, dx, dt);
        grid_particle(q, qdot, u, v, u_prev, v_prev, corner, mode, 0.95, dx);
}


inline void assignment_setup(int argc, char **argv, Eigen::VectorXd &q, Eigen::VectorXd &qdot, Eigen::VectorXd &offset, Eigen::MatrixXd &P1, Eigen::MatrixXd &P2) {
    mode = 0;
    if(argc > 1)
    {
        if(strcmp(argv[1], "0")){mode = 0;}
        else if(strcmp(argv[1], "1")){mode = 1;}
        else if(strcmp(argv[1], "2")){mode = 2;}
    }

    rest = 0;
    track = 0; 
    rho = 10;
    corner.setZero();
    corner << 6, 6;
    dx = 2./corner(0);
    gravity << 0., -9.8, 0.;
    init_state(q, qdot, corner, P1, P2, offset);
}

#endif

