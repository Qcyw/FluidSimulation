#include <iostream>
#include <thread>

#include <assignment_setup.h>
// #include <visualization.h>

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <imgui/imgui_internal.h>

Eigen::VectorXd q, q_disp;
Eigen::VectorXd qdot;
Eigen::MatrixXd P1, P2;

Eigen::VectorXd offset;

double t = 0;
double dt = 1E-5;

bool simulating = true;

bool simulation_callback(){
    while(simulating){
        simulate(q, qdot, dt, t);
        t += dt;
    }
    return false;
}

int main(int argc, char **argv){
    std::cout<<"Start Simulation\n";
    std::cout << "Default mode is FLIP, 1 is PIC, 2 is 0.95 FLIP/PIC" << std::endl;
    igl::opengl::glfw::Viewer viewer;
    P1.resize(4,3);
    P2.resize(4,3);
    assignment_setup(argc, argv, q, qdot, offset, P1, P2);
    q_disp.resize(q.size());
    int k = q.size()/3;
    
    std::thread simulation_thread(simulation_callback);
    simulation_thread.detach();
    
    viewer.callback_post_draw = [&](igl::opengl::glfw::Viewer & )->bool
    {  
        q_disp = q + offset;
        viewer.data().add_edges(P1,P2,Eigen::RowVector3d(1,0,0));
        viewer.data().set_points(Eigen::Map<Eigen::MatrixXd>((q_disp).data(), 3, k).transpose(), Eigen::RowVector3d(1,1,1));
        return false;
    };
    viewer.launch();
    return 1;
}
