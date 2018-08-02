//
// Vortexje -- Gmsh wing section construction example.
//
// Copyright (C) 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <cmath>
#include <iostream>

#include <vortexje/solver.hpp>
#include <vortexje/parameters.hpp>
#include <vortexje/surface-loaders/ply-surface-loader.hpp>
#include <vortexje/surface-writers/vtk-surface-writer.hpp>
#include <vortexje/field-writers/vtk-field-writer.hpp>
#include <vortexje/surface-builder.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

static const double pi = 3.141592653589793238462643383279502884;

class Board : public Surface
{
public:
    // Constructor:
    Board()
        : Surface("board")
    {
        // Create surface:      
        SurfaceBuilder surface_builder(*this);
        
        const double length = 0.1;
        const int n_layers = 3;
        
        vector<int> prev_nodes;
        
        // -0.03, 0.5,
        // -0.05, 0.2,
        // -0.05, 0.05,
        for (int i = 0; i < n_layers; i++) {
            vector<Vector3d, Eigen::aligned_allocator<Vector3d> > points;
            points.push_back(Vector3d(-0.03, 0.0, -0.05));
            points.push_back(Vector3d(0.13, 0.0, -0.05));
            points.push_back(Vector3d(0.23, 0.0, -0.05));
            for (int j = 0; j < (int) points.size(); j++)
                points[j](2) += i * length / (double) (n_layers - 1);
                 
            vector<int> nodes = surface_builder.create_nodes_for_points(points);
            
            if (i > 0)
                vector<int> airfoil_panels = surface_builder.create_panels_between_shapes(nodes, prev_nodes);
                
            prev_nodes = nodes;
        }

        surface_builder.finish();
    }
};

void load_fin(shared_ptr<LiftingSurface> surface, const std::string &file)
{    
    // Load Gmsh mesh file:
    PLYSurfaceLoader loader;
    loader.load(surface, file);

    // Set lifting surface metadata:
    // N.B.  The node and panel numbers must match those in the Gmsh file,
    // *except* that Gmsh counts from 1, whereas Vortexje counts from 0. 
    // The Gmsh node and panel numbers are one higher than the Vortexje numbers.
    int n_airfoils = 10;
    int n_nodes_per_airfoil_side = 10;
    
    surface->upper_nodes.resize(n_nodes_per_airfoil_side, n_airfoils);
    surface->lower_nodes.resize(n_nodes_per_airfoil_side, n_airfoils);
    for (int j = 0; j < n_airfoils; j++) {
        for (int i = 0; i < n_nodes_per_airfoil_side; i++) {
            surface->upper_nodes(i, j) = 2*j * n_nodes_per_airfoil_side + i;

            surface->lower_nodes(i, j) = ((2*j)+1) * n_nodes_per_airfoil_side + i;
        }
    }
    
    surface->upper_panels.resize(n_nodes_per_airfoil_side, n_airfoils - 1);
    surface->lower_panels.resize(n_nodes_per_airfoil_side, n_airfoils - 1);
    for (int j = 0; j < n_airfoils - 1; j++) {
        for (int i = 0; i < n_nodes_per_airfoil_side; i++) {
            surface->upper_panels(i, j) = 2*j * n_nodes_per_airfoil_side + i;
            
            surface->lower_panels(i, j) = ((2*j)+1) * n_nodes_per_airfoil_side + i;
        }
    }
    
    // Finish trailing edge setup:
    surface->finish_trailing_edge();
}

// Main:
int
main (int argc, char **argv)
{
    // Enable wake convection:
    Parameters::convect_wake = false;
    
    // Create lifting surface object:
    shared_ptr<LiftingSurface> fin(new LiftingSurface("main"));
    load_fin(fin, "ply-lifting-surface.ply");
    shared_ptr<LiftingSurface> mirror(new LiftingSurface("mirror"));
    load_fin(mirror, "ply-lifting-surface.ply");
    Eigen::Transform<double, 3, Eigen::Affine> m(Eigen::Scaling(1.0,-1.0,1.0));
    mirror->transform(m);
    
    // Prescribe angle of attack:
    double alpha = 5.0 / 180.0 * pi;
    fin->rotate(Vector3d::UnitY(), alpha);
    mirror->rotate(Vector3d::UnitY(), alpha);
    
    // Create surface body:
    shared_ptr<Body> body(new Body(string("fin-section")));
    body->add_lifting_surface(fin);
    //body->add_lifting_surface(mirror);
    shared_ptr<Board> board(new Board());
    body->add_non_lifting_surface(board);
    
    // Set up solver:
    Solver solver("ply-lifting-surface-log");
    solver.add_body(body);
    
    // 10m/s = 36kph
    Vector3d freestream_velocity(10, 0, 0);
    solver.set_freestream_velocity(freestream_velocity);
    
    // water density 997kg/m3
    double fluid_density = 997;
    solver.set_fluid_density(fluid_density);
    
    // Set up surface writer:
    VTKSurfaceWriter surface_writer;
    
    // Set up field writer:
    VTKFieldWriter field_writer(nstream::BINARY, nstream::FLOAT);
    
    // Run simulation:
    double t = 0.0;
    double dt = 0.01;
    int step_number = 0;
    
    solver.initialize_wakes(dt);
    while (step_number < 20) {
        // Solve:
        solver.solve(dt);
        
        // Log coefficients:
        solver.log(step_number, surface_writer);
        
        // Update wake:
        solver.update_wakes(dt);
        
        // Step time:
        t += dt;
        step_number++;
    }
        
        // Enable below to log the velocity field:
        field_writer.write_velocity_field(solver, "ply-lifting-surface-log/velocity-field.vtk",
                                          -0.03, 0.5,
                                          -0.05, 0.2,
                                          -0.05, 0.05,
                                          0.005, 0.005, 0.005);
    
    // Done:
    return 0;
}
