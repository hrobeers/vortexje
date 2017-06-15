//
// Vortexje -- Create Elliptic wing with NACA0012 airfoil gmsh file.
//
// Copyright (C) 2017 hrobeers.
//
// Authors: Hans Robeers <hansrobeers@gmail.com>
//

#include <cmath>

#include <vortexje/solver.hpp>
#include <vortexje/lifting-surface-builder.hpp>
#include <vortexje/shape-generators/airfoils/naca4-airfoil-generator.hpp>
#include <vortexje/surface-writers/gmsh-surface-writer.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

static const double pi = 3.141592653589793238462643383279502884;

// Cosine rule:
static double
cosine_rule(int n_points, int i)
{
    return 0.5 * (1 - cos(pi * i / (double) n_points));
}

int
main (int argc, char **argv)
{
    // Needed to avoid creating extra vertices
    Parameters::zero_threshold = 1;

    // Create wing:
    shared_ptr<Surface> wing(new Surface("main"));

    SurfaceBuilder surface_builder(*wing);

    const double AR = 6.0;
    const double chord = 0.5;
    const double span = AR * chord;

//    const int n_points_per_airfoil = 40;
//    const int n_airfoils = AR * n_points_per_airfoil / 2;

    int n_points_per_airfoil = 32;
    int n_airfoils = 21;

    int trailing_edge_point_id;
    vector<int> prev_airfoil_nodes;

    vector<vector<int> > node_strips;
    vector<vector<int> > panel_strips;

    for (int i = 0; i < n_airfoils; i++) {
        double y        = -span / 2.0 + span * cosine_rule(n_airfoils - 1, i);
        double chord_y  = chord * sqrt(1.0 - pow(2.0 * y / span, 2));
        double offset_x = (chord - chord_y) / 4.0;

        vector<Vector3d, Eigen::aligned_allocator<Vector3d> > airfoil_points;

//        if (i == 0 || i == n_airfoils - 1)
//            chord_y = span * 0.01;
        if (i == 0 || i == n_airfoils - 1) {
            Vector3d tip_point(0.0, 0.0, 0.0);
            for (int j = 0; j < n_points_per_airfoil; j++)
                airfoil_points.push_back(tip_point);

        } else
            airfoil_points = NACA4AirfoilGenerator::generate(0, 0, 0.12, true, chord_y, n_points_per_airfoil, trailing_edge_point_id);

        for (int j = 0; j < (int) airfoil_points.size(); j++) {
            airfoil_points[j](0) += offset_x;
            airfoil_points[j](2) += y;
        }

        vector<int> airfoil_nodes = surface_builder.create_nodes_for_points(airfoil_points);
        node_strips.push_back(airfoil_nodes);

        if (i > 0) {
            vector<int> airfoil_panels = surface_builder.create_panels_between_shapes(airfoil_nodes, prev_airfoil_nodes, trailing_edge_point_id);
            panel_strips.push_back(airfoil_panels);
        }

        prev_airfoil_nodes = airfoil_nodes;
    }

//    surface_builder.finish(node_strips, panel_strips, trailing_edge_point_id);

    string filename = "elliptic-wing.msh";

    GmshSurfaceWriter writer;
    const std::vector<std::string> view_names;
    const vector<MatrixXd, Eigen::aligned_allocator<MatrixXd> > view_data;
    writer.sort_vertices = true;
    writer.write(wing, filename, 0, 0, view_names, view_data);

    // Done.
    return 0;
}
