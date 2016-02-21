/**
* @file mass_spring.cpp
* Implementation of mass-spring system using Graph
*
* @brief Reads in two files specified on the command line.
* First file: 3D Points (one per line) defined by three doubles
* Second file: Tetrahedra (one per line) defined by 4 indices into the point
* list
*/

#include <fstream>
#include "CME212/SDLViewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"
#include "Graph.hpp"

// Gravity in meters/sec^2
static constexpr double grav = 9.81;


/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
};


/** Custom structure of data to store with Edges */
struct EdgeData {
  double L; // edge length
};


// Define your Graph type
typedef Graph<NodeData, EdgeData> GraphType;
typedef typename GraphType::node_type Node;
typedef typename GraphType::edge_type Edge;
struct FixedConstraint;

/** Change a graph's nodes according to a step of the symplectic Euler
*    method with the given node force.
* @param[in,out] g      Graph
* @param[in]     t      The current time (useful for time-dependent forces)
* @param[in]     dt     The time step
* @param[in]     force  Function object defining the force per node
* @return the next time step (usually @a t + @a dt)
*
* @tparam G::node_value_type supports ???????? YOU CHOOSE
* @tparam F is a function object called as @a force(n, @a t),
*           where n is a node of the graph and @a t is the current time.
*           @a force must return a Point representing the force vector on Node
*           at time @a t.
*/
template <typename G, typename F, typename C = FixedConstraint>
double symp_euler_step(G& g, double t, double dt, F force, C c = C()) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  // Set constraint
  c(g,t);
  return t + dt;
}


/** Force Function. */
struct Problem1Force {
  /** The force applied to @a n at time @a t combines
  *     mass-spring force and gravity.
  *  Points at (0, 0, 0) and (1, 0, 0) can't move. The forces applied
  *     to these points are always zero. */
  double K = 100;
  Point operator()(Node n, double t) {
    (void) t;     // silence compiler warnings

    // Zero force at (0,0,0) and (1,0,0).
    if (n.position() == Point(0 ,0 ,0) || n.position() == Point(1 ,0 ,0)) {
      return Point(0,0,0);
    } // Forces for other nodes
    else {
      Point springForce = Point(0,0,0);
      Point this_pos = n.position();
      for (GraphType::IncidentIterator it = n.edge_begin(); it != n.edge_end(); ++it) {
        Point incidentNode = (*it).node2().position();
        Point p_diff = this_pos - incidentNode;
        double dist = norm(p_diff);
        double L = (*it).value().L;
        springForce += (-1.0)* K * (p_diff)/dist * (dist - L);
      }
      return springForce + n.value().mass*Point(0,0,-1.0*grav);
    }
  }
};


/** Functor for Gravity Force */
struct GravityForce {
  double grav_;

  /** GravityForce Constructor.
  * @param[in] g Gravity with unit m/s^2.
  */
  GravityForce(double g = grav) : grav_(g) {}

  /** Compute Gravity Force
  * @param[in] n Valid node.
  * @param[in] t Valid time.
  * @return Point object, representing the gravity force F_g = m*g.
  */
  Point operator()(Node n, double t) {
    (void) t;     // silence compiler warnings
    return n.value().mass * Point(0,0,-1.0 * grav_);
  }
};


/** Functor for Spring Force */
struct MassSpringForce {
  double K_;
  /** MassSpringForce Constructor.
  * @param[in] K Spring constant with unit N/m.
  */
  MassSpringForce(double K = 100) : K_(K) {}

  /** Computes Mass Spring Force
  * @param[in] n Valid node.
  * @param[in] t Valid time.
  * @return Point object, representing the mass spring force.
  */
  Point operator()(Node n, double t) {
    (void) t;     // silence compiler warnings

    Point springForce = Point(0,0,0);
    Point this_pos = n.position();
    for (GraphType::IncidentIterator it = n.edge_begin(); it != n.edge_end(); ++it) {
      Point incidentNode = (*it).node2().position();
      Point p_diff = this_pos - incidentNode;
      double dist = norm(p_diff);
      double L = (*it).value().L;
      springForce += (-1.0)* K_ * (p_diff)/dist * (dist - L);
    }
    return springForce;
  }
};


/** Functor for Damping Force.
* @param[in] d_const Damping constant with unit N*s/m.
*/
struct DampingForce {
  double d_const_;

  /** DampingForce Constructor.
  * @param[in] d_const Damping constant with unit N*s/m
  */
  DampingForce(double d_const) : d_const_(d_const) {}

  /** Computes Damping Force
  * @param[in] n Valid node.
  * @param[in] t Valid time.
  * @return Point object, representing the damping force F = -d_const * v.
  */
  Point operator()(Node n, double t) {
    (void) t;     // silence compiler warnings
    return -1.0 * d_const_ * n.value().vel;
  }
};


/** Functor to combine two forces.
* @param[in] f1, f2 valid forces.
*/
template <typename F1, typename F2>
struct CombinedForce {
  F1 f1_;
  F2 f2_;

  /** CombinedForce Constructor.
  * @param[in] f1 First valid force.
  * @param[in] f2 Second valid force.
  */
  CombinedForce(F1 f1, F2 f2) : f1_(f1), f2_(f2) {}

  /** Compute Combined Forces
  * @param[in] n Valid node.
  * @param[in] t Valid time.
  * @return Point object that represents the sum of @a f1_ and @a f2_.
  */
  Point operator()(Node n, double t){
    return f1_(n,t) + f2_(n,t);
  }
};


/** Returns a combination of arbitrary forces, specified by template arguments.
* @param[in] f1, f2 valid forces
* @pre Valid forces that take as input arguments a node and a time.
* @return A CombinedForce object representing the sum of the given forces.
*/
template <typename F1, typename F2>
CombinedForce<F1, F2> make_combined_force(F1 f1, F2 f2){
  return CombinedForce<F1, F2> ({f1, f2});
};

/** Combine Force Function that returns a combination of forces
* @param[in] f1, f2, and f3 are valid forces.
* @pre Valid forces that take as input arguments a node and a time.
* @return A CombinedForce object representing the sum of the given forces.
*/
template <typename F1, typename F2, typename F3>
CombinedForce<CombinedForce<F1, F2>, F3>
make_combined_force(F1 f1, F2 f2, F3 f3){
  return CombinedForce<CombinedForce<F1, F2>, F3> (CombinedForce<F1, F2>(f1,f2), f3);
}


/** Null Constraint */
struct NullConstraint {
  /** Set Null Constraint.
  * @param[in] g Valid graph.
  * @param[in] t Valid time.
  * @return void
  */
  void operator()(GraphType& g, double t) {
    (void) t; // silence compiler warnings
    (void) g;
  }
};

/** Point(0,0,0) and Point(1,0,0) cannot be moved. */
struct FixedConstraint {
  /** Set Fixed Constraint.
  * @param[in] g Valid graph.
  * @param[in] t Valid time.
  * @post The velocity of Point(0,0,0) and Point(1,0,0) are 0.
  */
  void operator()(GraphType& g, double t) {
    (void) t;     // silence compiler warnings
    for (GraphType::NodeIterator it = g.node_begin(); it != g.node_end(); ++it) {
      Node n = (*it);
      if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
        n.value().vel =  Point(0,0,0);
      }
    }
  }
};


// Constraints

/** Horizontal Plane constraint, representing a solid plane parallel to the x-y plane.
*/
struct HPlane {
  double z_; //coordinate for the horizontal plane

  /** HPlane Constructor.
  * @param[in] z_constraint Sets the z-coordinate to define the horizontal plane.
  */
  HPlane(double z) : z_(z) {}

  /** Horizontal Constraint Setter
  * @param[in] g Valid graph.
  * @param[in] t Valid time.
  * @post Any node violating this constraint has:
  *       velocity in @a z_ direction set to 0
  *       position set to the nearest point to the horizontal plane defined by @a z_.
  */
  void operator()(GraphType& g, double t) {
    (void) t;     // silence compiler warnings
    for (GraphType::NodeIterator it = g.node_begin(); it != g.node_end(); ++it) {
      Node n = (*it);
      if (n.position().z < z_) {
        n.position().z = z_;
        n.value().vel.z = 0;
      }
    }
  }
};


/** Sphere constraint, representing solid sphere defined by a center and a r.
*/
struct Sphere {
  double r_;
  Point center_;

  /** Sphere Constructor.
  * @param[in] r Sets the r of the sphere constraint.
  * @param[in] center Sets the center of the sphere constraint.
  */
  Sphere(double r, Point center) :
  r_(r), center_(center) {}

  /** Sphere Constraint Setter
  * @param[in] g Valid graph.
  * @param[in] t Valid time.
  * @post Any node that moves into the sphere has:
  *       its position set to the nearest point on the sphere.
  *       The component of the node's velocity normal to the sphere's surface
  *       is set to v - (v*R)*R,
  *       where R = (x - center) / |x - center|
  */
  void operator()(GraphType& g, double t) {
    (void) t;     // silence compiler warnings
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      Node n = (*it);
      double dist = norm(n.position() - center_);

      // Reset nodes' position to be on the sphere's surface
      if (dist < r_) {
        n.position() = (n.position() - center_) * (r_ / dist) + center_;
        Point R = (n.position() - center_)/dist;
        n.value().vel -= dot(n.value().vel, R) * R;
      }
    }
  }
};


/** Functor to remove nodes that touch the sphere's surface.
*/
struct SphereRemove {
  double r_;
  Point center_;
  /** SphereRemove Constructor.
  * @param[in] r Sets the r of the remove sphere area.
  * @param[in] center Sets the center of the remove sphere area.
  */
  SphereRemove(double r, Point center) :
  r_(r), center_(center) {}

  /** Remove Sphere Setter.
  * @param[in] g Valid graph.
  * @param[in] t Valid time.
  * @post Nodes that move through the sphere are removed from the graph.
  */
  void operator()(GraphType& g, double t) {
    (void) t;     // silence compiler warnings
    for (GraphType::NodeIterator it = g.node_begin(); it != g.node_end(); ++it) {
      Node n = (*it);
      double dist = norm(n.position() - center_);
      // Remove violating nodes
      if (dist < r_){
        g.remove_node(n);
      }
    }
  }
};


/** Functor to combine constraints.
* @param[in] c1 and c2 are valid constraints.
*/
template <typename C1, typename C2>
struct CombinedConstraint {
  C1 c1_;
  C2 c2_;

  CombinedConstraint(C1 c1, C2 c2) : c1_(c1), c2_(c2) {}

  void operator() (GraphType& g, double t){
    c1_(g,t);
    c2_(g,t);
  }
};


/** Return a CombinedConstraint object.
* @param[in] c1, c2 are valid constraints.
* @pre Valid constraints that take as input a graph g and a time t.
* @return A CombinedConstraint object combining the given constraints.
*/
template <typename C1, typename C2>
CombinedConstraint<C1, C2> make_combined_constraint(C1 c1, C2 c2){
  return CombinedConstraint<C1, C2> ({c1, c2});
}


/** Return a CombinedConstraint object.
* @param[in] c1, c2 are valid constraints.
* @pre Valid constraints that take as input a graph g and a time t.
* @return A CombinedConstraint object combining the given constraints.
*/
template <typename C1, typename C2, typename C3>
CombinedConstraint<CombinedConstraint<C1, C2>, C3>
make_combined_constraint(C1 c1, C2 c2, C3 c3){
  return CombinedConstraint<CombinedConstraint<C1, C2>, C3> ({{c1, c2}, c3});
}



// =============================================================================
int main(int argc, char** argv) {
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct an empty graph
  GraphType graph;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  std::vector<typename GraphType::node_type> nodes;
  while (CME212::getline_parsed(nodes_file, p))
  nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CME212::getline_parsed(tets_file, t)) {
    graph.add_edge(nodes[t[0]], nodes[t[1]]);
    graph.add_edge(nodes[t[0]], nodes[t[2]]);
    #if 1
    // Diagonal edges
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
    #endif
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }



  // Set initial velocity and mass
  for (GraphType::NodeIterator it = graph.node_begin(); it != graph.node_end(); ++it) {
    Node n = *it;
    n.value().vel = Point(0,0,0);       // Initial velocity == 0
    n.value().mass = 1.0 / graph.num_nodes(); // graph has total mass == 1, constant density
  }

  // Set rest length for all of the Edges to their initial length
  for (GraphType::EdgeIterator ei = graph.edge_begin(); ei != graph.edge_end(); ++ei ) {
    (*ei).value().L = (*ei).length();
  }

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the SDLViewer
  CME212::SDLViewer viewer;
  auto node_map = viewer.empty_node_map(graph);
  viewer.launch();

  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  viewer.center_view();

  // Begin the mass-spring simulation
  double dt = 0.001;
  double t_start = 0;
  double t_end = 5.0;

  for (double t = t_start; t < t_end; t += dt) {
    // P1 ---------------------------------------------------------------------
    // symp_euler_step(graph, t, dt, Problem1Force());

    // P3 ----------------------------------------------------------------------
    //symp_euler_step(graph, t, dt, cf);

    // Create individual forces
    GravityForce g(grav);
    MassSpringForce msf(100);
    DampingForce d(1.0 / graph.num_nodes());

    // Combine the individual forces
    auto cf = make_combined_force(g, msf, d);

    // P4 ----------------------------------------------------------------------
    // Create individual constraints
    HPlane hp(-0.75);
    Sphere sp(0.15, Point(0.5,0.5,-0.5));
    SphereRemove sr(0.15, Point(0.5,0.5,-0.5));

    // Combined individual constraints
    // P4.1
    // auto c = make_combined_constraint(hp, FixedConstraint());
    // P4.2
    // auto c = make_combined_constraint(sp, FixedConstraint());
    // P4.3
    auto c = make_combined_constraint(sr, FixedConstraint());
    // Mixed constraints
    // auto c = make_combined_constraint(hp, sr, FixedConstraint());
    // auto c = make_combined_constraint(hp, sp, FixedConstraint());

    symp_euler_step(graph, t, dt, cf, c);

    viewer.clear();
    node_map.clear();
    viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
    viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

    // Update viewer with nodes' new positions
    viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
    viewer.set_label(t);

    // These lines slow down the animation for small graphs, like grid0_*.
    // Feel free to remove them or tweak the constants.
    if (graph.size() < 100)
    CME212::sleep(0.0001);
  }

  return 0;
}
