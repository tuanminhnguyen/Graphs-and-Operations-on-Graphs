/**
 * @file poisson.cpp
 * Test script for treating the Graph as a MTL Matrix
 * and solving a Poisson equation.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles.
 * Second file: Eges (one per line) defined by 2 indices into the point list
 *              of the first file.
 *
 * Launches an SDLViewer to visualize the solution.
 */

#include <fstream>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
#include "CME212/SDLViewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"
#include "CME212/BoundingBox.hpp"
#include "math.h"
#include "Graph.hpp"

// Define a GraphSymmetricMatrix that maps
// your Graph concept to MTL's Matrix concept. This shouldn't need to copy or
// modify Graph at all!
typedef Graph<char,char> GraphType;  //<  DUMMY Placeholder
typedef mtl::dense_vector<double> VectorType;
void tagboundary_nodes();

/** Represent an abstract matrix using an underlying graph.
 */
struct GraphSymmetricMatrix {
  GraphType& g;
  size_t n;
  GraphSymmetricMatrix(GraphType& graph) : g(graph), n(g.size()) {
    tagboundary_nodes();
  };

  /** Helper function to perform multiplication. Allows for delayed
   * evaluation of results.
   * Assign::apply(a, b) resolves to an assignment operation such as
   * a += b, a-= b, or a = b.
   * @pre @a size(v) == size(w)
   */
  template <typename VectorIn, typename VectorOut, typename Assign>
  void mult(const VectorIn& v, VectorOut& w, Assign) const {
    assert(size(w) == size(v));

    for(size_t i = 0; i < g.size(); ++i) {
      double total = 0;
      for(size_t j = 0; j < g.size(); ++j) {
	total += A(i, j) * v[j];
      }
      w[i] = Assign::apply( w[i], total );
    }
  }

  /** Implicit Poisson matrix.
   * @param[in] i, j row and column indices of entries in a matrix.
   * @return 1 if @a i == @a j && node @a i is a boundary node
   *         0 if @a i != @a j && (node @a i is a boundary node ||
   *                                    @a j is a boundary node)
   *         L(@a i, @a j), otherwise.
   *
   * Complexity: O(1)
   */
  double A(size_t i, size_t j) const {
    if(i == j && g.node(i).value())  return 1;
    else if(i != j && (g.node(i).value() || g.node(j).value()))  return 0;
    else  return L(i, j);
  }

  /** Laplacian operator.
   * @param[in] i, j row and column indices of an entry in a matrix.
   * @return -node(i).degree() if @a i == @a j.
   *         1 if has_edge(node(i), node(j)).
   *         0, otherwise.
   *
   * Complexity: O(1)
   */
  double L(size_t i, size_t j) const {
    if (i == j)  return -1.0 * g.node(i).degree();
    else if (g.has_edge(g.node(i), g.node(j)))  return 1;
    else  return 0;
  }

  /** Matrix-vector multiplication forwards to MTL's lazy mat_cvec_multiplier operator
   */
  template<typename Vector>
  mtl::vec::mat_cvec_multiplier<GraphSymmetricMatrix, Vector>
  operator*(const Vector& v) const {
    return mtl::vec::mat_cvec_multiplier
      <GraphSymmetricMatrix, Vector>(*this, v);
  }


  /** Tag nodes in a graph as either on boundary or not on boundary.
   *  @param[in] g reference to a graph.
   *  @post      for any node n in g,
   *               if norm_inf(n.position()) == 1 ||
   *                  norm_inf(n.position() - (+/-0.6, +/-0.6, 0)) < 0.2 ||
   *                  n.position() is in the bounding box defined by (-0.6, -0.2, -1), (0.6, 0.2, 1), then
   *                  n.value() == 1.
   *               n.value() == 0 otherwise.
   *
   *  Complexity: O(1).
   */
  void tagboundary_nodes() {
    GraphType::NodeIterator ni = g.node_begin();
    GraphType::NodeIterator ne = g.node_end();

    for (; ni != ne; ++ni) {
      Point pos = (*ni).position();
      if (norm_inf(pos) == 1) {
	(*ni).value() = 1;
      }

      else if ((norm_inf(pos - Point(0.6, 0.6, 0)) < 0.2)  ||
	       (norm_inf(pos - Point(-0.6, 0.6, 0)) < 0.2) ||
	       (norm_inf(pos - Point(-0.6, -0.6, 0)) < 0.2) ||
	       (norm_inf(pos - Point(0.6, -0.6, 0)) < 0.2)) {
	(*ni).value() = 1;
      }

      else if (Box3D(Point(-0.6, -0.2, -1), Point(0.6, 0.2, 1)).contains(pos)) {
	(*ni).value() = 1;
      }

      else {
	(*ni).value() = 0;
      }
    }
  }

};

/** The number of elements in the matrix. */
inline std::size_t size(const GraphSymmetricMatrix& A) {
  return A.n * A.n;
}

/** The number of rows in the matrix. */
inline std::size_t num_rows(const GraphSymmetricMatrix& A) {
  return A.n;
}

/** The number of columns in the matrix. */
inline std::size_t num_cols(const GraphSymmetricMatrix& A) {
  return A.n;
}

/** Traits that MTL uses to determine properties of our GraphSymmetricMatrix. */
namespace mtl {
  namespace ashape {

    /** Define GraphSymmetricMatrix to be a non-scalar type. */
    template<>
    struct ashape_aux<GraphSymmetricMatrix> {
      typedef nonscal type;
    };
  }

  /** GraphSymmetricMatrix implements the Collection concept
   *  with value_type and size_type. */
  template<>
  struct Collection<GraphSymmetricMatrix> {
    typedef char value_type;
    typedef unsigned size_type;
  };
}

/** Remove all the nodes in graph @a g whose position is within Box3D @a bb.
 * @post For all i, 0 <= i < @a g.num_nodes(),
 *        not bb.contains(g.node(i).position())
 */
void remove_box(GraphType& g, const Box3D& bb) {
  // HW3: YOUR CODE HERE
  GraphType::NodeIterator ni = g.node_begin();
  GraphType::NodeIterator ne = g.node_end();

  // for ( ; ni != ne; ++ni) {
  //   if (bb.contains((*ni).position())) g.remove_node(ni);
  // } This is incorrect. It won't keep up with invalidated iterators.

  while (ni != ne) {
    if (bb.contains((*ni).position())) ni = g.remove_node(ni);
    else ++ni;
  }
  return;
}

/** Function to fill in values for the right hand side b */
double func_f(Point& pos) {
  return 5 * cos(norm_1(pos));
}

/** Function to fill in values for the right hand side b */
double func_g(Point& pos) {

  if (norm_inf(pos) == 1) {
    return 0;
  }

  else if ((norm_inf(pos - Point(0.6, 0.6, 0)) < 0.2) ||
	   (norm_inf(pos - Point(-0.6, 0.6, 0)) < 0.2) ||
	   (norm_inf(pos - Point(-0.6, -0.6, 0)) < 0.2) ||
	   (norm_inf(pos - Point(0.6, -0.6, 0)) < 0.2)) {
    return -0.2;
  }

  else if (Box3D(Point(-0.6, -0.2, -1), Point(0.6, 0.2, 1)).contains(pos)) {
    return 1;
  }
}


/** NodeColor functor to color the nodes in the solution */
struct NodeColor {

  CME212::Color operator()(const GraphType::Node& n) const {
    if (cos(norm(n.position())) > 0.9) {
      return CME212::Color::make_hsv(0.8,1,1);
    }
    else if (cos(norm(n.position())) > 0.8) {
      return CME212::Color::make_hsv(0.7,0.5,1);
    }
    else if (cos(norm(n.position())) > 0.7) {
      return CME212::Color::make_hsv(0.3, 0.8,1);
    }
    else if (cos(norm(n.position())) > 0.6) {
      return CME212::Color::make_hsv(0.5,0.8,0.8);
    }
    else if (cos(norm(n.position())) > 0.5) {
      return CME212::Color::make_hsv(0.5,0.5,0.4);
    }
    else if (cos(norm(n.position())) > 0.4) {
      return CME212::Color::make_hsv(0.8,0.6,0.3);
    }
    else if (cos(norm(n.position())) > 0.3) {
      return CME212::Color::make_hsv(0.7, 0.5,0.7);
    }
    else if (cos(norm(n.position())) > 0.2) {
      return CME212::Color::make_hsv(0.9,0.9,0.9);
    }
    else if (cos(norm(n.position())) > 0.1) {
      return CME212::Color::make_hsv(0.9,0.5,0.9);
    }
    else if (cos(norm(n.position())) > 0.0) {
      return CME212::Color::make_hsv(0.7, 0.5,0.9);
    }
    else if (sin(norm(n.position())) > -0.1) {
      return CME212::Color::make_hsv(1,1,0.5);
    }
    else {
      return CME212::Color::make_hsv(0.9,0.5,0.8);
    }
  }
};


/** NodePosition function resets the z-coordinate of the nodes to their solution values. */
template<typename Vector>
struct NodePosition {
  Vector u_;

  /** Constructor.
   * @param[in] u solution vector.
   */
  NodePosition(const Vector& u) : u_(u) {};

  /** Operator().
   * @param[in] n Valid node.
   * @return a Point with the z-coordinate set to the solution at its node.
   */
  Point operator()(const GraphType::Node& n) const {
    return Point(n.position().x, n.position().y, u_[n.index()]);
  }
};


namespace itl {
  /** visual_iteration class */
  template<typename Real, typename ColorType, typename PositionType>
  class visual_iteration : public cyclic_iteration<Real> {
    typedef cyclic_iteration<Real> super;

    CME212::SDLViewer& v_;
    GraphType& g_;
    VectorType& u_;   // solution
    ColorType& c_;    // color functor
    PositionType& p_; // position functor
    std::map<typename GraphType::node_type, unsigned> node_map_;

    /** Viewer initialization */
    void initViewer() {
      node_map_ = v_.empty_node_map(g_);
      v_.launch();
      updateViewer();
    }

    /** Update viewer using color and position functors */
    void updateViewer() {
      v_.add_nodes(g_.node_begin(), g_.node_end(), c_, NodePosition<VectorType>(u_), node_map_);
      v_.add_edges(g_.edge_begin(), g_.edge_end(), node_map_);
      v_.center_view();
    }

  public:
    /** Constructor
     *  @param[in] r, max_iter, tol, tol2, c parameters for cyclic_iteration.
     *  @param[in] v, g, u Viewer, graph and solution vector.
     */
    template<typename Vector>
    visual_iteration(const Vector& r, int max_iter, Real tol, Real tol2, int cp,
		     CME212::SDLViewer& v, GraphType& g, VectorType& u,
		     ColorType c, PositionType p) : super(r, max_iter, tol, tol2, cp), v_(v), g_(g), u_(u), c_(c), p_(p) {
      initViewer();
    }

    /** Update viewer
     *  @return True if super::finished().
     */
    bool finished() {
      updateViewer();
      return super::finished();
    }

    /** Overloaded finished()
     *  @param[in] r
     *  @return True if super::finished().
     */
    template<typename T>
    bool finished(const T& r) {
      updateViewer();
      bool ret = super::finished(r);
      return ret;
    }
  };
}


// =============================================================================
int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Define an empty Graph
  GraphType graph;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);

  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  std::vector<typename GraphType::node_type> node_vec;
  Point p;
  while (CME212::getline_parsed(nodes_file, p))
    node_vec.push_back(graph.add_node(2*p - Point(1,1,0)));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CME212::getline_parsed(tets_file, t)) {
    graph.add_edge(node_vec[t[0]], node_vec[t[1]]);
    graph.add_edge(node_vec[t[0]], node_vec[t[2]]);
    graph.add_edge(node_vec[t[1]], node_vec[t[3]]);
    graph.add_edge(node_vec[t[2]], node_vec[t[3]]);
  }

  // Get the edge length, should be the same for each edge
  auto it = graph.edge_begin();
  assert(it != graph.edge_end());
  double h = norm((*it).node1().position() - (*it).node2().position());

  // Make holes in our Graph
  remove_box(graph, Box3D(Point(-0.8+h,-0.8+h,-1), Point(-0.4-h,-0.4-h,1)));
  remove_box(graph, Box3D(Point( 0.4+h,-0.8+h,-1), Point( 0.8-h,-0.4-h,1)));
  remove_box(graph, Box3D(Point(-0.8+h, 0.4+h,-1), Point(-0.4-h, 0.8-h,1)));
  remove_box(graph, Box3D(Point( 0.4+h, 0.4+h,-1), Point( 0.8-h, 0.8-h,1)));
  remove_box(graph, Box3D(Point(-0.6+h,-0.2+h,-1), Point( 0.6-h, 0.2-h,1)));


  // HW3: YOUR CODE HERE
  // Define b using the graph, f, and g.
  // Construct the GraphSymmetricMatrix A using the graph
  // Solve Au = b using MTL.
  GraphSymmetricMatrix A(graph);

  // start with x == 0
  VectorType x(graph.size(), 0.0), b(graph.size(), 0.0), c(graph.size(), 0.0);

  // Fill in values of the right hand side b
  GraphType::NodeIterator ni = graph.node_begin();
  GraphType::NodeIterator ne = graph.node_end();

  for (int i = 0; ni != ne; ++ni, ++i) {
    GraphType::Node n = (*ni);
    if (n.value()) {
      b[i] = func_g(n.position());
    }
    else {
      b[i] = h*h*func_f(n.position());
      for (GraphType::IncidentIterator ii = n.edge_begin(); ii != n.edge_end(); ++ii) {
	GraphType::Node temp = (*ii).node2();
	if ((*ii).node2().value()) {
	  b[i] -= func_g(temp.position());
	}
      }
    }
  }

  ni = graph.node_begin();

  // Create an identity preconditioner
  itl::pc::identity<GraphSymmetricMatrix> P(A);

  // Termination criterion: r < 1.e-6
  // itl::cyclic_iteration<double> iter(b, 100, 1.e-10, 0.0, 50);

  // Launch the SDLViewer
  CME212::SDLViewer viewer;
  // viewer.launch();
  // auto node_map = viewer.empty_node_map(graph);
  // viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  // viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  NodePosition<VectorType> nodePos(x);

  // update viewer as solution converges
  itl::visual_iteration<double, NodeColor, NodePosition<VectorType>>
    viter(b, 100, 1.e-10, 0.0, 50, viewer, graph, x, NodeColor(), nodePos);

  // Solve Ax == b with preconditioner P
  cg(A, x, b, P, viter);
  return 0;
}
