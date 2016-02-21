/**
 * @file shortest_path.cpp
 * Test script for using our templated Graph to determine shortest paths.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <vector>
#include <fstream>

#include "CME212/SDLViewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"

#include "Graph.hpp"

#include <queue>

typedef Graph<int, char> GraphType;

/** Comparator that compares the distance from a given point p.
 */
struct MyComparator {
  Point point_;
  /** Construct a valid MyComparator.
   *  @param[in] point Valid point.
   */
  MyComparator(const Point& point) : point_(point) {
  };

  /** Determine if node1 is closer to point_ than node2.
   *  @param[in] node1 Valid node.
   *  @param[in] node2 Another valid node.
   *  @return True if the 2-norm distance between point_ and node1 is less than
   *          that between point_ and node2. False aotherwise.
   */
  template <typename NODE>
  bool operator()(const NODE& node1, const NODE& node2) const {
    // Squared-2-norm gives same result as 2-norm so no need to take sqrt
    // Use Point::normSq()
    return normSq(node1.position() - point_) <
      normSq(node2.position() - point_);
  }
};


/** Calculate shortest path lengths in @a g from the nearest node to @a point.
 * @param[in,out] g Input graph
 * @param[in] point Point to find the nearest node to.
 * @post Graph has modified node values indicating the minimum path length
 *           to the nearest node to @a point
 * @post Graph nodes that are unreachable to the nearest node to @a point have
 *           the value() -1.
 * @return The maximum path length found.
 *
 * Finds the nearest node to @a point and treats that as the root node for a
 * breadth first search.
 * This sets node's value() to the length of the shortest path to
 * the root node. The root's value() is 0. Nodes unreachable from
 * the root have value() -1.
 */
int shortest_path_lengths(GraphType& graph, const Point& point) {

  // Find the root node
  GraphType::NodeIterator root_it = std::min_element(graph.node_begin(), graph.node_end(), MyComparator(point));

  // Initialize all node values to -1
  for (GraphType::NodeIterator node_iter = graph.node_begin(); node_iter != graph.node_end(); ++node_iter) {
    (*node_iter).value() = -1;
  }

  // Initialize data structures for breath-first search
  int longest_path = 0;
  (*root_it).value() = 0;
  std::queue<GraphType::Node> nodes_queue;
  // push the root node into the queue
  nodes_queue.push(*root_it);

  // Breath-first search
  while (!nodes_queue.empty()) {
    GraphType::Node curr_node = nodes_queue.front();
    for (GraphType::IncidentIterator in_iter = curr_node.edge_begin();
	 in_iter != curr_node.edge_end(); ++in_iter) {
      if ((*in_iter).node2().value() == -1) { // this node has not been visited
        // Add one edge length to the path length
        (*in_iter).node2().value() = curr_node.value() + 1;
        // Push unvisited  node into the queue
        nodes_queue.push((*in_iter).node2());
      }
      if ((*in_iter).node2().value() )
	longest_path = (*in_iter).node2().value();
    }
    nodes_queue.pop();
  }
  return longest_path;
}

/** Create a functor to visualize in colors the shortest path lengths in the skull.
 *  @param[in] Longest_path The longest path which is output from shortest_path_lengths().
 *  @return A functor which returns a heat map, with the following color code: blue for longer path lengths and red for shorter path lengths.
 */
typedef CME212::Color color;
struct ColorFunctor {
  float longest_path_;
  ColorFunctor(const float& longest_path) : longest_path_(longest_path) {};

  CME212::Color operator()(const GraphType::Node& node) const {
    return CME212::Color::make_heat(node.value()/longest_path_);
  }
};



int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct a Graph

  GraphType graph;
  std::vector<GraphType::node_type> nodes;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  while (CME212::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CME212::getline_parsed(tets_file, t))
    for (unsigned i = 1; i < t.size(); ++i)
      for (unsigned j = 0; j < i; ++j) {
	graph.add_edge(nodes[t[i]], nodes[t[j]]);
      }

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;
  std::cout << std::endl;
  // Launch the SDLViewer
  CME212::SDLViewer viewer;
  viewer.launch();

  // Use shortest_path_lengths to set the node values to the path lengths
  int longest_path = shortest_path_lengths(graph, Point(-1, 0, 1));
  std::cout << "Longest path length = " << longest_path << std::endl;

  auto node_map = viewer.empty_node_map(graph);
  viewer.add_nodes(graph.node_begin(), graph.node_end(), ColorFunctor(longest_path),  node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  return 0;
}
