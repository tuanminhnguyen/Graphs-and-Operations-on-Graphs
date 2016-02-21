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
#include "math.h"

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

template <typename Pred, typename It>
class filter_iterator : private equality_comparable<filter_iterator<Pred,It>> {
public:
  // Get all of the iterator traits and make them our own
  typedef typename std::iterator_traits<It>::value_type        value_type;
  typedef typename std::iterator_traits<It>::pointer           pointer;
  typedef typename std::iterator_traits<It>::reference         reference;
  typedef typename std::iterator_traits<It>::difference_type   difference_type;
  typedef typename std::input_iterator_tag                     iterator_category;

  typedef filter_iterator<Pred,It> self_type;

  /** Constructor for filter_iterator.
  *  @param[in] pred Predicate to keep certain nodes if it returns True.
  *  @param[in] first Iterator pointing to the first position on which @a pred is called.
  *  @param[in] last  Iterator pointing to the end position in @a first.
  */
  filter_iterator(const Pred& pred, const It& first, const It& last)
  : pred_(pred), it_(first), end_(last) {
    fix();
  }

  /** Deference a filter_iterator.
  *  @pre filter_iterator != @a first end.
  *  @return Valid value type pointed to by the fiter_iterator.
  *
  *  Complexity: O(*first)
  */
  value_type operator*() const{
    return *it_;
  }

  /** Increment filter_iterator.
  *  @post filter_iterator points the next valid position, as determined by a predicate.
  *  @return Reference of the filter_iterator that points to the next valid position.
  *
  *  Complexity: O(++first).
  */
  self_type& operator++() {
    ++it_;
    fix();
    return *this;
  }

  /** Compare another filter_operator with this fiter_operator.
  *  @param[in] b filter_iterator in the same graph, to be compared with this filter_operator.
  *  @return True if this filter_iterator and @a b refer to the same iterator.
  *
  *  Complexity: O(==first).
  */
  bool operator==(const self_type& b) const {
    return (it_ == b.it_);
  }


private:
  /** Representation Invariant:  it_ == end_ || p(*it_).
  */
  Pred pred_;
  It it_;
  It end_;

  /** Correct the Representation Invariant.
  *  @post filter_iterator points to a valid position.
  */
  void fix() {
    while (it_ != end_ && !pred_(*it_)) ++it_;
  }
};

/** Helper function for constructing filter_iterators.
*
* Usage:
* // Construct an iterator that filters odd values out and keeps even values.
* std::vector<int> a = ...;
* auto it = make_filtered(a.begin(), a.end(), [](int k) {return k % 2 == 0;});
*/
template <typename Pred, typename Iter>
filter_iterator<Pred,Iter> make_filtered(const Iter& it, const Iter& end,
  const Pred& pred) {
    return filter_iterator<Pred,Iter>(pred, it, end);
  }


/** Test predicate
*  Carve out a mohawk of the skull.
*/
struct MoHawkPredicate {
  template <typename NODE>
  bool operator()(const NODE& node) {
    auto x = node.position().x;
    auto z = node.position().z;
    Point p1 = Point(0.3, 0, 0.15); // center for left ear ball
    Point p2 = Point(-0.3, 0, 0.15);// center right ear ball
    double dist = 0.1;
    return
    ((x < -0.09 || (x > -0.05 && x < -0.02) || (x < 0.05 && x > 0.02) || x > 0.09) || // nodes that are not removed
    (x >= -0.02 && x <= 0.02 && z < 0.15) || // middle hair line
    (x >= 0.05 && x <= 0.09 && z < 0.17) ||  // right hair line
    (x >= -0.09 && x <= -0.05 && z < 0.17)) && // left hair line
    (z < 0.1 || z > 0.12) &&                  // headband
    (sqrt(normSq(node.position() - p1)) > dist) && // left ear
    (sqrt(normSq(node.position() - p2)) > dist);  // right ear
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

  typedef MoHawkPredicate predFunc;
  filter_iterator<predFunc, GraphType::NodeIterator> first =
  make_filtered(graph.node_begin(), graph.node_end(), predFunc());
  filter_iterator<predFunc, GraphType::NodeIterator> last =
  make_filtered(graph.node_end(), graph.node_end(), predFunc());

  // auto node_map = viewer.empty_node_map(graph);
  // viewer.add_nodes(graph.node_begin(), graph.node_end(), ColorFunctor(longest_path),  node_map);
  // viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  auto node_map = viewer.empty_node_map(graph);
  viewer.add_nodes(first, last, ColorFunctor(longest_path), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  return 0;
}
