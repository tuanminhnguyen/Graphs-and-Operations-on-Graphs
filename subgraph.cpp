/**
* @file subgraph.cpp
* Test script for viewing a subgraph from our Graph
*
* @brief Reads in two files specified on the command line.
* First file: 3D Points (one per line) defined by three doubles
* Second file: Tetrahedra (one per line) defined by 4 indices into the point
* list
*/

#include <fstream>
#include <iterator>

#include "CME212/SDLViewer.hpp"
#include "CME212/Util.hpp"

#include "Graph.hpp"
#include "math.h"

/** An iterator that skips over elements of another iterator based on whether
* those elements satisfy a predicate.
*
* Given an iterator range [@a first, @a last) and a predicate @a pred,
* this iterator models a filtered range such that all i with
* @a first <= i < @a last and @a pred(*i) appear in order of the original range.
*/
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

// Specify and write an interesting predicate on the nodes.
// Explain what your predicate is intended to do and test it.
// If you'd like you may create new nodes and tets files.

/**  Carve out a mohawk of the skull.
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
    (sqrt(normSq(node.position() - p2)) > dist) ;  // right ear
  }
};

int main(int argc, char** argv) {
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct a Graph
  typedef Graph<int,char> GraphType;
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
  for (unsigned j = 0; j < i; ++j)
  graph.add_edge(nodes[t[i]], nodes[t[j]]);

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the SDLViewer
  CME212::SDLViewer viewer;
  viewer.launch();

  // Use the filter_iterator to plot an induced subgraph.

  // Initializing filter_iterator objects
  //typedef SlicePredicate predFunc;
  typedef MoHawkPredicate predFunc;
  filter_iterator<predFunc, GraphType::NodeIterator> first =
  make_filtered(graph.node_begin(), graph.node_end(), predFunc());
  filter_iterator<predFunc, GraphType::NodeIterator> last =
  make_filtered(graph.node_end(), graph.node_end(), predFunc());

  // Add nodes and edges to the viewer, using the filter_iterator objects
  auto node_map = viewer.empty_node_map(graph);
  viewer.add_nodes(first, last, node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  return 0;
}
