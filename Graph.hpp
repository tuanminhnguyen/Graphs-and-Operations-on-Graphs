#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <iostream>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph {
private:
  struct node_elem;

public:

  // **************************************************************************
  // PUBLIC TYPE DEFINITIONS
  // **************************************************************************

  /** Type of this graph. */
  typedef Graph graph_type;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  typedef Node node_type;
  typedef V node_value_type;
  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  typedef Edge edge_type;
  typedef E edge_value_type;

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  typedef NodeIterator node_iterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  typedef EdgeIterator edge_iterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  typedef IncidentIterator incident_iterator;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  typedef unsigned size_type;

  // **************************************************************************
  // CONSTRUCTORS AND DESTRUCTOR
  // **************************************************************************

  /** Construct an empty graph. */
  Graph() : nodes_(), edges_adj_(), edges_vals_(), num_edges_(0), i2u_() {}
  /** Default destructor */
  ~Graph() = default;

  /** Return the size of the graph.
   *  @return The number of existing nodes in this graph, indicated by the
   *          size of i2u_.
   *  Complexity: O(1).
   */
  size_type size() const {
    return i2u_.size();
  }


  /** Remove all nodes and edges from this graph.
   *  @post num_nodes() == 0 && num_edges() == 0
   *
   *  Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_adj_.clear();
    edges_vals_.clear();
    i2u_.clear();
    num_edges_ = 0;
  }


  // **************************************************************************
  // NODES
  // **************************************************************************

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   * ========================================================================*/
  class Node : private totally_ordered<Node> {
  public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Graph::node_type x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */
    Node() : graph_(nullptr), uid_(0) {}


    /** Return this node's position.
     *  @return This node's Point.
     *  Complexity: O(1).
     */
    const Point& position() const {
      return graph_->nodes_[uid_].point_;
    }


    /** Return or modify this node's position.
     * @return The node's Point.
     *
     * Complexity O(1).
     */
    Point& position() {
      return graph_->nodes_[uid_].point_;
    }


    /** Return this node's index, a number in the range [0, graph_size).
     * @return The node's index i, s.t. 0 <= i < num_nodes()
     * Complexity O(1)
     * This index corresponds to the index in i2u_ (index, not value, which is uid).
     * */
    size_type index() const {
      return graph_->nodes_[uid_].uid_;
    }


    /** Test whether this node and @a n are equal.
     *  @param[in] @a n A node to be compared with this node.
     *  @return True if node and this node have the same graph pointer and the same index.
     *  Complexity O(1).
     */
    bool operator==(const Node& n) const {
      return graph_ == n.graph_ && uid_ == n.uid_;
    }


    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     * @param[in] @a n A node to be compared with this node.
     * return True if this node's graph pointer has address less than that of n, or
     *             if both nodes have the same graph pointer and this node's index
     *             is less than that of n.
     * Complexity O(1).
     */
    bool operator<(const Node& n) const {
      return (graph_ < n.graph_ || (graph_ == n.graph_ && uid_ < n.uid_));
    }


    /** Modify this node's value.
     *  @return A reference to this node's node_value_type new value.
     *  Complexity O(1).
     */
    node_value_type& value() {
      return const_cast<node_value_type&>(static_cast<const Node*>(this)->value());
    }


    /** Get the value stored at this node.
     *  @return The value stored at the current node.
     *  @post This node's value has not been modified.
     *  Complexity O(1).
     */
    const node_value_type& value() const{
      return graph_->nodes_[uid_].val_;
    }


    /** Get the number of edges incident on this node.
     *  @return the number of edges incident on this node.
     *  It must be the case that 0 <= degree() < num_nodes().
     *  Complexity O(1).
     */
    size_type degree() const {
      return graph_->edges_adj_[uid_].size();
    }


    /** Set an incident_iterator to the first of this node's incident edges.
     *  @return incident_iterator pointing to the first of this node's incident edges
     *          (as stored in edges_adj_).
     *  Complexity O(1).
     */
    incident_iterator edge_begin() const {
      return incident_iterator(graph_, uid_, 0);
    }


    /** Set an incident iterator to one past the last of this node's incident edges.
     *  @return incident_iterator pointing to one past the last of this node's incident edges.
     *          This pointer is invalid and must not be dereferenced or incremented.
     *  Complexity O(1).
     */
    incident_iterator edge_end() const {
      return incident_iterator(graph_, uid_, degree());
    }


  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    /** graph_ pointer to the graph container.
     *  uid_ this node's uid in i2u_, which is this node's index in nodes_.
     *  Representation Invariants:
     *  graph_ != nullptr
     *  0 <= uid_ < num_nodes_total()
     */
    Graph* graph_;
    size_type uid_;


    /** Construct a valid node.
     *  @param[in] graph The parent graph of this node.
     *  @param[in] uid The index of this node in nodes_, also its uid in i2u_.
     *  Complexity O(1).
     */
    Node(const Graph* graph, size_type uid) :
      graph_(const_cast<Graph*>(graph)), uid_(uid) {}
  };


  /** Return the number of nodes in the graph.
   * @return size_type The number of nodes in the graph, a nonnegative number.
   * Complexity: O(1).
   */
  size_type num_nodes() const {
    return i2u_.size();
  }


  /** Return the number of nodes in nodes_.
   * @return size_type non-negative number at least num_nodes(), but
   *         may be larger if at least one node has been removed.
   * Complexity: O(1).
   */
  size_type num_nodes_total() const {
    return nodes_.size();
  }


  /** Add a node to the graph, returning the added node.
   * @param[in] pos The new node's position.
   * @param[in] val value to be stored at this node. Default value
   *            is provided by the default constructor of node_value_type.
   * @return A copy of the node that was added to the graph.
   * @post new num_nodes() == old num_nodes() + 1.
   * @post new node.index() == old num_nodes().
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& pos,
		const node_value_type& val = node_value_type()) {
    nodes_.emplace_back(pos, val, i2u_.size());
    i2u_.push_back( nodes_.size() - 1);

    edges_adj_.resize( nodes_.size() );
    edges_vals_.resize( nodes_.size() );

    return Node(this, nodes_.size() - 1);
  }


  /** Determine if this Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (this == n.graph_ && n.uid_ < num_nodes_total());
  }


  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes().
   * @post result_node.index() == i.
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, i2u_[i]);
  }


  /** Remove the node and its associated edges from the graph.
   *  param[in] n Valid node of this Graph.
   *  @pre g.has_node(@a n).
   *  @pre g.node(i).index() == i.
   *
   *  @post !g.has_node(@a n).
   *  @post new num_nodes() = old num_nodes() - 1.
   *  @post g.node(i).index() == i.
   *  @post g.node(@a n.index()) == @a n.
   *  @post Edges e with e.node1() == @a n OR e.node2() == @a n are invalid and
   *        have been removed from the Graph.
   *  @post new num_edges() = old num_edges() - old @a n.degree().
   *  @post Invalidates existing NodeIterators.
   *  @post Invalidates existing EdgeIterators.
   *  @post Invalidates existing IncidentIterators currently iterating through
   *        incident edges on @a n.
   *  @return size_type position of where the node was in i2u_ when it
   *        was removed from i2u_.
   *  Complexity: O(n.degree()^2)).
   */
  size_type remove_node(const Node& n) {
    unsigned j = 0;
    // while (i2u_[j] != n.uid_ && j < i2u_.size()) {
    //   ++j;
    // }
    // ++j;
    // The while loop will not work if there is only 1 node in the graph.

    for(unsigned i = 0; i < i2u_.size(); ++i) {
      if(i2u_[i] == n.uid_) {
        j = i;
	break;
      }
    }

    // Rmove this node from i2u_
    i2u_.erase(i2u_.begin() + j);

    // Remove the incident edges, by calling remove_edge we only
    // need to do it for one of (a,b) and (b,a), for example.
    while( !edges_adj_[n.uid_].empty() ) {
      remove_edge(n, Node(n.graph_, edges_adj_[n.uid_][0]));
    }

    // Shift the nodes previously behind n forward (towards beginning of nodes_)
    for(unsigned i = j; i < i2u_.size(); ++i) {
      nodes_.at(i2u_[i]).uid_ = i;
    }

    return j;
  }


  /** Remove the node and its associated edges from the graph.
   *  param[in] n Valid node of this Graph.
   *  @pre g.has_node(@a n).
   *  @pre g.node(i).index() == i.
   *
   *  @post !g.has_node(@a n).
   *  @post new num_nodes() = old num_nodes() - 1.
   *  @post g.node(i).index() == i.
   *  @post g.node(@a n.index()) == @a n.
   *  @post Edges e with e.node1() == @a n OR e.node2() == @a n are invalid and
   *        have been removed from the Graph.
   *  @post new num_edges() = old num_edges() - old @a n.degree().
   *  @post Invalidates existing NodeIterators.
   *  @post Invalidates existing EdgeIterators.
   *  @post Invalidates existing IncidentIterators currently iterating through
   *        incident edges on @a n.
   * @return NodeIterator pointing to the position of where the node was when
   *       it was removed from i2u_.
   *  If the removed node was the last in i2u_, the returned result must be
   *        new i2u_.size().
   *  Complexity: O(n.degree()^2).
   */
  node_iterator remove_node(node_iterator n_it) {
    return NodeIterator(this, remove_node(*n_it));
  }


  // **************************************************************************
  // EDGES
  // **************************************************************************

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   * ========================================================================*/
  class Edge : private totally_ordered<Edge> {
  public:
    /** Construct an invalid Edge. */
    Edge() : graph_(nullptr), n1_uid_(0), n2_uid_(0) {}

    /** Return the first node of this Edge.
     *  @return This edge's first node.
     *  Complexity O(1).
     */
    Node node1() const {
      return Node(graph_, n1_uid_);
    }


    /** Return the other node of this Edge.
     *  @return This edge's second node.
     *  Complexity O(1).
     */
    Node node2() const {
      return Node(graph_, n2_uid_);
    }


    /** Test whether this edge and @a e are equal.
     *  @pre @a e is a valid edge in the graph.
     *  @return True if both edges' graph pointers point to the same graph, and
     *          the edges' nodes are the same.
     *  Complexity O(1).
     */
    bool operator==(const Edge& e) const {
      return (graph_ == e.graph_ &&
	      ((n1_uid_ == e.n1_uid_ && n2_uid_ == e.n2_uid_) ||
	       (n1_uid_ == e.n2_uid_ && n2_uid_ == e.n1_uid_)));
    }


    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     * @param[in] e A valid edge.
     * @return True if any of this edge's graph pointer address, n1_uid_,
     *         n2_uid_, is less than the corresponding datum in e.
     * Complexity O(1).
     */
    bool operator<(const Edge& e) const {
      return (std::tie(graph_, n1_uid_, n2_uid_) <
	      std::tie(e.graph_, e.n1_uid_, e.n2_uid_));
    }


    /** Return the length of this Edge.
     * @pre End nodes of this Edge have valid positions.
     * @return double Length as measured by the Euclidean distance between the
     *         positions of the end nodes.
     *
     * Complexity: O(1).
     */
    double length() const {
      return norm(node1().position() - node2().position());
    }


    /** Get or modify the value stored at this Edge.
     *  @return Reference to the value stored at the current Edge.
     *  Complexity O(1).
     */
    edge_value_type& value() {
      Edge e;
      if(node1().uid_ < node2().uid_) {
	e = Edge(graph_, node1().uid_, node2().uid_);
      }
      else{ // node1().uid_ > node2().uid_ because node1() != node2()
	e = Edge(graph_, node2().uid_, node1().uid_);
      }


      IncidentIterator it = e.node1().edge_begin();
      IncidentIterator ee = e.node1().edge_end();
      while(((*it).n2_uid_ != e.n2_uid_) && (it != ee)) ++it;

      return graph_->edges_vals_[e.n1_uid_][it.n2_idx_];
    }


    /** Get the value stored at this Edge.
     *  @return Reference to the value stored at the current Edge.
     *  @post This Edge's value has not been modified.
     *  Complexity O(1).
     */
    const edge_value_type& value() const {
      Edge e;
      if(n1_uid_ < n2_uid_) {
	e = edge(n1_uid_, n2_uid_);
      }
      else{ // n1_uid_ > uid2, since n1_uid_ != n2_uid_
	e = edge(n2_uid_, n1_uid_);
      }

      IncidentIterator it = e.node1().edge_begin();
      IncidentIterator ee = e.node1().edge_end();
      while( ((*it).n2_uid_ != e.n2_uid_) && (it != ee) ) {
      	++it;
      }

      return graph_->edges_vals_[e.n1_uid_][it.n2_idx_];
    }


  private:
    friend class Graph;

    /** Representation Invariant:
     *  graph_ != nullptr
     *  0 <= n1_uid_ < num_nodes_total();
     *  0 <= n2_uid_ < num_nodes_total();
     */
    Graph* graph_;
    size_type n1_uid_;
    size_type n2_uid_;

    /** Construct a valid Edge object.
     *  @param[in] graph pointer to the Edge's parent graph.
     *  @param[in] n1_uid index of this node in nodes_
     *  @param[in] n2_uid index of this node in nodes_
     */
    Edge(const Graph* graph, size_type n1_uid, size_type n2_uid) :
      graph_(const_cast<Graph*>(graph)), n1_uid_(n1_uid), n2_uid_ (n2_uid) {}
  };


  /** Return the total number of edges in the graph.
   * @return The number of edges,
   *         satisfying 0 <= num_edges() <= nchoosek(num_nodes(), 2),
   *         where the upper limit is the maximum number of edges in a simple graph.
   *
   * Complexity: O(1).
   */
  size_type num_edges() const {
    return num_edges_;
  }


  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: O(@a a.degree()).
   */
  Edge add_edge(const Node& a, const Node& b) {
    if(!has_edge(a, b)){
      // Add edges (a,b) and (b,a) to the adjacency vector
      edges_adj_[a.uid_].push_back(b.uid_);
      edges_adj_[b.uid_].push_back(a.uid_);
      // Initialize edge values to default
      edges_vals_[a.uid_].push_back(edge_value_type());
      edges_vals_[b.uid_].push_back(edge_value_type());
      ++num_edges_;
      return Edge(this, a.uid_, b.uid_);
    }

    return Edge(this, a.uid_, b.uid_);
  }


  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity O(@a a.degree()).
   */
  bool has_edge(const Node& a, const Node& b) const {
    size_type n1_idx = a.uid_;
    size_type n2_idx = b.uid_;

    // only check existence of (a,b) since if (a,b) exists then (b,a) exists
    // by the way add_edge() works
    for(unsigned i = 0; i < edges_adj_[n1_idx].size(); i++) {
      if (edges_adj_[n1_idx][i] == n2_idx) return true;
    }

    return false;
  }


  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges().
   * @return A copy of Edge object having index i.
   *
   * Complexity O(num_edges()).
   */
  Edge edge(size_type i) const {
    EdgeIterator it = edge_begin();
    for( ; i != 0; --i) {
      ++it;
    }

    return Edge(this, (*it).n1_uid_, (*it).n2_uid_);
  }


  /** Remove an Edge from this Graph.
   * @param[in] Nodes @a a and @a b.
   * @pre @a a and @a b are valid nodes of this graph.
   * @post new num_edge() = old num_edge() - 1 if has_Edge(a, b).
   * @post new num_edge() = old num_edge() if !has_Edge(a, b).
   * @post existing EdgeIterators are invalidated.
   * @post existing IncidentIterators are invalidated.
   * @return size_type, as follows:
   *      if has_edge(a,b) return n2_idx_, correspond to b's index
   *         in edges_adj_[a.uid_],
   *      if !has_edge(a,b) return @a old a.degree().
   *
   * Complexity: O(@a a.degree() + @a b.degree()).
   */
  size_type remove_edge(const Node& a, const Node& b) {
    size_type result = a.degree();
    // avoid using has_edge(a,b) to improve performance, assuming in most cases,
    // there exists an edge between a and b.
    bool hasEdgeab = false;

    // delete (a,b)
    for(IncidentIterator iit = a.edge_begin(); iit != a.edge_end(); ++iit) {
      if((*iit).n2_uid_ == b.uid_) {
    	edges_adj_[a.uid_].erase(edges_adj_[a.uid_].begin() + iit.n2_idx_);
    	edges_vals_[a.uid_].erase(edges_vals_[a.uid_].begin() + iit.n2_idx_);
    	result = iit.n2_idx_;
	hasEdgeab = true;
    	break;
      }
    }

    // delete (b,a)
    for(IncidentIterator iit = b.edge_begin(); iit != b.edge_end(); ++iit) {
      if( (*iit).n2_uid_ == a.uid_) {
    	edges_adj_[b.uid_].erase(edges_adj_[b.uid_].begin() + iit.n2_idx_);
    	edges_vals_[b.uid_].erase(edges_vals_[b.uid_].begin() + iit.n2_idx_);
    	break;
      }
    }

    if (hasEdgeab) --num_edges_;
    return result;
  }

  /** Remove an Edge from this Graph.
   * @param[in] Edge e.
   * @pre Valid edge @a e of this graph.
   * @post new num_edge() = old num_edge() - 1.
   * @post existing EdgeIterators are invalidated.
   * @post existing IncidentIterators are invalidated.
   * @return size_type  e.node2()'s index in edges_adj_[e.node1().uid_].
   *
   * Complexity: O(@a e.node1().degree() + @a e.node2().degree()).
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }


  /** Remove Edge from this Graph.
   * @param[in] EdgeIterator.
   * @pre the EdgeIterator is valid for this graph.
   * @post new num_edge() = old num_edge() - 1.
   * @post existing EdgeIterators are invalidated.
   * @post existing IncidentIterators are invalidated.
   * @return EdgeIterator that pointing to the position of the next valid
   * edge in edges_adj_. This is graph.edge_end() if the removed edge is the
   * last edge or if there was no edge in the graph.
   *
   * Complexity: O((*eit).node1().degree() + (*eit).node1().degree()).
   */
  edge_iterator remove_edge(edge_iterator eit) {
    return EdgeIterator(this, remove_edge((*eit).node1(), (*eit).node2()));
  }


  // **************************************************************************
  // NODE ITERATOR
  // **************************************************************************

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {
  public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Node value_type;
    /** Type of pointers to elements. */
    typedef Node* pointer;
    /** Type of references to elements. */
    typedef Node& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid NodeIterator. */
    NodeIterator() : graph_(nullptr), idx_(0) {}


    /** Prefix-increment a node_iterator.
     *  @return A reference to the incremented node_iterator.
     *  @pre The argument does not point past the last node.
     *  @post Increases the @a idx_ to the next node's idx_.
     *  Complexity O(1).
     */
    Node operator*() const {
      return Node(graph_, graph_->i2u_[idx_]);
    }


    /** Prefix-increment a node_iterator.
     *  @return A reference to the incremented node_iterator.
     *  @pre The argument does not point past the last node.
     *  @post Increases the @a idx_ to the next node's idx_.
     *  Complexity O(1).
     */
    node_iterator& operator++() {
      idx_++;
      return *this;
    }


    /** Compare another node_iterator with this node_iterator.
     *  @param[in] b other node_iterator to be compared.
     *  @return True if b has the same graph pointer and index as this node_iterator.
     *  @pre b.idx_ is a valid node index.
     *  @post b is not modified.
     *  Complexity O(1).
     */
    bool operator==(const NodeIterator& b) const {
      return (graph_ == b.graph_ && idx_ == b.idx_);
    }


  private:
    friend class Graph;

    /* Representation Invariants
     * g_ != nullptr
     * 0 <= idx_ < num_nodes()
     */
    Graph* graph_;
    size_type idx_; // index in i2u_

    NodeIterator(const Graph* graph, size_type idx) :
      graph_(const_cast<Graph*>(graph)), idx_(idx) {}
  };


  /** Get an instant of  node_iterator that points to the first Node.
   *  @return a node_iterator that points to the first Node.
   */
  node_iterator node_begin() const {
    return node_iterator(this, 0);
  }


  /** Get an instant of  node_iterator that points to one past the last Node.
   *  @return a node_iterator that points to one past the last Node.
   */
  node_iterator node_end() const {
    return node_iterator(this, num_nodes());
  }


  // **************************************************************************
  // EDGE ITERATOR
  // **************************************************************************

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
  public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() : graph_(nullptr), n1_uid_(0), n2_idx_(0) {}


    /** Dereference an EdgeIterator.
     *  @pre EdgeIterator != edge_end().
     *  @post the Edge object pointed to by this EdgeIterator is not modified.
     *  @return the Edge object pointed to by the EdgeIterator object.
     *  Complexity O(1).
     */
    Edge operator*() const {
      return Edge(graph_, n1_uid_, graph_->edges_adj_[n1_uid_][n2_idx_]);
    }


    /** Prefix-increment an EdgeIterator object.
     *  @pre EdgeIterator != edge_end().
     *  @post The old and new Edges pointed to by this Edge Iterator are not
     *        modified.
     *  @return a reference to the EdgeIterator that points to the next Edge.
     */
    EdgeIterator& operator++() {
      //increment_n_index1();
      ++n2_idx_;
      fix();
      return *this;
    }


    /** Compare another EdgeIterator with this EdgeIterator.
     *  @param b a reference to the other EdgeIterator.
     *  @pre b refers to a valid EdgeIterator.
     *  @post no EdgeIterator is modified.
     *  @return true if b is identical to this EdgeIterator. False otherwise.
     */
    bool operator==(const EdgeIterator& b) const {
      return graph_ == b.graph_ &&
	n1_uid_ == b.n1_uid_ &&
	n2_idx_ == b.n2_idx_;
    }


  private:
    friend class Graph;

    /** Representation Invariant:
     *  graph_ != nullptr
     *  0 <= n1_uid_ < graph_->edges_adj_.size()
     *  0 <= n2_idx_ < graph_->edges_adj_[n1_uid_].size()
     */
    Graph* graph_;     // pointer to the parent graph
    size_type n1_uid_; // index of the first node, as in edges_adj_
    size_type n2_idx_; // index of the second node, as in edges_adj_[n1_uid_]


    /** Fix an invalid edge_iterator (i.e. invalid if node1().index() > node2().index().
     *   Move the edge_iterator forward until  node1().index() < node2().index().
     */
    void fix() {
      while (n1_uid_ < graph_->edges_adj_.size()) {
        while (n2_idx_ < graph_->edges_adj_[n1_uid_].size()) {
       	  if (graph_->edges_adj_[n1_uid_][n2_idx_] > n1_uid_){
      	    return;
          }
      	  ++n2_idx_;
        }
        n2_idx_ = 0;
        ++n1_uid_;
      }
    }


    /** Construct a valid EdgeIterator.
     *  @param[in] graph this EdgeIterator's pointer to the parent graph.
     *  @param[in] n1_uid Index of the first node of the edge pointed to by this EdgeIterator.
     *  @param[in] n2_idx Index of the second node of the edge pointed to by this EdgeIterator.
     */
    EdgeIterator(const Graph* graph, size_type n1_uid, size_type n2_idx) :
      graph_(const_cast<Graph*>(graph)), n1_uid_(n1_uid), n2_idx_(n2_idx) {
      fix();
    }
  };


  /** Set an edge_iterator to the first edge of this graph's edges.
   *  @return edge_iterator pointing to the first edge of this graph's edges,
   *          i.e. graph->edges_adj_[0][0].
   *  edge_iterator == edge_end() if num_node() == 0.
   *  Complexity O(1).
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0, 0);
  }


  /** Set an edge_iterator to one past the last edge of this graph's edges.
   *  @return edge_iterator pointing to one past the last edge of this graph's edges,
   *          here we choose the representation graph_->edges_adj_[num_nodes()][0].
   *          This edge_iterator is invalid and must not be dereferenced.
   *  edge_iterator == edge_end() if num_node() == 0.
   *  Complexity O(1).
   */
  edge_iterator edge_end() const {
    return edge_iterator(this, edges_adj_.size(), 0);
  }



  // **************************************************************************
  // INCIDENT ITERATOR
  // **************************************************************************

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
  public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() : graph_(nullptr), n1_uid_(0), n2_idx_(0) {}


    /** Dereference an IncidentIterator.
     *  @pre IncidentIterator != node.edge_end().
     *  @return The Edge pointed to by the IncidentIterator.
     *  Complexity O(1).
     */
    Edge operator*() const {
      return Edge(graph_, n1_uid_, graph_->edges_adj_[n1_uid_][n2_idx_]);
    }


    /** Increment an IncidentIterator to the next position.
     *  @pre IncidentIterator != node.edge_end().
     *  @post IncidentIterator points to the next incident edge (possibly invalid edge).
     *  @return A reference of the updated IncidentIterator.
     *  Complexity O(1).
     */
    IncidentIterator& operator++() {
      ++n2_idx_;
      return *this;
    }


    /** Compare another IncidentIterator with this IncidentIterator.
     *  @param[in] b the other IncidentIterator
     *  @return True if b and this IncidentIterator have the same graph pointer
     *  and end-nodes' indices.
     *  Complexity O(1).
     */
    bool operator==(const IncidentIterator& b) const {
      return (graph_ == b.graph_ &&
	      n1_uid_ == b.n1_uid_ &&
	      n2_idx_ == b.n2_idx_);
    }


  private:
    friend class Graph;

    /** Representation Invariant:
     * graph_ != nullptr
     * 0 <= n1_uid_ < num_nodes()
     * 0 <= n2_idx_ < graph_->edges_adj_[n1_uid_].size()
     */
    Graph* graph_;
    size_type n1_uid_;
    size_type n2_idx_;

    /** Construct a valid IncidentIterator.
     *  @param[in] graph Pointer to the parent graph container.
     *  @param[in] n1_uid Index of the first node, in edges_adj_
     *  @param[in] n2_idx Index of the second node, in edges_adj_[n1_uid]
     */
    IncidentIterator(const Graph* graph, size_type n1_uid, size_type n2_idx) :
        graph_(const_cast<Graph*>(graph)), n1_uid_(n1_uid), n2_idx_(n2_idx) {}
  };


// =============================================================================
// =============================================================================
private:
  struct node_elem {
    Point point_;
    node_value_type val_;
    size_type uid_;
    node_elem(const Point& point, const node_value_type& val, size_type uid) :
      point_(point), val_(val), uid_(uid) {}
  };

  // Node vector
  std::vector<node_elem > nodes_;
  // Adjacency vector; outer indices correspond to node1's uid in i2u_,
  // stored values correspond to node 2's uid in i2u_.
  std::vector<std::vector<size_type>> edges_adj_;
  // Edge value vector
  std::vector<std::vector<edge_value_type>> edges_vals_;
  // Number of edges in the graph
  size_type num_edges_;
  // The value stored at index i is the i-th node's index in nodes_.
  std::vector<size_type> i2u_;
};

#endif
