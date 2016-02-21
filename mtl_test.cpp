/**
 * @file mtl_test.cpp
 * Test script for interfacing with MTL4 and it's linear solvers.
 */


#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>


// Define a IdentityMatrix that interfaces with MTL
struct IdentityMatrix {
  int n;
  IdentityMatrix(int dim) : n(dim) {};

  /** Helper function to perform multiplication. Allows for delayed
   * evaluation of results.
   * Assign::apply(a, b) resolves to an assignment operation such as
   * a += b, a-= b, or a = b.
   * @pre @a size(v) == size(w) */
  template <typename VectorIn, typename VectorOut, typename Assign>
  void mult(const VectorIn& v, VectorOut& w, Assign) const {
    w = v;
  }

  /** Matrix-vector multiplication forwards to MTL's lazy mat_cvec_multiplier operator
   */
  template<typename Vector>
  mtl::vec::mat_cvec_multiplier<IdentityMatrix, Vector>
  operator*(const Vector& v) const {
    return mtl::vec::mat_cvec_multiplier
      <IdentityMatrix, Vector>(*this, v);
  }
};

/** The number of elements in the matrix. */
inline std::size_t size(const IdentityMatrix& A) {
  return A.n * A.n;
}

/** The number of rows in the matrix. */
inline std::size_t num_rows(const IdentityMatrix& A) {
  return A.n;
}

/** The number of columns in the matrix. */
inline std::size_t num_cols(const IdentityMatrix& A) {
  return A.n;
}

/** Traits that MTL uses to determine properties of our IdentityMatrix. */
namespace mtl {
  namespace ashape {

    /** Define IdentityMatrix to be a non-scalar type. */
    template<>
    struct ashape_aux<IdentityMatrix> {
      typedef nonscal type;
    };
  }

  /** IdentityMatrix implements the Collection concept
   *  with value_type and size_type. */
  template<>
  struct Collection<IdentityMatrix> {
    typedef double value_type;
    typedef unsigned size_type;
  };
}


int main()
{  
  // Construct an IdentityMatrix and "solve" Ix = b
  // using MTL's conjugate gradient solver

  const int size = 10, N = size * size;
  // Set up the identity matrix
  typedef IdentityMatrix im;
  im I(N);

  // Create an identity preconditioner
  itl::pc::identity<im> P(I);

  // Set b such that x == 1 is solution; start with x == 0
  mtl::dense_vector<double> x(N, 1.0), b(N);
  b = I * x;
  x = 0;

  // Termination criterion: r < 1.e-6
  itl::cyclic_iteration<double> iter(b, 100, 1.e-6);

  // Solve Ix == b with preconditioner P
  cg(I, x, b, P, iter);
  std::cout << x << std::endl;

  return 0;
}
