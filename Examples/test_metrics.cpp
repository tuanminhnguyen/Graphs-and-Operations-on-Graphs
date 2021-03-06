#include "Metrics.hpp"

int main()
{
  using namespace Metrics;

  // Type deduction and enforcement
  Mass m = 3;

  m *= 2.1;       // Can multiply by a scalar
  //m += 1;       // Cannot add a scalar
  m += Mass(1.2); // Can add Mass
  m += 1.2_g;     // Can add a scalar with explicit mass units

  Acceleration3 a = Point(0, 1, 2);
  Force3 f = m * a;
  std::cout << "f = m*a" << std::endl;
  std::cout << f << " = " << m << " * " << a << std::endl;
  std::cout << std::endl;
  std::cout << "a / m = " << a / m << std::endl;
  std::cout << std::endl;

  // Use of auto and user-defined literals
  Acceleration grav1 = 9.81;
  auto grav2 = 9.81_m / (1_s * 1_s);

  std::cout << "grav1 = " << grav1 << std::endl;
  std::cout << "grav2 = " << grav2 << std::endl;
  std::cout << "grav1 == grav2 : "
            << std::boolalpha << (grav1 == grav2) << std::endl;
  std::cout << std::endl;
}
