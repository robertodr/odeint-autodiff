#include <iostream>
#include <cmath>

#include <taylor.hpp>

#include <boost/numeric/odeint.hpp>

template <typename T> class pendulum {

private:
  /// Length of the pendulum
  double length_;
  /// Gravitational constant (little g)
  double g_;

public:
  pendulum(double l, double g = 9.81) : length_(l), g_(g) {}

  void operator()(const T &theta, T &dtheta_dt, const double /* t */) {
    dtheta_dt[0] = theta[1];
    dtheta_dt[1] = -(g_ / length_) * theta[0];
  }
};

class pendulum_exact {
private:
  /// Length of the pendulum
  double length_;
  /// Initial position (angle)
  double theta0_;
  /// Gravitational constant (little g)
  double g_;

public:
  pendulum_exact(double l, double theta0, double g = 9.81) : length_(l), theta0_(theta0), g_(g) {}
  double position(double t) {
    return theta0_ * std::cos(std::sqrt(g_/length_) * t);
  }
  double speed(double t) {
    return (-theta0_ * std::sqrt(g_/length_) * std::sin(std::sqrt(g_/length_) * t));
  }
};

int main() {
  using namespace std;
  using namespace boost::numeric::odeint;

  typedef std::vector< taylor<double, 1, 1> > state_type;

  //[ state_initialization
  state_type x(2);
  x[0] = 1.0; // start at x=1.0, p=0.0
  x[1] = 0.0;
  //]

  double length = 9.81;
  pendulum<state_type> approx(9.81);
  pendulum_exact exact(9.81, 1.0);

  runge_kutta4< state_type> stepper;
  integrate_const( stepper, approx, x , 0.0 , 10.0 , 0.01 );

  std::cout << "x.size() " << x.size() << std::endl;

  std::cout << exact.position(10.0) << "   " << exact.speed(10.0) << std::endl;
  std::cout << x[0] << std::endl;
  std::cout << x[1] << std::endl;

  return 0;
}
