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
    dtheta_dt[0][0] = theta[0][1];
    dtheta_dt[0][1] = -(g_ / length_) * theta[0][0];
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
  x[0][0] = 1.0; // start at x=1.0, p=0.0
  x[0][1] = 0.0; // start at x=1.0, p=0.0
  x[1][0] = 0.0;
  x[1][1] = 0.0;
  //]

  double length = 9.81;
  pendulum<state_type> approx(9.81);
  pendulum_exact exact(9.81, 1.0);

  runge_kutta4< state_type> stepper;

  const double dt = 0.01;
  const double tstart = 0.0;
  const double tend = 1.0;
  for( double t = tstart ; t < tend ; t += dt )
  {
    std::cout << t << "   " << std::setprecision(16) << exact.position(t) << "   " << std::setprecision(16) << exact.speed(t) << std::endl;
    std::cout << t << "   " << std::setprecision(16) << x[0]  << "   " << std::setprecision(16) << x[1] << std::endl;
    stepper.do_step( approx , x , t , dt );
  }

  return 0;
}
