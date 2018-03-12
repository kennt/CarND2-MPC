#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;


// TODO: Set the timestep length and duration
// Projects : N*dt seconds ahead
//  N is # iterations
//  dt is seconds/iteration
constexpr size_t N = 12;
constexpr double dt = 0.11;
constexpr double target_speed = 120;
constexpr int latency = 100;			// in milliseconds

// Some parameters
constexpr double latency_secs = latency/1000.0;	// in seconds
constexpr double latency_steps = dt / latency_secs;


// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;


// State variables
// ========
// Length of state vector
//  x, y, psi, v, cte, epsi
constexpr size_t STATE_SIZE = 6;

// Length of the actuator vector (# of actuators)
//  delta, a
constexpr size_t ACTUATOR_SIZE = 2;


class MPC {
 public:
  MPC(double target_velocity);

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);

  double m_target_v;
};

#endif /* MPC_H */
