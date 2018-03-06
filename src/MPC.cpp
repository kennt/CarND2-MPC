#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the timestep length and duration
// Projects : N*dt seconds ahead
//  N is # iterations
//  dt is seconds/iteration
constexpr size_t N = 10;
constexpr double dt = 0.1;

// State variables
// ========
// Length of state vector
//  x, y, psi, v, cte, epsi
constexpr size_t STATE_SIZE = 6;

// Length of the actuator vector (# of actuators)
//  delta, a
constexpr size_t ACTUATOR_SIZE = 2;

//
// position variables (index into the vars array) based on N
// NOTE: This assumes a certain layout of the vars array.
//
constexpr size_t x_start = 0;
constexpr size_t y_start = x_start + N;
constexpr size_t psi_start = y_start + N;
constexpr size_t v_start = psi_start + N;
constexpr size_t cte_start = v_start + N;
constexpr size_t epsi_start = cte_start + N;
constexpr size_t delta_start = epsi_start + N;
constexpr size_t a_start = delta_start + N - 1;



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

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  double target_v;

  FG_eval(double target_v, Eigen::VectorXd coeffs)
  {
    this->coeffs = coeffs;
    this->target_v = target_v;
  }

  // Calculates the value of f(x)
  // Using the coeffs passed in durng construction
  // Functon copied from main.cpp
  AD<double> evaluate_at(AD<double>& x)
  {
    AD<double> result = 0.0;
    for (int i = 0; i < coeffs.size(); i++) {
      result += coeffs[i] * CppAD::pow(x, i);
    }
    return result;
  }

  AD<double> evaluate_derivative_at(AD<double>& x)
  {
    AD<double> result = 0.0;
    for (int i = 1; i < coeffs.size(); i++) {
      result += i * coeffs[i] * CppAD::pow(x, i-1);
    }
    return result;
  }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.

    // Calculate the cost function

    // Initialize the cost
    fg[0] = 0.0;

    // CHECK: Should the errors be normalized?  Can they be normalized?

    for (int t=0; t < N; t++)
    {
      // Penalize failure to meet reference state
      // Add a penalty for any CTE
      fg[0] += 1500 * CppAD::pow(vars[cte_start+t], 2);

      // Add a penalty for deviation in the angle
      fg[0] += 30 * CppAD::pow(vars[epsi_start+t], 2);

      // Add a penalty for failing to meet the target speed
      fg[0] += CppAD::pow(vars[v_start+t] - target_v, 2);
    }

    // Actuator cost penalty
    for (int t=0; t < N-1; t++)
    {
      fg[0] += 1 * CppAD::pow(vars[delta_start + t], 2);
      fg[0] += 1 * CppAD::pow(vars[a_start + t], 2);
    }

    // Try to minimize the actuator jumps
    for (int t=0; t < N-2; t++)
    {
      fg[0] += 20 * CppAD::pow(vars[delta_start+t+1] - vars[delta_start+t], 2);
      fg[0] += CppAD::pow(vars[a_start+t+1] - vars[a_start+t], 2);
    }


    // Initialize the rest of the fg array
    // Remember that fg[0] is the cost, so the indexes have to be
    // adjusted by 1.

    // Initialize the values for t=0
    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];

    // Now setup for t > 0
    for (int t=1; t < N; t++)
    {
      // Calculate the state at time t
      // Only need the previous state at t-1 due to Markov assumption
      //
      // NOTE: The use of `AD<double>` and use of `CppAD`!
      // This is also CppAD can compute derivatives and pass
      // these to the solver.

      // State at time t-1
      AD<double> x0 = vars[x_start + t - 1];
      AD<double> y0 = vars[y_start + t - 1];
      AD<double> psi0 = vars[psi_start + t - 1];
      AD<double> v0 = vars[v_start + t - 1];
      AD<double> cte0 = vars[cte_start + t - 1];
      AD<double> epsi0 = vars[epsi_start + t - 1];

      // State at time t
      AD<double> x1 = vars[x_start + t];
      AD<double> y1 = vars[y_start + t];
      AD<double> psi1 = vars[psi_start + t];
      AD<double> v1 = vars[v_start + t];
      AD<double> cte1 = vars[cte_start + t];
      AD<double> epsi1 = vars[epsi_start + t];

      // Actuator actions at time t-1
      AD<double> delta0 = vars[delta_start + t - 1];
      AD<double> a0 = vars[a_start + t - 1];

      // Calculate the position
      // Model equations
      // ========
      // (note xt' = xt at t+1)
      //
      // xt' = xt + vt*cos(psit)*dt
      // yt' = yt + vt*sin(psit)*dt
      // psit' = psit - (vt*deltat*dt / Lf)
      // vt' = vt + at*dt
      //
      // Change in the calculation for cte
      // See https://discussions.udacity.com/t/mistake-in-lessons-motion-equation-for-cte/458082/13
      // ctet' = yt - f(xt) + vt*sin(epsit)*dt
      //
      // epsit' = psit - psidest + (vt*deltat*dt / Lf)
      //
      // at in [-1. 1]
      // deltat in [-25deg, +25deg]
      //

      // Setup our constraints
      // Our constaints are that the difference in values from time t and t-1
      // should be equal to 0.  For example, if we have
      //    xt' = xt + vt*cos(psit)*dt
      // then we set this to 0, so
      //    0 = xt' - (xt + vt*cos(psit)*dt
      // this is our constraint
      //
      AD<double> f0 = evaluate_at(x0);
      AD<double> psides0 = CppAD::atan(evaluate_derivative_at(x0));

      fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + psi_start + t] = psi1 - (psi0 - v0 * delta0 / Lf * dt);
      fg[1 + v_start + t] = v1 - (v0 + a0 * dt);
      fg[1 + cte_start + t] = cte1 - ((y0 - f0) + v0 * CppAD::sin(epsi0) * dt);

      fg[1 + epsi_start + t] = epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC(double desired_v)
  : m_target_v(desired_v)
{}

MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd initial_state, Eigen::VectorXd coeffs) {
  bool ok = true;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  double x = initial_state[0];
  double y = initial_state[1];
  double psi = initial_state[2];
  double v = initial_state[3];
  double cte = initial_state[4];
  double epsi = initial_state[5];

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9

  // So the equation for the number of vars is:
  //   state_size*N + actuator_size*(N-1)
  size_t n_vars = STATE_SIZE * N + ACTUATOR_SIZE * (N-1);

  // TODO: Set the number of constraints
  size_t n_constraints = STATE_SIZE * N;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }
  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

  // TODO: Set lower and upper limits for variables.
  // All non-actuator limits to max negative/positive values
  // allow them full range
  for (int i=0; i<delta_start; i++)
  {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  // Limit the delta to += 25 degrees
  // Values below are in radians
  for (int i = delta_start; i < a_start; i++)
  {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }

  // Limit acceleration to += 1.0
  for (int i = a_start; i < n_vars; i++)
  {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  // Set lower and upper bound constraints for the
  // initial state.
  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;

  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;

  // object that computes objective and constraints
  FG_eval fg_eval(m_target_v, coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  return { solution.x[delta_start], solution.x[a_start] };
}
