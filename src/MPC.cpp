#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// The timestep length and duration
size_t N = 20;
double dt = 0.05;

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
const double LF = 2.67;

// NOTE: feel free to play around with this
// or do something completely different
double REF_F = 40;

// Constants to help locate variables in a singular vector.
const size_t X_START = 0;
const size_t Y_START = X_START + N;
const size_t PSI_START = Y_START + N;
const size_t V_START = PSI_START + N;
const size_t CTE_START = V_START + N;
const size_t EPSI_START = CTE_START + N;
const size_t DELTA_START = EPSI_START + N;
const size_t A_START = DELTA_START + N - 1;  //TODO why -1?

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs)
  {
    this->coeffs = coeffs;
  }

  typedef CPPAD_TESTVECTOR(AD<double>)ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // First in the fg vector is the Cost
    // The part of the cost based on the reference state.
    for (size_t t = 0; t < N; t++) {
      fg[0] += CppAD::pow(vars[CTE_START + t], 2);
      fg[0] += CppAD::pow(vars[EPSI_START + t], 2);
      fg[0] += CppAD::pow(vars[V_START + t] - REF_F, 2);
    }

    // Minimize the use of actuators.
    for (size_t t = 0; t < N - 1; t++) {
      fg[0] += CppAD::pow(vars[DELTA_START + t], 2);
      fg[0] += CppAD::pow(vars[A_START + t], 2);
    }

    // Minimize the value gap between sequential actuations.
    for (size_t t = 0; t < N - 2; t++) {
      fg[0] += 200 * CppAD::pow(vars[DELTA_START + t + 1] - vars[DELTA_START + t], 2);
      fg[0] += CppAD::pow(vars[A_START + t + 1] - vars[A_START + t], 2);
    }

    // Initial constraints
    fg[1 + X_START] = vars[X_START];
    fg[1 + Y_START] = vars[Y_START];
    fg[1 + PSI_START] = vars[PSI_START];
    fg[1 + V_START] = vars[V_START];
    fg[1 + CTE_START] = vars[CTE_START];
    fg[1 + EPSI_START] = vars[EPSI_START];

    // The rest of the constraints
    for (size_t t = 1; t < N; t++) {
      AD<double> x1 = vars[X_START + t];
      AD<double> y1 = vars[Y_START + t];
      AD<double> psi1 = vars[PSI_START + t];
      AD<double> v1 = vars[V_START + t];
      AD<double> cte1 = vars[CTE_START + t];
      AD<double> epsi1 = vars[EPSI_START + t];

      AD<double> x0 = vars[X_START + t - 1];
      AD<double> y0 = vars[Y_START + t - 1];
      AD<double> psi0 = vars[PSI_START + t - 1];
      AD<double> v0 = vars[V_START + t - 1];
      AD<double> cte0 = vars[CTE_START + t - 1];
      AD<double> epsi0 = vars[EPSI_START + t - 1];

      AD<double> delta0 = vars[DELTA_START + t - 1];
      AD<double> a0 = vars[A_START + t - 1];
      AD<double> f0 = coeffs[0] + coeffs[1] * x0;
      AD<double> psi_desired0 = CppAD::atan(coeffs[1]);

      fg[1 + X_START + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + Y_START + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + PSI_START + t] = psi1 - (psi0 + v0/LF * delta0 * dt);
      fg[1 + V_START + t] = v1 - (v0 + a0 * dt);
      fg[1 + CTE_START + t] = cte1 - ((f0 - y0) + v0 * CppAD::sin(epsi0) * dt);
      fg[1 + EPSI_START + t] = epsi1 - ((psi_desired0 - psi0) - v0/LF * delta0 * dt);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC()
{
}
MPC::~MPC()
{
}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs)
{
  bool ok = true;
  typedef CPPAD_TESTVECTOR(double)Dvector;

  // Set the number of model variables (includes both states and inputs).
  size_t n_vars = N * 6 + (N - 1) * 2;
  // Set the number of constraints
  size_t n_constraints = N * 6;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (size_t i = 0; i < n_vars; i++) {
    vars[i] = 0.0;
  }

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

  // Set all non-actuators upper and lower limits
  // to the max negative and positive values.
  for (size_t i = 0; i < DELTA_START; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  for (size_t i = DELTA_START; i < A_START; i++) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }

  // Acceleration/decceleration upper and lower limits.
  for (size_t i = A_START; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (size_t i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

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
  CppAD::ipopt::solve<Dvector, FG_eval>(options, vars, vars_lowerbound,
                                        vars_upperbound, constraints_lowerbound,
                                        constraints_upperbound, fg_eval,
                                        solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // Return the first actuator values
  return {solution.x[X_START + 1], solution.x[Y_START + 1],
    solution.x[PSI_START + 1], solution.x[V_START + 1],
    solution.x[CTE_START + 1], solution.x[EPSI_START + 1],
    solution.x[DELTA_START], solution.x[A_START]};
}
