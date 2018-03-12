#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include <algorithm>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial (with coeffs) at a point x
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Evaluate the first deriviative of a polynomial (with coeffs) at a point x
double deriveval(Eigen::VectorXd coeffs, double x){
  double result = 0.0;
  for (int i = 1; i < coeffs.size(); i++) {
    result += i * coeffs[i] * pow(x, i-1);
  }
  return result;
}
// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

int main() {
  uWS::Hub h;

  // MPC is initialized here!
  MPC mpc(target_speed);

  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v = j[1]["speed"];
          double delta = j[1]["steering_angle"];
          double a = j[1]["throttle"];

          /*
          * TODO: Calculate steering angle and throttle using MPC.
          *
          * Both are in between [-1, 1].
          *
          */

          /*
          * Data Preprocessing
          * Convert the waypoints from global coordinates into
          * car coordinates.
          */
          // Note that we are also converting the containers from
          // STL std::vector to Eigen::VectorXd
          Eigen::VectorXd car_ptsx(ptsx.size());
          Eigen::VectorXd car_ptsy(ptsy.size());
          double cospsi = cos(-psi);
          double sinpsi = sin(-psi);
          for (size_t i=0; i<ptsx.size(); i++)
          {
            double x = ptsx[i] - px;
            double y = ptsy[i] - py;
            car_ptsx[i] = x * cospsi - y * sinpsi;
            car_ptsy[i] = x * sinpsi + y * cospsi;
          }

          /*
          * Polynomial fitting
          * Fit a polynomial to this vector to get the initial coefficients
          */
          Eigen::VectorXd coeffs =  polyfit(car_ptsx, car_ptsy, 2);
          Eigen::VectorXd state(6);

          /*
          * Incorporate latency into the model
          * Build the initial state vector
          * The latency effects are incorporated into the initial state vector.
          */
          double x_latency = v * cos(delta) * latency_secs;
          double y_latency = v * sin(delta) * latency_secs;
          double psi_latency = -v * delta * latency_secs / Lf;
          double cte_latency = -polyeval(coeffs, x_latency);
          double epsi_latency = -atan(deriveval(coeffs, x_latency));

          state << x_latency, y_latency, psi_latency, v, cte_latency, epsi_latency;

          /*
          * Run the solver.
          */
          std::vector<double> params = mpc.Solve(state, coeffs);

          json msgJson;
          double steer_value = params[0];

          // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
          // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
          msgJson["steering_angle"] = steer_value/deg2rad(25);

          //Extract the x,y values and display the MPC predicted trajectory
          vector<double> mpc_x_vals(&params[2], &params[2+N]);
          vector<double> mpc_y_vals(&params[2+N], &params[2+2*N]);
          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;

          //Display the waypoints/reference line
          vector<double> next_x_vals(&car_ptsx[0], car_ptsx.data() + car_ptsx.cols() * car_ptsx.rows());
          vector<double> next_y_vals(&car_ptsy[0], car_ptsy.data() + car_ptsy.cols() * car_ptsy.rows());

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Yellow line
          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

          double throttle_value = params[1];
#if 0
          // Use this to turn down the throttle in a tight turn
          if (v > 20) {
            // Depending on the steer value, adjust the throttle
            // (the tighter the turn the less throttle)
            // throttle-steer curve
            //                            0     1     2     3     4     5     6     7      8       9       10
            double throttle_adjust[] = { 1.00, 1.00, 0.70, 0.50, 0.30, 0.10, 0.03, -0.001, -0.002, -0.003, -0.007 };
            double constrained_steer_value = steer_value/deg2rad(25);
            throttle_value *= throttle_adjust[static_cast<int>(fabs(constrained_steer_value*10))];
          }
#endif
          msgJson["throttle"] = throttle_value;

          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          std::cout << msg << std::endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          this_thread::sleep_for(chrono::milliseconds(100));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
