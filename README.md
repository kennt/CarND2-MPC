# Model Predictive Control (MPC) Project Writeup
Self-Driving Car Engineer Nanodegree Program

---

This is the implementation of an MPC controller for the Udacity Self-Driving Car Engineer Nanodegree Program.

## Model Description

### State variables
The state vector contains six variables:

* __x__ - the x position of the car, relative to the current position of the car.
* __y__ - the y position of the car, relative to the current position of the car
* __psi__ - this is a measurement of the yaw, relative to the current position of the car.
* __v__ - this is the velocity of the car
* __cte__ - this is the Cross Track Error. For this model, this is the difference in the y coordinate between the current position and the desired position.
* __epsi__ - this is the error in psi (the yaw) between the current yaw and the desired yaw.


### Actuator variables
We also have two actuators:

* __delta__ - this is the change applied to the steering (__psi__)
* __a__ - this is the change applied to the throttle (__v__)


### Model parameters
Some other model parameters (these are hyperparameters that are defined outside of the model).

* __dt__ - the elapsed time for each model step
* __N__ - the total number of steps
* __Lf__ - distance from the center of mass of the car to the front axle (experimentally derived for this project and was provided to us).

### Model equations
The equations used to update the model
(Note that values such as x0 indicate the state at time t and x1 indicate the state at time t+1).

```
// Some helper values
// f0 is the value of the desired curve, evaluated at x0
f0 = f(x0)

// deriv0 is the value of the first derivative of the desired
// curve, evaluated at x0
deriv0 = f'(x0)
psides0 = atan(deriv0)

x1 = x0 + v0 * cos(psi0) * dt
y1 = y0 + v0 * sin(psi0) * dt
psi1 = psi0 - v0 * delta0 / Lf * dt
v1 = v0 + v0 * a0 * dt
cte1 = (y0 - f0) + v0 * sin(epsi0) * dt
epsi1 = (psi0 - psides0) + v0 * delta0 / Lf * dt;
```

This are the model equations.  The solver will attempt to minimize a
cost function by manipulating the actuators (__delta__ and __a__). The cost function used is discussed later.

## Cost Function
The basic cost function is provided here.  Note that the individual components of the cost function are weighted.  The initial values for the weights were taken from the Q&A session for this project.

Initial values:

* __CTE__ : 2000
* __EPSI__ : 2000
* __v__ : 1
* __delta__ : 5
* __a__ : 5
* __delta change__ : 200
* __a change__ : 10

These actually worked fine (at a target speed of 60 for my model). But as I tried higher speeds the car would start to miss some turns, so I increased the CTE and EPSI weights.  So I increased those two values to 3000. 

```
    // Initialize the cost
    fg[0] = 0.0;

    for (int t=0; t < N; t++)
    {
      // Penalize failure to meet reference state
      // Add a penalty for any CTE
      fg[0] += 3000 * CppAD::pow(vars[cte_start+t], 2);

      // Add a penalty for deviation in the angle
      fg[0] += 3000 * CppAD::pow(vars[epsi_start+t], 2);

      // Add a penalty for failing to meet the target speed
      fg[0] += 1 * CppAD::pow(vars[v_start+t] - target_v, 2);
    }

    // Actuator cost penalty
    for (int t=0; t < N-1; t++)
    {
      fg[0] += 5 * CppAD::pow(vars[delta_start + t], 2);
      fg[0] += 5 * CppAD::pow(vars[a_start + t], 2);
    }

    // Try to minimize the actuator jumps
    for (int t=0; t < N-2; t++)
    {
      fg[0] += 200 * CppAD::pow(vars[delta_start+t+1] - vars[delta_start+t], 2);
      fg[0] += 10 * CppAD::pow(vars[a_start+t+1] - vars[a_start+t], 2);
    }

```

## Polynomial Fitting

Given the waypoints, a polynomial was fitted to the waypoints using the provide `polyfit()` function.  Initially I used a cubic polynomial.  This would fit pretty well most of the time, but on some of the curves the polynomial would be very unstable (This usually happens on the 2nd or 3rd lap of the track).

[Unstable turn with cubics](./videos/cubic.mp4)


To reduce this problem, I instead used a quadratic polynomial for fitting.  This does result in a more inexact fit but is more stable and reduced the instability compared to the cubic.

Also, the impact of the other hyperparameters (__N__ and __dt__) also helped to smooth out the curves as well as making the car take sharper turns.

## Timestep Length (N) and Elapsed Duration (dt)

My initial values were

* __N__ = 10 and __dt__ = 0.1

I didn't have that much problem at this point. But I tried other settings to test the impact:

* __N__ = 20 and __dt__ = 0.1

With a longer timestep length, this helped to smooth the curve that the car took, but also increased the amount of processing needed.

* __N__ = 10 and __dt__ = 0.05

This __dt__ is less than the latency and the car was very unstable, causing the car to go off the track pretty early in the run.

* __N__ = 20 and __dt__ = 0.2

Given the duration length (__dt__), this caused the car to take some of the turns sharper (the car turned into the turn earlier).  On some of the sharper turns this would cause the car to run over the edge.

After some trial and error, I found that __N__ = 12 and __dt__ = 0.11 works well.  I increased __N__ to smooth out the curve.  I also increased __dt__ slightly so that the car would be a little more aggressive when taking the turns (0.15 was still too aggressive, but 0.11 seems to work well).

## MPC Preprocessing

The only preprocessing we do is to convert the waypoints from global (map) coordinates into local (car) coordinates.

It's easier to do the calculations in car coordinates.  Doing this makes a lot of the starting values 0, because the value are relative to the car's position (for example, x0 = y0 = psi0 = 0).

We do a translation by (-px, -py), px and py are the global coordinates of the car at the current time.  We then do a rotation by -psi to orient the car along the x axis.

```
          double cospsi = cos(-psi);
          double sinpsi = sin(-psi);
          for (size_t i=0; i<ptsx.size(); i++)
          {
            double x = ptsx[i] - px;
            double y = ptsy[i] - py;
            car_ptsx[i] = x * cospsi - y * sinpsi;
            car_ptsy[i] = x * sinpsi + y * cospsi;
          }
```

## Latency

We also incorporate the actuator latency into the system.

We do this by pushing the state of the system 100ms (the latency) into the future.  From our starting point (x=0, y=0, from the car's point of view), we assume a constant velocity, recalculate our initial state, and then solve the model from that point on.

Basically, we're taking our model equations given above, and running them for a 100ms, and using that as our initial starting point.

```
          double x_latency = v * cos(delta) * latency_secs;
          double y_latency = v * sin(delta) * latency_secs;
          double psi_latency = -v * delta * latency_secs / Lf;
          
          // evaluate the polynomial at x_latency
          double cte_latency = -polyeval(coeffs, x_latency);
          
          // evaluate the 1st derivative of the polynomial at x_latency
          double epsi_latency = -atan(deriveval(coeffs, x_latency));
```

## Results

I've managed to reach a top speed of 108mph. This is achieved on the second lap of a run (right at the end of the first bridge). This occurs  after the first run, since the car has to start from a high enough speed to be able to accelerate to 100mph before slowing down for the multiple turns.

Video of car at a target speed of 120 (top speed reached was ~109).  This is a video of a little more than 2 full laps.  The top speed is reached about 46 secs into the video.

[Target speed = 120](./videos/target120.mp4)
