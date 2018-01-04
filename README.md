# CarND-Controls-MPC Overview

This repository contains the code needed to complete the MPC project for Udacity's Self-Driving Car Nanodegree.  When run, the executable created in this project feeds the steering angle to a simulated vehicle driving in the [Term 2 Simulator](https://github.com/udacity/self-driving-car-sim/releases).

---

## Dependencies

* cmake >= 3.5
* make >= 4.1(mac, linux), 3.81(Windows)
* [uWebSockets](https://github.com/uWebSockets/uWebSockets)
* [Udacity Term 2 Simulator](https://github.com/udacity/self-driving-car-sim/releases) in the classroom.


## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./pid`.

---

## MPC Information

### The Model

The kinematic model of the Model Predictive Control used in this project incorporates six states:
`[x, y, psi, velocity, CTE, psi error]`
The first four states and their updates equations are the same global kinematic model first presented in **Lesson 18: Vehicle Models**.  CTE is the predicted distance of the vehicle from the trajectory.  Psi error is the predicted difference between vehicle orientation and trajectory orientation.  Both are updated by adding in the delta error due to movement to the current value.
The actuators in this model are steering angle and throttle.  The model limits the steering angle (-25 to 25) and throttle (-1 to 1).
The solver of this project is setup to constrain each of the N calculated values of the model states to zero.

### Cost Function

The cost function used for the optimization of this control problem places the following importance, in order, on the states and actuators:
1. Minimize the cross-track error
2. Minimize the psi error
3. Lessen the use of the steering angle
4. Lessen the use of the throttle
5. Minimize the change between sequential steering and throttle actuations
6. Avoid stopping or driving too fast

The aggregate cost should keep the car on the road, with a smooth ride.

### Latency

Because the car is in motion, the actuation adjustments based on the most recent readings will be dated when they take effect within the system.  This latency can be accounted for by predicting the position of the car at the time when the actuation commands finish propagating through the system.  The vehicle model, described above, can be used to make the predictions that account for movement:

```
const double LATENCY = 0.1;
x += v*cos(psi)*LATENCY;
y += v*sin(psi)*LATENCY;
psi -= v*delta/LF*LATENCY;
v += acceleration*LATENCY;
```
*Note: thanks to the orientation of the simulator, a negative psi is a left turn.*

### Tuning

For the MPC hyperparameters, I started with the guidelines set in the lecture:
1. To maintain logical predictions, the prediction horizon (T) needs to be less than 2 seconds.
2. To minimize the computational cost, keep the number of timesteps (N) small.
3. To avoid discretization error, try smaller timesteps (dt).

An interesting impact of latency is that it provides a lower bound for dt.  Therefore to maintain the goal of T, the following values were used:
`N = 10; dt = 0.11;`

[iCTE]: ./data/cte.jpg "CTE Plot"
[iepsi]: ./data/epsi.jpg "epsi Plot"

![iCTE]

![iepsi]

Looking at the errors recorded in the plots above, the values are correctly hovering around zero.