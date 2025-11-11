# Numerical Methods for Differential Equations (ODEs)

This repository contains the MATLAB code and final report for a laboratory assignment on numerical methods for solving Ordinary Differential Equations (ODEs).

## 📝 Overview of Exercises

The project is divided into six main exercises, each exploring a different aspect of numerical methods.

### Exercise 1: 2-Step Simpson Method
* [cite_start]**Objective:** To solve the test equation $y' = -5y$ [cite: 7] [cite_start]using the implicit 2-step Simpson method[cite: 5].
* [cite_start]**Analysis:** This multi-step method requires a second initial value ($y(1)$) to start[cite: 28]. This value is computed in three different ways for comparison:
    1.  [cite_start]The exact solution ($e^{-5h}$) [cite: 29]
    2.  [cite_start]A 4th-order Runge-Kutta (RK4) step [cite: 30]
    3.  [cite_start]A Forward Euler step [cite: 31]
* [cite_start]**Result:** The error and stability of the solution resulting from each of the three starting values are analyzed and plotted[cite: 43, 85].

### Exercise 2: 4th-Order Runge-Kutta (RK4)
* [cite_start]**Objective:** To solve the equation $y' = -10y^2$ [cite: 90] [cite_start]using the 4-stage RK4 method[cite: 89].
* [cite_start]**Analysis:** The convergence of the method is studied[cite: 93]. [cite_start]The error at the final time $T=2$ is plotted on a log-log scale as a function of the step size $h$ (for $h=2^{-k}$, $k=5...10$)[cite: 93, 94].
* [cite_start]**Result:** The slope of the resulting line in the log-log plot confirms the expected theoretical convergence order of $O(h^4)$ for the RK4 method[cite: 112, 113].

### Exercise 3: Richardson Extrapolation
* [cite_start]**Objective:** To derive a formula for Richardson extrapolation [cite: 117] [cite_start]as applied to the 4-stage Runge-Kutta method[cite: 116].
* [cite_start]**Analysis:** Starting from the local truncation error formula for RK4, which is $O(h^5)$[cite: 122, 125], a system of equations is solved to cancel this leading error term.
* [cite_start]**Result:** The derived formula for the more accurate, extrapolated solution is $y^* = \frac{-y_{n+1}^h + 16y_{n+2}^{h/2}}{15}$[cite: 137].

### Exercise 4: BDF2 Method
* [cite_start]**Objective:** To derive the BDF2 (2-step Backward Differentiation Formula)[cite: 145].
* [cite_start]**Analysis:** The formula is derived by finding the Lagrange interpolating polynomial $P(t)$ that passes through the points $(t_n, y_n)$, $(t_{n+1}, y_{n+1})$, and $(t_{n+2}, y_{n+2})$[cite: 152, 153, 156]. [cite_start]The derivative $y'(t_{n+2})$ is then approximated by $P'(t_{n+2})$[cite: 162, 163].
* [cite_start]**Result:** The BDF2 formula is derived, and its local truncation error [cite: 171] [cite_start]and region of absolute stability [cite: 181] are analyzed.

### Exercise 5: Linear ODE Systems
* [cite_start]**Objective:** To solve the stiff linear system of ODEs $y'(t) = -Ay(t)$[cite: 198].
* **Analysis:** Three different methods are implemented and compared:
    1.  [cite_start]`ode45` (MATLAB's built-in adaptive solver) [cite: 225]
    2.  [cite_start]Crank-Nicolson (CN) method [cite: 227]
    3.  [cite_start]BDF3 (3-step BDF) method [cite: 236]
* [cite_start]**Result:** The implicit methods (CN and BDF3) are solved using the Preconditioned Conjugate Gradient (PCG) method[cite: 227, 233, 245]. [cite_start]The infinity-norm error and total CPU time are compiled in a table for comparison[cite: 247, 248].

### Exercise 6: Predator-Prey Model (Lotka-Volterra)
* [cite_start]**Objective:** To solve the Lotka-Volterra equations, a classic model for predator-prey dynamics[cite: 257, 259].
* [cite_start]**Parameters:** $\alpha=0.2$, $\beta=0.01$, $\gamma=0.004$, $\delta=0.07$[cite: 261].
* [cite_start]**Method:** The 4th-order Runge-Kutta (RK4) method is used to solve the system of two coupled non-linear ODEs[cite: 261, 262].
* [cite_start]**Result:** A plot is generated showing the cyclical behavior of the prey (x) and predator (y) populations over time[cite: 273, 293].

## 📂 Repository Structure

* The folders (`exercises 1...`, `ex3...`, etc.) contain the various MATLAB scripts (`.m`) used to run the simulations and generate the plots for each exercise.
