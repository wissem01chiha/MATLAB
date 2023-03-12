Least Recursive Square Method and Kalman Filter for First-Order Systems

This MATLAB code implements two algorithms for estimating the parameters of a first-order system: the least recursive square (LRS) method and the Kalman filter. Both algorithms use a set of input and output data to estimate the system's time constant and gain.

Usage
To use this code, simply run the main.m script in MATLAB. This script reads in the input and output data from the file data.csv, and then applies both the LRS method and Kalman filter to estimate the system parameters. The estimated parameters are then plotted on a graph along with the original data.

Requirements
This code requires MATLAB and the Control System Toolbox.

LRS Method
The LRS method estimates the system parameters recursively using a set of input and output data. The estimated parameters at each time step are used to update the estimate for the next time step. This process continues until all data has been processed.

Kalman Filter
The Kalman filter is an algorithm that estimates the state of a system based on a set of noisy input and output data. In this implementation, the Kalman filter is used to estimate the time constant and gain of the first-order system.

Results
The main.m script generates a plot that shows the original data along with the estimated time constant and gain using both the LRS method and Kalman filter. The plot also shows the error between the estimated parameters and the true parameters of the system.

License
This code is released under the MIT license. See the LICENSE file for more information.
