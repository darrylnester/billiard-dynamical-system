# billiard-dynamical-system
Python code to simulate perfectly elastic collisions of 2 or more masses ("billiard balls") on a one-dimensional table, with reflection (a "bumper") at one end.

This code was written to create plots and animations for a contributed-paper presentation at the 2023 Miami University Mathematics Conference. The conference theme was "Differential Equations and Dynamical Systems and Their Applications," and my talk was entitled "A dynamical system with two (or three, or more) billiard balls." It was an extension of the system described by Dr. Gregory Galperin in his 2003 paper "Playing pool with π (the number π from a billiard point of view)"(https://www.researchgate.net/publication/266523056_Playing_pool_with_p_the_number_p_from_a_billiard_point_of_view).

I am uploading the current version to GitHub in case it might prove useful to others who want to explore this dynamical system.

billiard.py defines the function run_bounces, and several other helper functions.
* Given a system (specified as a list of [mass, position, velocity] lists), run_bounces simulates the collisions, returning a dict containing the final state of the system.

billiard_plot.py uses matplotlib, and includes functions to create several different visualizations of these systems
* two_ball_animation creates an animation of the position plots for a two-ball system (mass A = 1, mass B = m)
* the other functions assume a three-ball system (masses 1, m2 and m3, located at x1, 1, and 1.1,  with initial velocities 0, 0, and -1
    * position_plot and position_animation show the positions of the balls over time
    * show_bounces animates the positions (ending at the final collision)
    * phase_plane (and related functions) show a phase plot of the normalized velocities of the first two masses.
