# modelling_tank_wall_reflections

This is a collection of functions which can be used to model the expected acoustic signal at a receiver placed in a rectangular tank given a signal produced at a source within the tank. The method of images is used to calculate the various possible reflection paths. The calculation is based on the theoretical forumlation laid out in Allen, J. B., & Berkley, D. A. (1979). Image method for efficiently simulating small‐room acoustics. The Journal of the Acoustical Society of America, 65(4), 943-950.

The file `tank_propagation_example.m` provides a simple example use case, in which the expected received signals are computed for several different source waveforms with a fixed receiver and source geometry.
