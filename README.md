First repository, so please don't get mad at me. I'm really sensitive.

"Squigglies" is a short, friendly, simple name for a somewhat more complex and specific discipline referred to as instability-assisted fused deposition modeling. A single piece of code is now left, which elaborately solves a series of step operations in one go. This code is part of my bachelor thesis: *Numerical evaluation of periodic forcing influence in the Fluid Mechanical Sewing Machine problem*. For further information on the topic, check the PDF file 'TFG'.

- ~~"Sinusoidal.m" contains the MATLAB code used to numerically solve the Geometrical Model of the movement produces by a viscous thread falling onto a moving belt, while being transversally perturbated in a sinusoidal-like motion. The movement is dependent on three dimensionless parameters: quotients of speed, time and length. The code includes the declaration of the ODE function, as well as the plotting of the results.~~
17/04/2020
"Sinusoidal.m" does not exist anymore. It was out of date. The expected approach -and a greater insight- is performed in the "Complete_FFT.m" file.
29/05/2020

- "Complete_FFT.m" contains the MATLAB code used to numerically solve the standard Geometrical Model; that is, the one established in "Sinusoidal.m" when forcing is null (LR = 0). It also provides its FFT scheme and determines the SR-pattern frequency (defined as the first primary Y-axis frequency). Lastly, while using this said frequncy to compare with the perturbation frequency, it solves the forced ODE -as done in "Sinusoidal.m". Again, a FFT graph is shown.
25/05/2020
