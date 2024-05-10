Aero-Propulsive Panel Code in Two Dimensions (APPC-2D)
Authors: Himavath Jois, Phillip J. Ansell
Inquiries: hjois2@illinois.edu

/////////////////////////////////////////////////////////////////////

Published and distributed under the MIT Open License (LICENSE.txt).

Copyright © 2024 Himavath Jois, Phillip J. Ansell.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

/////////////////////////////////////////////////////////////////////

Changelog:
v1.0.1 (May 10th, 2024) - suite of changes made to boost performance and increase compatibility with GNUOctave
v1.0.0 (January 20th, 2024) - initial public release

/////////////////////////////////////////////////////////////////////

This code calculates the performance of an aero-propulsive airfoil system by modeling the propulsive actuator as a set of powered wake boundaries emanating from the trailing edge of two airfoil surfaces. The details of the implementation are found in the comments of the code as well as the following paper:

Jois, H., and Ansell, P. J., “Analytical Framework for Design of Aero-Propulsive Geometries with Powered Wakes,” AIAA SCITECH 2023 Forum, American Institute of Aeronautics and Astronautics, 2023. https://doi.org/10.2514/6.2023-1754.

Disclaimers:
The code is designed to run in MATLAB R2021a and is not guaranteed to run in older versions of MATLAB or in GNUOctave. The code is not optmized for performance but should run within 30 seconds on modern consumer-grade laptop computers.
There are currently no plans to develop this code further, though feedback would be greatly appreciated, particularly on the scientific aspects of the implementation. Feedback and/or requests for bug fixes and/or feature additions can be made to the email address at the top of this file, but the authors make no guarantee of these services. If available, future code versions will be posted to the research group website: https://ansell.aerospace.illinois.edu/. 

/////////////////////////////////////////////////////////////////////

Instructions:
A user should first download all files and begin working in designDriver.m. The following inputs are required to begin an analysis:
1. Airfoil coordinates - at least two sets of airfoil coordinates are required for the code to run properly, since a powered wake must be bounded from both top and bottom. Beyond that, you can have as many surfaces as your computer can handle. Example coordinates are given in the "Airfoils" folder.
2. System angle of attack (AoA, alpha)
3. Propulsion thrust coefficient (CT)

Once the three inputs are provided in designDriver.m, users should run the first section of code to complete the analysis. Typical configurations, such as the ones given in the examples, should converge within 10-15 iterations (iter). A user may then choose to plot streamlines by running the second section, and/or pressure distributions from the third section. Further discussion about what the code is doing can be found in the comments and the aforementioned paper. 