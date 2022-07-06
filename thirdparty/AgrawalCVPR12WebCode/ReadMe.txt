

CVPR 2012 Paper Matlab Code and Supplementary Materials

Title: A Theory of Multi-Layer Flat Refractive Geometry

Authors: Amit Agrawal, Srikumar Ramalingam, Yuichi Taguchi, Visesh Chari
 

1. Supplementary.pdf

Contains detailed derivation for

a) Estimation of unknown refractive index, layer thicknesses and translation along axis for 
	Case 1, Case 2,  and Case 3. 
	Case 1 and 2 results in a 6th degree equation. Case 3 was too difficult to solve as shown and we could not obtain an analytical equation.

b) Detailed derivation of axis estimation for calibartion using planar grid. We show how the third column of E matrix can be estimated from the first two columns using Demazure constraints.


c) Detailed derivation of forward projection equations for 
	Case 1, Case 2,  and Case 3. 
	Case 1 and 2 results in a 4th degree equation. Case 3 results in a 12th degree equation


2. Matlab Codes for forward projection equations
	Requires Matlab symbolic toolbox. Tested on Windows XP with Matlab 7.9.0 (R2009b)
	If Matlab give _mlans as output, run the file again

a. At matlab command prompt

>ForwardProjectionCase1
This code outputs the forward projection equation for Case 1

>ForwardProjectionCase2
This code outputs the forward projection equation for Case 2

>ForwardProjectionCase3
This code outputs the forward projection equation for Case 3


b. We also provide Matlab code that computes the 6th degree equation for Case 1 to find unknown refractive index, layer thickness and translation along axis

At matlab command prompt

>UnknownMuCase1

gives the 6th degree equation for Case 1


3. Test codes that demonstrate how to use the forward projection equations to compute the 2D projection of a 3D point

At matlab command prompt

>TestForwardProjectionCase1  test projection of 1000 3D points for case1


>TestForwardProjectionCase2  test projection of 1000 3D points for case2

>TestForwardProjectionCase3  test projection of 1000 3D points for case3



4. Complete Code for Real dataset for calibration













 

