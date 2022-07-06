Generalized Pose-and-Scale Solver
=================================

### Authors
Jonathan Ventura, Clemens Arth, Gerhard Reitmayr and Dieter Schmalstieg

### Introduction
This repository contains MATLAB code for solving the generalized pose-and-scale problem.
The problem and our solution are described in our CVPR 2014 paper.

### Usage
	GENPOSEANDSCALE Solves for the pose and scale of a generalized camera.
	  [Psolns,ssolns] = genposeandscale(p,d,q)

	  p is a 3xN matrix of ray origins
	  d is a 3xN matrix of ray directions
	  q is a 3xN matrix of anchor points

### Acknowledgment
If you use this algorithm in your work, please cite the following paper:

J. Ventura, C. Arth, G. Reitmayr, and D. Schmalstieg,
A Minimal Solution to the Generalized Pose-and-Scale Problem, in CVPR '14: Proceedings of the 2014 IEEE Conference on Computer Vision and Pattern Recognition, 2014.

### License
Copyright (c) 2014 Graz University of Technology.

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted (subject to the limitations in the disclaimer
below) provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

 * Neither the name of Graz University of Technology nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY THIS
LICENSE.  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
