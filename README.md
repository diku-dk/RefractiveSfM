# RefractiveSfM

This repository contains code for refractive structure-from-motion.
## Introduction

Example code for verify  the algorithms described in:

> Absolute and Relative Pose Estimation in Refractive Multi View
> In Submission to ICCV 2021

## Side contribution
1. We extend Agrawal's open-source implementation to handle single refraction and double refraction with three different media. 
2. We implement the algorithms described in the following articles. Furthermore, we accelerate them using the virtual camera transformation.
> Chadebecq, François, et al. "Refractive two-view reconstruction for underwater 3d vision." International Journal of Computer Vision 128.5 (2020): 1101-1117.
> Chadebecq, François, et al. "Refractive structure-from-motion through a flat refractive interface." Proceedings of the ieee international conference on computer vision. 2017.

## How-to-use

`main_absolute_pose_toy_example.m` : generate synthetic data for evaluating absolution pose estimation.

`main_relative_pose_toy_example.m` : generate synthetic data for evaluating relative pose estimation.

### Dependencies:

Some certain functions have dependencies on [mexopencv](https://github.com/kyamagu/mexopencv) library.

Please following the instructions mentioned in the  [mexopencv](https://github.com/kyamagu/mexopencv) repository for compilation and installation.

UPnP and GPnP rely on installation of the OpenGV library.

## Acknowledgements

In folder `thirdparty`, we include some third-party software made by:

> 1. Amit Agrawal, Srikumar Ramalingam, Yuichi Taguchi, and Visesh Chari. A theory of multi-layer flat refractive geometry. In Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition, pages 3346–3353. IEEE, 2012.
> 2. Sebastian Haner and Kalle Astrom. Absolute pose for cameras under flat refractive interfaces. In Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition, pages 1428–1436, 2015.
>3. J. Ventura, C. Arth, G. Reitmayr, and D. Schmalstieg, A Minimal Solution to the Generalized Pose-and-Scale Problem, in CVPR '14: Proceedings of the 2014 IEEE Conference on Computer Vision and Pattern Recognition, 2014.
> 4. P. Miraldo, T. Dias, S. Ramalingam (2018), A Minimal Closed-Form Solution for Multi-Perspective Pose Estimation using Points and Lines, European Conf. Computer Vision (ECCV) [arXiv];


## Copyright

> This program is free software: you can redistribute it and/or modify
>  it under the terms of the version 3 of the GNU General Public License
>  as published by the Free Software Foundation.
>
>  This program is distributed in the hope that it will be useful, but
>  WITHOUT ANY WARRANTY; without even the implied warranty of
>  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
>  General Public License for more details.       
>  You should have received a copy of the GNU General Public License
>  along with this program. If not, see <http://www.gnu.org/licenses/>.





