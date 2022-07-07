# RefractiveSfM

This repository contains code for refractive structure-from-motion representing implementations of the methods presented in the following scientific publications:

> Xiao Hu, Francois Lauze, Kim Steenstrup Pedersen. Refractive Pose Refinement：Generalising the geometric relation between camera and refractive interface.
> International Journal of Computer Vision (IJCV), 2022.

and

> Xiao Hu, Francois Lauze, Kim Steenstrup Pedersen, Jean Mélou. Absolute and Relative Pose Estimation in Refractive Multi View.
> 1st Workshop on Traditional Computer Vision in the Age of Deep
Learning (TradiCV), Proceedings of the 2021 IEEE/CVF International Conference on Computer Vision Workshop (ICCVW). IEEE, p. 2569-2578, 2021.

If you use this code in a scientific context, we appreciate if you cite either of the two publications above (which ever is appropriate for your purpose).

## Dataset
For the experiments in the two papers listed above we used the *DIKU Refractive Scenes Dataset 2022* which is available at:

[http://doi.org/10.17894/ucph.5d1b9bea-b105-4d43-aefb-c53df7806c2a]( http://doi.org/10.17894/ucph.5d1b9bea-b105-4d43-aefb-c53df7806c2a)


## Requirements
The code has been implemented using Matlab and should work with version R2016b and above.


### Dependencies:

Some certain functions have dependencies on [mexopencv](https://github.com/kyamagu/mexopencv) library.

Please following the instructions mentioned in the  [mexopencv](https://github.com/kyamagu/mexopencv) repository for compilation and installation.

UPnP and GPnP rely on installation of the OpenGV library.


## Side contributions
1. We extend Agrawal's open-source implementation to handle single refraction and double refraction with three different media. 
2. We implement the algorithms described in the following articles and furthermore, we accelerate them using the virtual camera transformation.
> Chadebecq, François, et al. "Refractive two-view reconstruction for underwater 3d vision." International Journal of Computer Vision 128.5 (2020): 1101-1117.

> Chadebecq, François, et al. "Refractive structure-from-motion through a flat refractive interface." Proceedings of the ieee international conference on computer vision. 2017.

## How-to-use

`main_absolute_pose_toy_example.m` : generate synthetic data for evaluating absolution pose estimation.

`main_relative_pose_toy_example.m` : generate synthetic data for evaluating relative pose estimation.


## Acknowledgements

In folder `thirdparty`, we include some third-party software made by:

> 1. Amit Agrawal, Srikumar Ramalingam, Yuichi Taguchi, and Visesh Chari. A theory of multi-layer flat refractive geometry. In Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition, pages 3346–3353. IEEE, 2012.
> 2. Sebastian Haner and Kalle Astrom. Absolute pose for cameras under flat refractive interfaces. In Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition, pages 1428–1436, 2015.
>3. J. Ventura, C. Arth, G. Reitmayr, and D. Schmalstieg, A Minimal Solution to the Generalized Pose-and-Scale Problem, in CVPR '14: Proceedings of the 2014 IEEE Conference on Computer Vision and Pattern Recognition, 2014.
> 4. P. Miraldo, T. Dias, S. Ramalingam (2018), A Minimal Closed-Form Solution for Multi-Perspective Pose Estimation using Points and Lines, European Conf. Computer Vision (ECCV) [arXiv];


## Copyright

Copyright 2022, Xiao Hu, Francois Lauze, Kim Steenstrup Pedersen

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License. 
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.


