
Linking and Clustering Dikes
===========================================================

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

[![DOI](https://zenodo.org/badge/272334230.svg)](https://zenodo.org/badge/latestdoi/272334230)


About
-----

Linking and clustering Dike segments uses the Hough Transform and Agglomerative clustering to link together line segments spread over Cartesian space in slope and intercept space. This software was designed to link dissected segments of dikes in Large Igneous Provinces mainly the Columbia River Basalt Group and the Deccan Traps. However, any type of line segment can be linked and clustered and this software has applications in geosciences and planetary sciences. 


Installation instructions
-------------------------

This software uses the conda system to manage python packages. You must have conda, mamba or another env manager before you can run these files. We prefer mamba, see documentation [here](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html).
Download the github repository using the web interface or commandline. Then make the python environment in the repository directory using the command line command: 

```
mamba create dikes_linking --file dikes_linking.yaml
```

Running and extending Linking and Clustering Dikes
---------------------------------------

You can start to use this code using the demostration file "DemoFile_DikeLinking.py". 

New projects should begin by placing a copy of this repository into your project repository and modifying the contents of each file as appropriate.


Citing Linking and Clustering Dikes
------------------------
This code was used to create the results of [Kubo Hutchison et al., 2023]() published in G Cubed. Please refer to the paper for the interpretation of dike linking and application to Large Igneous Provinces. 

Linking and Clustering Dikes is free to use and does not require citation or acknowledgement.

Contributing to our software
----------------------

Linking and Clustering Dikes is a community project that lives by the participation of its
members — i.e., including you!

Feedback and support
----------------------

For support, please make a comment on the repository or email me at akubo@uoregon.edu

To contribute to this project, see [CONTRIBUTING.md](CONTRIBUTING.md).

License
-------
This work is licensed under the MIT License, see LICENSE for details.

This SOFTWARE_TEMPLATE is published under the [MIT license](LICENSE).
