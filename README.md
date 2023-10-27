
# Welcome to LinkingLines!

[![DOI](https://zenodo.org/badge/272334230.svg)](https://zenodo.org/badge/latestdoi/272334230)

[![PyPI](https://img.shields.io/pypi/v/LinkingLines.svg)](https://pypi.org/project/LinkingLines/)

[![ReadtheDocs](https://readthedocs.org/projects/LinkingLines/badge/)](https://linkinglines.readthedocs.io/)

[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

#Read the Full documentation on [ReadtheDocs!](https://linkinglines.readthedocs.io/en/latest/)#

## 1. Introduction
Welcome to the documentation for our Python module that performs the Hough
Transform on line data from a CSV, clusters it using Agglomerative Clustering,
and provides functionality to export the results into a CSV file.
This module also includes custom plotting scripts and feature extraction
methods to help you analyze and visualize your data effectively.

This code was used to create the results published in
[Kubo Hutchison et al., 2023](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2022GC010842).
Initially, it was designed to link together mapped dike segments in Cartesian space
to find their true lengths. This code can be applied to any linear features including
roads, fractures, and other types of linear data.

- **Data Clustering**: Apply Agglomerative Clustering to group similar data points, this
can be used for data reduction, analysis, and mapping .

- **Data Visualization**: Custom plotting scripts help you visualize and analyze
your data, making it easier to identify patterns and anomalies.

- **Feature Extraction**: Extract meaningful features from clustered data to
perform further analysis, such as linear or radial type features.


Full documentation can be found on [ReadTheDocs](https://linkinglines.readthedocs.io/en/latest/)

## 2. Installation
To use this module, make sure you have Python installed (preferably Python 3.x).
You can install the required packages using pip:

```bash
pip install linkinglines
```

## 3. Quick Start

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import linkinglines as ll

data=pd.read_csv('path/to/data')
theta,rho,xc,yc=ll.HoughTransform(data)
data['theta']=theta
data['rho']=rho

dtheta=2 #degrees
drho=500 #meters

dikeset, Z=ll.AggCluster(data)
lines,evaluation=examineCluster(data)
fig,ax=DotsLines(lines, ColorBy='AvgTheta')

```

Follow this indepth [tutorial](DemoLinkingLines.md) to get started!

You are now ready to utilize the power of Hough Line Transform, Agglomerative Clustering, and custom plotting in your data analysis projects. If you have any questions or need further assistance, please refer to the detailed documentation or contact our support team.

Happy coding!

