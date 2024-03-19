
# Welcome to LinkingLines!
 [![status](https://joss.theoj.org/papers/64eeef828a1100bfba74052d89314758/status.svg)](https://joss.theoj.org/papers/64eeef828a1100bfba74052d89314758) [![DOI](https://zenodo.org/badge/272334230.svg)](https://zenodo.org/badge/latestdoi/272334230) [![PyPI](https://img.shields.io/pypi/v/LinkingLines.svg)](https://pypi.org/project/LinkingLines/) [![ReadtheDocs](https://readthedocs.org/projects/linkinglines/badge/)](https://linkinglines.readthedocs.io/) [![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Read the Full documentation on [ReadtheDocs!](https://linkinglines.readthedocs.io/en/latest/)

## 1. Introduction
Welcome to the documentation for our Python module that performs the Hough
Transform on lines from geospatial data, clusters it using Agglomerative Clustering.
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

data=ll.readFile('path/to/data')

dtheta=2 #degrees
drho=500 #meters

dikeset, Z=ll.AggCluster(data)
lines,evaluation=examineCluster(data)
fig,ax=DotsLines(lines, ColorBy='AvgTheta')

```

We have three examples:
1. Indepth tutorial with hough transform, clustering, and feature extraction using Spanish Peaks Data CSV file.
2. Hough Transform and feature extraction on Venus lineament data shape file.
3. Hough transform on fracture data geoJSON.

Data from: 
1. https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2022GC010842
2. https://pubs.usgs.gov/sim/3121/
3. https://doi.org/10.5281/zenodo.7919843

You are now ready to utilize the power of Hough Line Transform, Agglomerative Clustering, and custom plotting in your data analysis projects. If you have any questions or need further assistance, please refer to the detailed documentation or contact our support team.


## 4. Contributing Guidelines
Thank you for your interest in contributing to `linkinglines`. Please feel free to open up issues with bugs or requested features. Any contributions you make will benefit everybody else and are greatly appreciated.


If you would like to contribute code please do so in a seperate branch and open up an issue describing your contribution.

```
git clone git@github.com:USER/LinkingLines.git
git checkout my-development-branch
```

We recommend using a virutal environment to manage packages. We use `poetry` to manage dependencies and building.
See more about [poetry](https://python-poetry.org/).


```
pipx install poetry # you may need to install pipx first

cd linkinglines # go to the repo

poetry install --with test,dev # install in editable mode

# add your code 

poetry run pytest  # test code locally

```


Before submitting your pull request please verify the following:

1. Code is documented in [NumPy Docstring Style](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html)
2. Code is tested and passes test 
    - To run the tests please go to "/tests" and run `poetry run pytest` or `pytest`
    - Add your test code to any file with the name `test`
    - More here on [pytest and testing practices](https://docs.pytest.org/en/8.0.x/)
3. Open an issue and pull request 
4. After your pull request the code will be reviewed by maintainers. 
5. After passing review and automated tests it will be added to the next release and published to pypi.
