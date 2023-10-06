.. linkinglines documentation master file, created by
   sphinx-quickstart on Wed Oct  4 11:52:32 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive

linkinglines: hough transform for clustering and feature extraction
===================================================================

.. image:: https://zenodo.org/badge/272334230.svg
   :target: https://zenodo.org/badge/latestdoi/272334230
   :alt: DOI


.. image:: https://img.shields.io/pypi/v/LinkingLines.svg
   :target: https://pypi.org/project/LinkingLines/
   :alt: PyPI


.. image:: https://readthedocs.org/projects/LinkingLines/badge/
   :target: https://linkinglines.readthedocs.io/
   :alt: ReadtheDocs


.. image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://opensource.org/licenses/MIT
   :alt: License


*Introduction*
    Welcome to the documentation for our Python module that performs the Hough
    Transform on line data from a CSV, clusters it using Agglomerative Clustering,
    and provides functionality to export the results into a CSV file.
    This module also includes custom plotting scripts and feature extraction
    methods to help you analyze and visualize your data effectively.

    This code was used to create the results published in
    `Kubo Hutchison et al., 2023 <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2022GC010842>`_.
    Initially, it was designed to link together mapped dike segments in Cartesian space
    to find their true lengths. This code can be applied to any linear features including
    roads, fractures, and other types of linear data.


    *
      **Data Clustering**\ : Apply Agglomerative Clustering to group similar data points, this
      can be used for data reduction, analysis, and mapping .

    *
      **Data Visualization**\ : Custom plotting scripts help you visualize and analyze
      your data, making it easier to identify patterns and anomalies.

    *
      **Feature Extraction**\ : Extract meaningful features from clustered data to
      perform further analysis, such as linear or radial type features.

1. Installation
---------------
To use this module, make sure you have Python installed (preferably Python 3.x).
    You can install the required packages using pip:

    .. code-block:: bash

       pip install linkinglines

2. Quick Start
--------------

    .. code-block:: python

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

    Follow this indepth `tutorial <>`_ to get started!

    You are now ready to utilize the power of Hough Line Transform, Agglomerative Clustering, and custom plotting in your data analysis projects. If you have any questions or need further assistance, please refer to the detailed documentation or contact our support team.

    Happy coding!

Indepth Tutorials
-----------------
.. toctree::
   maxdepth:2 

   DemoLinkingLines

Module Documentation
--------------------

.. toctree::
   :maxdepth: 4
   :caption: Contents:

   modules 

