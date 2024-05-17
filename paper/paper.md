---
title: 'LinkingLines: Using the Hough Transform to Cluster Line Segments and
Mesoscale Feature Extraction'
tags:
  - python
  - hough transform
  - earth science
  - clustering

authors:
  - name: Allison Kubo Hutchison
    orcid: 0000-0002-1378-361X
    affiliation: 1
  - name: Leif Karlstrom
    orcid: 0000-0002-2197-2349
    affiliation: 1
  - name: Tushar Mittal
    orcid: 0000-0002-8026-0018
    affiliation: 2
affiliations:
 - name: University of Oregon
 - name: Department of Earth Sciences, University of Oregon, Eugene, OR, USA
   index: 1
 - name: Department of Geosciences, Pennsylvania State University, University Park, PA, USA
   index: 2
date: 5 October 2023
bibliography: paper.bib
header-includes:
  - \usepackage{amsmath}
---

# Summary

Linear feature analysis plays a fundamental role in geospatial applications, from detecting infrastructure networks to characterizing geological formations. In this paper, we introduce `linkinglines`, an open-source Python package tailored for the clustering, and feature extraction of linear structures in geospatial data. Our package leverages the Hough Transform, commonly used in image processing, performs clustering of line segments in the Hough Space then provides unique feature extraction methods and visualization. `linkinglines` empowers researchers, data scientists, and analysts across diverse domains to efficiently process, understand, and extract valuable insights from linear features, contributing to more informed decision-making and enhanced data-driven exploration. We have used `linkinglines` to map dike swarms with thousands of segments associated with Large Igneous Provinces in @kubo2023.

# Statement of Need

The `linkinglines` Python package addresses the need for quantitative and automated line clustering and analysis in geospatial data. In geology and related fields, there is a need to analyze and extract patterns in overlapping and sometimes incomplete datasets of linear features such as fractures, dikes, or roads.

The primary needs that the `linkinglines` package fulfills include:

1. **Dissected Line Extraction or Data Reduction**: In areas where land cover or data availability affects the complete mapping of linear features, this package can link together similarly oriented segments into lines in a quantitative automated way. This is also an data reduction technique.

2. **Feature Extraction and Analysis**: The `linkinglines` package provides functions to extract and then compute essential metrics and statistics on extracted linear, radial or circumferential type features which are commonly seen in dike swarms or fracture networks.

3. **Custom Plotting and Visualization**: Effective visualization is critical for data interpretation. The package offers custom plotting scripts, making it easier to visualize results and communicate findings.

## State of the Field

This package was originally developed to tackle the issue of mapping dike segments. Rugged terrain, vegetation cover, and a large area made it impossible to accurately map dikes. The length, density, and structure of the dike swarm affects how magma is transported and erupted. Scaling analysis indicated that for segments of widths of 10 m, dikes could be 10s to 100s of kilometers long; however observed segments were two orders of magnitude lower [@morriss2020]. Additionally, the complex overlapping structure of the dike swarm was difficult to analyze. We designed `linkinglines` to extract not only lines from line segments but also help analyze the mesoscale structure of the dike swarm. Using the unique properties of the Hough Transform we can extract several unique mesoscale structures within a group of lines.

Issues with data coverage, incompleteness, and complexity occur in a wide range of Earth and Planetary science data products. As an example we show how this work can also be applied to lineaments on Venus[@venus] and fracture networks[@fractures] in our documentation. Clustering and machine learning methods are commonly applied geosciences within the subdisciplines of seismology, geochemistry, and planetary science [@li2023machine]. There is currently available code to work with geospatial data clustering such as `geopandas`, `PySAL`, and others which perform the Hough Transform on images such as `scikit-image` or `OpenCV`, however this seems to be the first union of the two methods and specific to geoscience data applications [@gpd; @scikit; @opencv; @pysal]. Related packages which are specific to geosciences such as `fractopo` which focuses on network analysis of fractures could be utilized in parallel with `linkinglines` for additional analysis [@fractopo].

# Algorithm

`linkinglines` can read in any common geospatial data format using `geopandas` such as ESRI Shapefiles, GeoJSON, or Well Known Text [@gpd]. Preprocessing is applied to the dataset so that only straight lines are considered in the algorithm. This is done by loading in the points from each vector object and performing linear regression, only considering those that yield a line with a $p>0.05$. The data is then formatted into a `pandas` DataFrame. This package heavily uses `pandas` as the database structure for ease of use, data manipulation, and integration with `numpy` and `scipy`[@pandas].

![Dike linking algorithm using the Hough Transform. First, raw data in Cartesian space is converted into Hough space (a and b). Agglomerative clustering is then performed on the data in Hough coordinates (d). In this example, there are four dikes total and two (red and blue) clusters. The clusters are redrawn by connecting the endpoints of the segments in the cluster (c).](houghexamplefig1.png)

The Hough Transform is an image processing technique used for detecting straight lines and other patterns in images[@hough; @ballard1981generalizing]. After loading in the data, it is assumed to be already line structures, so the accumulator array of the Hough Transform is skipped. First, the angle of the line segment is found. In many methods, it is the left-hand corner of the image [@ballard1981generalizing], but we choose the average midpoint of the line segments (Figure 1B). Other origins can be specified in certain functions using the `xc` and `yc` arguments. After the coordinate transform, $\rho$ and $\theta$ become the basis for the next step, where we utilize Scipy's clustering algorithm [@scipy].

After labels are assigned in the clustering portion of the algorithm, new lines are drawn using the endpoints of the clustered lines (Figure 1D). For each cluster, the nearest neighbors of the segment midpoints are calculated in Cartesian space, allowing for an analysis of the Cartesian spatial clustering of the lines. We also introduce an optional filtering step, which analyzes the maximum nearest neighbor difference of midpoints normalized by the total cluster length. This filters out clusters with segments that are not evenly clustered in Cartesian space.

## Feature Extraction

Additionally, we leverage the unique properties of the Hough Transform to combine clustering with feature extraction. In the original usage case of overlapping complex dike swarms, two potential end members of swarm types are linear and radial or circumferential swarms (Figure 2). Using the equation of the Hough Transform, we can apply curve-fitting to quantitatively fit the data to a radial pattern.

![Synthetic dike swarms in a Cartesian space (gray background, uppercase label) and Hough Transform space (white background, lowercase label). (A) Shows a simple linear swarm oriented at 30°. (B) Shows three linear swarms at −30°, 30°, 75°. (C) Shows three radial swarms aligned at a −45° angle. The angle at which radial swarms intersect in the Hough space (HS) is the angle of their relative orientation in Cartesian space. (D) Shows a circumferential swarm with the lines extending to show how it converges to Equation 8. The radius of the circumferential swarm is equal to the spacing of the parallel two curves in HS. (E) Shows three circumferential swarms with the same radius aligned at a -45° angle.](SyntheticsInkscape.png)

For extraction of linear features we can apply the Hough accumulator array, a 2D histogram of $\theta$ and $\rho$. You set the size of bins in the histogram and if clusters fall within those boxes they can be thought of as mesoscale clusters. We allow for flexibility of cutoffs for these mesoscale feature extraction, so it can be tailored to each research or engineering application.


# Future Work

This package takes geospatial or other types of line segment data and clusters them based on their orientation. Currently, the implementation only works with linear features however it could be generalized to arbitrary shapes for more flexibility [@ballard1981generalizing]. Additionally, future work could incorporate other shapes or patterns in the Hough Space and could extend the feature extraction methods laid out here. We invite collaboration to increase the capabilities of this code.

# Acknowledgements

This paper was made possible by the "Crafting Quality Research Software and Navigating Publication in Software Journals" held by the Computational Infrastructure for Geodynamics in September 2023.

# References
