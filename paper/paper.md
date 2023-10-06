---
title: 'LinkingLines: Using the Hough Transform to Cluster Line Segments and for
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
affiliations:
 - name: University of Oregon
   index: 1
date: 5 October 2023
bibliography: paper.bib
---

# Summary

# Statement of Need

The `linkinglines` Python package addresses the critical need for efficient and accurate line clustering and analysis in geospatial and image data processing in addition to adding feature extraction capabilities. As the volume of data continues to grow across various domains, including remote sensing, computer vision, and geographic information systems (GIS), there is an increasing demand for tools that can simplify the extraction and analysis of linear features such as dikes, fractures, roads, rivers, and infrastructure.

This package was originally developed to tackle the issue of mapped dike segments. Rugged terrain, vegetation cover, and a large area made it impossible to accurately map dikes in high density areas called dike swarms. The length, density, and structure of the dike swarm is important to considering how magma is transported and erupted. Scaling analysis indicated that for segments of widths of 10 m as observed by __, dikes could be 10s of kilometers long however were an order of magnitude lower.

Additionally, the complex overlapping structure of the dike swarm was difficult to analyze. We designed `linkinglines` to extract not only lines from line segments but also help analyze the mesoscale structure of the dike swarm. Using the unique properties of the Hough Transform we an extract several unique mesoscale structures within a group of lines.

The primary needs that the `linkinglines` package fulfills include:

1. **Dissected Line Extraction or Data Reduction**: In areas where land cover or data availabilty effects the complete mapping of linear features, this package can link together similarly oriented segments into lines. This is also an effective data reduction technique. The `linkinglines` package offers algorithms and utilities for identifying and extracting such features from complex data.

2. **Feature Extraction and Analysis**: The `linkinglines` package provides feature extraction capabilities, allowing users to compute essential metrics and statistics on extracted linear, radial or circumferential type features, which is valuable in various scientific and engineering applications.

3. **Geospatial and Image Data Integration**: Geospatial data often involves complex relationships between various linear features. `linkinglines` integrates seamlessly with popular geospatial libraries like `pyproj` and `matplotlib` to facilitate georeferenced data analysis and visualization.

4. **Custom Plotting and Visualization**: Effective visualization is critical for data interpretation. The package offers custom plotting scripts, making it easier to visualize results and communicate findings effectively.

5. **Cross-Domain Applicability**: The package is versatile and adaptable to multiple domains, including geospatial analysis, image processing, and infrastructure monitoring, making it suitable for researchers, data scientists, and analysts across various disciplines.

In summary, the `linkinglines` Python package addresses a growing need for advanced line analysis and clustering tools in data-rich environments. Its capabilities empower users to efficiently process and analyze linear features in geospatial and image data, facilitating meaningful insights and supporting informed decision-making in a wide range of applications.

# Clustering and Feature Extraction

# Acknowledgements

We acknowledge contributions from . This paper was made possible by the "Crafting Quality Research Software and Navigating Publication in Software Journals" held by the Computational Infrastructure for Geodynamics in September 2023.

# References
