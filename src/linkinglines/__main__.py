# LinkingLines Package
# Written by aikubo
# __main__.py

import sys


def main():
    """
    LinkingLines - Command-line Utility for Line Data Processing

    Usage:
        python -m linkinglines [options]

    Options:
        -h, --help         Show this help message and exit.
        -v, --version      Show the version of the LinkingLines module.

    Description:
        The LinkingLines module is a command-line utility for processing line data stored in CSV files.
        It provides various features for analyzing and visualizing line data, including:

        - Hough Line Transform: Detect and extract lines from input data.
        - Agglomerative Clustering: Cluster similar lines for further analysis.
        - Custom Plotting: Generate visualizations of the clustered lines.
        - Feature Extraction: Extract line features for downstream analysis.
        - Export Results: Save the processed data and features to a CSV file.

    Examples:
        1. Perform Hough Line Transform and Agglomerative Clustering:
           python -m linkinglines process -i input_data.csv

        2. Generate custom plots of clustered lines:
           python -m linkinglines plot_clusters -i clustered_data.csv

        3. Extract line features from clustered data:
           python -m linkinglines extract_features -i clustered_data.csv

        4. Export results to a CSV file:
           python -m linkinglines export_results -i clustered_data.csv -o output_results.csv

        For detailed usage instructions and examples, refer to the documentation.

    Author:
        aikubo
        akubo@uoregon.edu
    """
    # Your main script logic goes here
    print(main.__doc__)

if __name__ == '__main__':
    main()
