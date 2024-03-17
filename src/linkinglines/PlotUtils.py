# LinkingLines Package
 # Written by aikubo
 # Version: 2.1.0
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 12:49:51 2021

@author: akh
"""

import numpy as np
import math
import pandas as pd
from scipy import *
from matplotlib import cm
#from skimage.transform import probabilistic_hough_line as ProbHough
#from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import PercentFormatter

from pyproj import Proj
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, Normalize
import seaborn as sns
import matplotlib.colors as mcolors
from .HT import HT_center
from .FitRectangle import *
from .PrePostProcess import whichForm, FilterLines, getCartLimits
import scipy.cluster.hierarchy as sch
import labellines
import matplotlib.gridspec as gridspec
import string
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.colors import LogNorm

np.random.seed(5)

import matplotlib.ticker as mticker
class FixCartesianLabels:
    """
    Class for adjusting axis labels by moving the offset of axis tick labels to the axis label.

    This class is designed to assist in adjusting the display of axis labels in plots by incorporating the
    offset (scale factor) of the axis tick labels into the axis label itself. This adjustment is particularly
    useful in cases where the offset notation is preferred to be part of the axis label, enhancing the
    readability and presentation of plots.

    Parameters
    ----------
    ax : matplotlib.axis.Axis
        The axis object for which the labels are to be adjusted.

    Methods
    -------
    update(ax, lim)
        Updates the axis labels by incorporating the offset into the axis label text.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> ax.plot([1, 2, 3], [1e6, 2e6, 3e6])
    >>> y_label_fixer = FixCartesianLabels(ax.yaxis)
    >>> plt.show()

    Note
    ----
    The `update` method is connected to the axis's 'Cartesian Plots Updated' event and is triggered
    whenever the axis is updated, ensuring that the label adjustments are applied automatically.
    """

    def __init__(self, ax):
        """
        Initializes the FixCartesianLabels instance for the specified axis.

        Parameters
        ----------
        ax : matplotlib.axis.Axis
            The axis for which the labels should be updated.
        """
        ax.callbacks.connect('Cartesian Plots Updated', self.update)
        ax.figure.canvas.draw()
        self.update(ax, None)

    def update(self, ax, lim):
        """
        Updates the axis labels by moving the offset to the axis label.

        This method adjusts the axis labels by incorporating the offset (scale factor) from the axis tick labels
        into the axis label itself. The update is applied to both the x and y axes.

        Parameters
        ----------
        ax : matplotlib.axis.Axis
            The axis for which the labels are to be updated.
        lim : unused
            An unused parameter, included for compatibility with callback requirements.

        Returns
        -------
        None
        """
        for i, l in zip([ax.yaxis, ax.xaxis], ['Y', 'X']):
            fmt = i.get_major_formatter()
            i.offsetText.set_visible(False)
            i.set_label_text(l + " (" + fmt.get_offset() + " m)")


def get_ax_size_inches(ax):
    """
    Calculate the size of a matplotlib axis in inches.

    This function determines the size (width and height) of a specified matplotlib axis in inches,
    taking into account the current figure's DPI settings. It's useful for precise layout adjustments
    or when needing to scale other plot elements relative to the axis size.

    Parameters
    ----------
    ax : matplotlib.axes._subplots.AxesSubplot
        The axis for which to calculate the size.

    Returns
    -------
    width : float
        The width of the axis in inches.
    height : float
        The height of the axis in inches.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> width, height = get_ax_size_inches(ax)
    >>> print(f"Width: {width} inches, Height: {height} inches")
    """
    fig = ax.get_figure()
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    return bbox.width, bbox.height


def FixAxisAspect(ax1, ax2):
    """
    Adjust the aspect ratios of two matplotlib axes to match each other.

    This function aligns the aspect ratios of two provided axes, ensuring that the graphical representation
    in both axes appears with consistent scaling. It's particularly useful in comparative visualizations where
    matching scales and proportions are crucial.

    Parameters
    ----------
    ax1 : matplotlib.axes._subplots.AxesSubplot
        The first axis, to which the second axis's aspect ratio will be adjusted.
    ax2 : matplotlib.axes._subplots.AxesSubplot
        The second axis, which will be adjusted to match the aspect ratio of the first axis.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> fig, ax1 = plt.subplots()
    >>> ax2 = plt.axes([0.2, 0.2, 0.4, 0.4])
    >>> FixAxisAspect(ax1, ax2)
    >>> plt.show()

    Note
    ----
    The aspect ratio adjustment is based on the size of the axes in inches and their positioning within the figure.
    This method may alter the ylim or xlim of `ax2` to ensure that the aspect ratios match.
    """
    # Implementation omitted for brevity.


def labelSubplots(ax, labels=None, **kwargs):
    """
    Add alphabetical labels to the corner of each subplot.

    This function labels each subplot in a figure with an alphabetical character, starting from "A". It can
    enhance the clarity of figures containing multiple subplots by providing a simple reference system.

    Parameters
    ----------
    ax : list, dict, or ndarray
        A collection of AxesSubplot objects to be labeled. If a single AxesSubplot is provided, it will be
        converted into a list automatically.
    labels : list of str, optional
        Custom labels to use instead of the default alphabetical labels. The list must have the same length
        as the `ax` collection. If None, labels will default to "A", "B", "C", etc.
    **kwargs : dict
        Additional keyword arguments to be passed to the label positioning function.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> fig, axs = plt.subplots(2, 2)
    >>> labelSubplots(axs.flatten())
    >>> plt.show()

    Note
    ----
    The function supports lists, dictionaries, or arrays of AxesSubplot objects. In the case of a single
    AxesSubplot, it will be handled appropriately by converting to a list format internally.
    """


    if type(ax) is list or len(ax) > 1:
        if labels is None:
            labels=list(string.ascii_lowercase)[:len(ax)]
            ax_dict=dict(zip(labels,ax))
        else:
            ax_dict=dict(zip(labels,ax))
        identify_axes(ax_dict, **kwargs)
    elif type(ax) is dict:
        identify_axes(ax, **kwargs)
    else:
        ax_list=ax.tolist()
        if labels is None:
            labels=list(string.ascii_lowercase)[:len(ax)]
            ax_dict=dict(zip(labels,ax))
        else:
            ax_dict=dict(zip(labels,ax))
        identify_axes(ax_dict, **kwargs)





def identify_axes(ax_dict, fontsize=12, **kwargs):
    """
    Helper to identify the Axes in the examples below.

    Draws the label in a large font in the center of the Axes.

    Parameters
    ----------
    ax_dict : dict[str, Axes]
        Mapping between the title / label and the Axes.
    fontsize : int, optional
        How big the label should be
    **kwargs : dict
        Additional arguments are passed to `ax.text`.
    
    Returns
    -------
    None

    """
    kw = dict(ha="center", va="center", fontsize=fontsize, fontstyle="oblique", color="black", **kwargs)
    for k, ax in ax_dict.items():
        ax.text(0.95, 0.05, k, transform=ax.transAxes, **kw)

from operator import sub

def get_aspect(ax):
    """
    Calculate the aspect ratio of a matplotlib axes.

    This function calculates the aspect ratio of a given matplotlib axes, taking into account both the aspect ratio of the figure and the aspect ratio of the data displayed in the axes.

    Parameters:
    ----------
    ax : matplotlib.axes._subplots.AxesSubplot
        The axes for which to calculate the aspect ratio.

    Returns:
    -------
    float
        The aspect ratio of the axes.

    Example:
        import matplotlib.pyplot as plt

        # Create a sample plot
        fig, ax = plt.subplots()

        # Calculate and print the aspect ratio
        aspect_ratio = get_aspect(ax)
        print("Aspect Ratio:", aspect_ratio)
    """
    # Total figure size
    figW, figH = ax.get_figure().get_size_inches()
    # Axis size on figure
    _, _, w, h = ax.get_position().bounds
    # Ratio of display units
    disp_ratio = (figH * h) / (figW * w)
    # Ratio of data units
    # Negative over negative because of the order of subtraction
    data_ratio = sub(*ax.get_ylim()) / sub(*ax.get_xlim())

    return disp_ratio / data_ratio


def RGBtoHex(vals, rgbtype=1):
  """Converts RGB values in a variety of formats to Hex values.

    Parameters:
    ----------
    vals : tuple
        An RGB/RGBA tuple
    rgbtype : int, default 1
        Valid valus are:
                    1 - Inputs are in the range 0 to 1
                    256 - Inputs are in the range 0 to 255

    Returns:
    -------
    list
         A hex string in the form '#RRGGBB' or '#RRGGBBAA'
    """

  if len(vals)!=3 and len(vals)!=4: 
    raise Exception("RGB or RGBA inputs to RGBtoHex must have three or four elements!")
  if rgbtype!=1 and rgbtype!=256:
    raise Exception("rgbtype must be 1 or 256!")

  #Convert from 0-1 RGB/RGBA to 0-255 RGB/RGBA
  if rgbtype==1:
    vals = [255*x for x in vals]

  #Ensure values are rounded integers, convert to hex, and concatenate
  return '#' + ''.join(['{:02X}'.format(int(round(x))) for x in vals])

def RGBArraytoHexArray(c):
    """
    Convert an array of RGB or RGBA values to an array of Hex values.

    This function takes an array of RGB or RGBA tuples and converts them to an array of their corresponding hexadecimal color representations.

    Parameters:
    ----------
    c : list 
        A list of RGB or RGBA tuples.

    Returns:
    -------
    list
    A list of hex strings in the form '#RRGGBB' or '#RRGGBBAA'.

    Example:
    -------
    >>> # Convert an array of RGB (0-1 range) to Hex
    >>> rgb_colors = [(0.2, 0.5, 0.8), (0.9, 0.1, 0.4)]
    >>> hex_colors = RGBArraytoHexArray(rgb_colors)
    >>> print("Hex Colors:", hex_colors)
    """
    return [RGBtoHex(i) for i in c]


def StringColors(values, palette="turbo"):
    """
    Map a list of strings to colors using a specified color palette.

    This function takes a list of strings and maps each unique string to a unique color from a specified color palette.

    Parameters:
    ----------
    values : list
        A list of strings to be mapped to colors.
    palette : str, default "turbo"
        The name of the color palette to use.

    Returns:
    -------
    tuple
        A tuple containing two elements:
        - color_idx (numpy.ndarray): An array of indices representing the colors for each string.
        - cm (matplotlib.colors.LinearSegmentedColormap): The colormap used for mapping the strings to colors.

    Example:
    -------
    >>> # Map a list of categories to colors
    >>> categories = ["Category A", "Category B", "Category C", "Category A"]
    >>> color_indices, colormap = StringColors(categories, palette="viridis")
    >>> print("Color Indices:", color_indices)
    """
    if type(values[0]) is not str:
        raise ValueError("Must pass a list of strings.")
    if len(values) < 2:
        raise ValueError("Must pass a list with more than 2 elements.")

    labels = np.unique(values)
    n_colors = len(labels)

    colors = [RGBtoHex(x) for x in sns.color_palette(palette)]

    cm = LinearSegmentedColormap.from_list("StringCM", colors, N=n_colors)
    color_idx = np.array([np.where(i == labels)[0][0] for i in values])

    return color_idx, cm

def StringCbar(c, fig, ax, values):
def StringCbar(c, fig, ax, values):
    """
    Create a colorbar for string-based categorical data.

    This function generates a colorbar specifically for visualizing categorical data that has been
    encoded with colors through a mapping process. It takes as input the collection of plot elements
    colored according to their category, the figure and axis objects containing the plot, and a list
    of the categorical values. The function then produces a colorbar that reflects this categorical
    color mapping, facilitating the interpretation of the plot's color scheme in relation to the
    categorical data.

    Parameters
    ----------
    c : matplotlib.collections.Collection
        The collection of colored elements in the plot, such as scatter plot points, which have been
        colored based on categorical data.
    fig : matplotlib.figure.Figure
        The figure object containing the plot. This is used to ensure that the colorbar is correctly
        sized and positioned within the plot.
    ax : matplotlib.axes.Axes
        The axes object on which the plot and the colorbar will be drawn. This specifies the plotting
        area for both the main plot and the associated colorbar.
    values : list of str
        A list of strings representing the categories in the categorical data. Each category is associated
        with a specific color in the plot.

    Returns
    -------
    colorbar : matplotlib.colorbar.Colorbar
        The colorbar object created for the categorical data. This colorbar shows the color associated
        with each category, facilitating the interpretation of the plot.

    Examples
    --------
    >>> # Assuming StringColors function has been used to map 'values' to colors
    >>> values = ["Category A", "Category B", "Category C", "Category A"]
    >>> color_indices, colormap = StringColors(values, palette="viridis")
    >>> scatter = ax.scatter(x, y, c=color_indices, cmap=colormap)
    >>> colorbar = StringCbar(scatter, fig, ax, values)
    >>> plt.show()

    Note
    ----
    This function is designed to work with plots where categorical data has been mapped to colors
    using a specific encoding function (e.g., `StringColors`). It assumes that such a mapping
    process has already been applied to create the input collection `c`.
    """

    labels = np.unique(values)
    n_colors = len(labels)
    c_ticks = np.arange(n_colors)

    # Clean labels by removing any prefixes (e.g., "Category:")
    for i in range(len(labels)):
        if ":" in labels[i]:
            labels[i] = labels[i].split(":")[1]

    # Create the colorbar and set tick labels
    cbar = fig.colorbar(c, ticks=c_ticks, ax=ax)
    cbar.ax.set_yticklabels(labels, rotation='vertical')

    return cbar


def fontItems(fig, ax):
    """
    Generate a list of items with fonts to change in a matplotlib figure.

    This function iterates over a given figure and axes (or a list of axes) to compile a list of all text items
    whose font properties can be modified. This includes titles, axis labels, and tick labels for each axis provided.
    Additionally, if the figure contains a legend, the fonts of the legend's texts are also included.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The matplotlib figure object.
    ax : matplotlib.axes.Axes or list of matplotlib.axes.Axes
        Single axis object or a list of axis objects contained within `fig`.

    Returns
    -------
    fontItems : list
        A list of matplotlib text objects for which the font properties can be changed. This includes titles,
        axis labels, tick labels, and legend texts within the figure and the specified axes.

    """
    if isinstance(ax, list):
        items = []
        for a in ax:
            items.append([a.title, a.xaxis.label, a.yaxis.label] + a.get_xticklabels() + a.get_yticklabels())
        fontItems = [item for sublist in items for item in sublist]
    else:
        fontItems = [ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()

    legend = fig.get_legend()
    if legend is not None:
        fontItems += legend.get_texts()
    return fontItems

def jgrSize(fig, ax, finalSize, units="mm"):
    """
    Resize a matplotlib figure to a Journal of Geophysical Research (JGR) appropriate size and set dpi to 600.

    This function adjusts the size of a matplotlib figure to conform to the specifications for JGR submissions,
    such as "quarter", "half", "full", or a specific size given by the user. It also sets the figure's DPI to 600,
    ensuring high-resolution output suitable for publication.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The matplotlib figure object to be resized.
    ax : matplotlib.axes.Axes or list of matplotlib.axes.Axes
        The axes contained within `fig` which will be adjusted alongside the figure resizing.
    finalSize : str or float or tuple or list
        The target size for the figure. Can be "quarter", "half", "full", a float representing a fraction of the full size,
        or a tuple/list specifying the width and height in the chosen units.
    units : str, optional
        The units of measurement for specifying custom sizes, by default "mm". Other valid options are "cm" and "inches".

    Returns
    -------
    fig : matplotlib.figure.Figure
        The resized matplotlib figure object.

    """
    # Resize first
    [wo,lo]=fig.get_size_inches()

    if wo > lo:
        orientation='landscape'
    elif wo < lo:
        orientation='portrait'
    else:
        orientation='landscape'

    if finalSize=='quarter':
        w=95/25.4 ## cm/(cm/in)
        l=115/25.4
    elif finalSize=='half':
        w=190/2/25.4
        l=230/2/25.4
    elif finalSize=='full':
        w=190/25.4
        l=230/25.4
    elif type(finalSize)==float:
        w=190/25.4*finalSize
        l=230/25.4*finalSize
    elif type(finalSize)==list or type(finalSize)==tuple:
        if units=='mm':
            w=finalSize[0]/25.4
            l=finalSize[1]/25.4
        elif units=='cm':
            w=finalSize[0]/2.54
            l=finalSize[1]/2.54
        elif units=='inches':
            w=finalSize[0]
            l=finalSize[1]
    else:
        w=190/2/25.4
        l=230/2/25.4

    if orientation=='landscape':
        fig.set_size_inches(l,w, forward=True)
    if orientation=='portrait':
        fig.set_size_inches(w,1, forward=True)


    #change font sizes
    items=fontItems(fig, ax)
    #default size is close to 8.5
    #there's no way to seperately set the super script and subscripts
    #it uses the default TeX settings of superscript = 70% of main font
    #so to maintain a MINIMUM of 6 pts in the superscripts
    #fonts must be set to 8.5
    DEFAULT_SIZE=6/0.7

    for i in items:
        i.set_fontsize(DEFAULT_SIZE)


    fig.set_dpi(600)


    return fig

def SetupAGUFig(finalSize, orientation, units='mm'):
    """
    Set up a Matplotlib figure for creating a American Geophysical Union Journal style plot.

    This function initializes a Matplotlib figure with specific dimensions and font settings suitable for
    JGR submission requirements. The figure size can be specified as a standard page size (quarter, half, full)
    or in custom dimensions. The function also adjusts global font sizes to ensure readability at the specified
    figure size.

    Parameters
    ----------
    finalSize : str or float or tuple or list
        The desired size of the figure. Can be 'quarter', 'half', 'full' for standard JGR sizes, a float
        indicating a fraction of the full page size, or a tuple/list specifying custom dimensions (width, height).
    orientation : str
        The orientation of the figure, either 'landscape' or 'portrait'.
    units : str, optional
        The units for the custom dimensions ('mm', 'cm', 'inches'), by default 'mm'.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The initialized Matplotlib figure object with the specified dimensions and font settings.

    Examples
    --------
    >>> fig = SetupJGRFig('half', 'landscape')
    >>> plt.plot(x, y)
    >>> plt.xlabel('X Label')
    >>> plt.ylabel('Y Label')
    >>> plt.title('AGU-style Plot')
    >>> plt.savefig('AGU_plot.png', dpi=300, bbox_inches='tight')
    >>> plt.show()

    Notes
    -----
    - This function adjusts the global matplotlib font settings which may affect other plots. Consider
      using matplotlib's context manager (`with plt.rc_context()`) if this is a concern.
    - The specified figure size and orientation are intended to match AGU submission standards, ensuring
      that figures fit well within the page layout of the journal.
    """
    # Set font sizes to ensure readability at JGR publication standards
    plt.rc('font', size=8)  # sets default font size
    plt.rc('axes', titlesize=8)  # fontsize of the axes title
    plt.rc('axes', labelsize=8)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=8)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=8)  # fontsize of the tick labels
    plt.rc('legend', fontsize=8)  # legend fontsize
    plt.rc('figure', titlesize=8)  # fontsize of the figure title
    plt.rcParams['legend.title_fontsize'] = 8  # fontsize of the legend title

    # Create the figure with the appropriate size
    fig = plt.figure()
    # Set figure size based on specified parameters
    if finalSize in ['quarter', 'half', 'full']:
        size_map = {
            'quarter': (95 / 25.4, 115 / 25.4),
            'half': (190 / 25.4 / 2, 230 / 25.4) if orientation == 'landscape' else (190 / 25.4, 230 / 25.4 / 2),
            'full': (190 / 25.4, 230 / 25.4),
        }
        w, l = size_map[finalSize]
    elif isinstance(finalSize, float):
        w, l = 190 / 25.4 * finalSize, 230 / 25.4 * finalSize
    elif isinstance(finalSize, (list, tuple)):
        unit_conversion = {'mm': 25.4, 'cm': 2.54, 'inches': 1}
        w, l = finalSize[0] / unit_conversion[units], finalSize[1] / unit_conversion[units]
    else:
        w, l = 190 / 2 / 25.4, 230 / 2 / 25.4

    # Apply orientation
    if orientation == 'landscape':
        fig.set_size_inches(l, w)
    elif orientation == 'portrait':
        fig.set_size_inches(w, l)

    return fig

def combinePlots(fig1,fig2, path):
    """
    Combine two figures into one and save

    Parameters
    ----------
    fig1 : matplotlib figure.Figure
        lefthand figure
    fig2 : matplotlib figure.Figure
        righthand figure
    path: string
        path to save combined figure

    Returns
    -------
    None.

    """
    #fig1=jgrSize(fig1, fig1.get_axes(), 'quarter')
    #fig2=jgrSize(fig2, fig2.get_axes(), 'quarter')


    backend = mpl.get_backend()
    mpl.use('agg')

    dpi = 600


    c1 = fig1.canvas
    c2 = fig2.canvas

    c1.draw()
    c2.draw()

    a1 = np.array(c1.buffer_rgba())
    a2 = np.array(c2.buffer_rgba())
    a = np.hstack((a1,a2))

    mpl.use(backend)
    fig,ax = plt.subplots(dpi=dpi)
    fig.subplots_adjust(0, 0, 1, 1)
    ax.set_axis_off()
    ax.matshow(a)
    fig.savefig(path, dpi=dpi)



def clustered_lines(xs, ys, theta, length, xmid=None, ymid=None):
    """
    Calculate the coordinates of two points to represent a line segment based on clustering.

    This function computes the coordinates of two points defining a line segment, given a set of x and y coordinates,
    a specified angle, and a length. The line segment is designed to have one end clustered around a central point,
    defined by (xmid, ymid). If the central point is not provided, it is calculated as the mean of the input coordinates.

    Parameters
    ----------
    xs : array-like
        An array-like object containing x-coordinates of data points.
    ys : array-like
        An array-like object containing y-coordinates of data points.
    theta : float
        The angle of the line segment relative to the horizontal, in degrees.
    length : float
        The length of the line segment.
    xmid : float, optional
        The x-coordinate of the central point around which one end of the line segment is clustered. If None, the
        mean of `xs` is used.
    ymid : float, optional
        The y-coordinate of the central point around which one end of the line segment is clustered. If None, the
        mean of `ys` is used.

    Returns
    -------
    tuple
        A tuple (x1, y1, x2, y2) representing the coordinates of the two endpoints of the line segment.

    Examples
    --------
    >>> xs = [1, 2, 3, 4, 5]
    >>> ys = [2, 3, 4, 5, 6]
    >>> theta = 45  # Angle in degrees
    >>> length = 3
    >>> xmid, ymid = 3, 4  # Central point
    >>> x1, x2, y1, y2 = clustered_lines(xs, ys, theta, length, xmid, ymid)
    >>> print(f'Point 1: ({x1}, {y1})')
    >>> print(f'Point 2: ({x2}, {y2})')

    Notes
    -----
    - The central point `(xmid, ymid)` serves as the midpoint of the line segment if not otherwise specified.
    """
 
    xstart = np.max(xs)
    ystart = np.max(ys)

    xend = np.min(xs)
    yend = np.min(ys)

    if xmid is None or ymid is None:
        print('Calculating xmid and ymid')
        xmid = (xstart + xend) / 2
        ymid = (ystart + yend) / 2

    a = np.cos(np.deg2rad(theta))
    b = np.sin(np.deg2rad(theta))

    x0 = xmid
    y0 = ymid
    x1 = int(x0 + length / 2 * (-b))
    y1 = int(y0 + length / 2 * (a))
    x2 = int(x0 - length / 2 * (-b))
    y2 = int(y0 - length / 2 * (a))

    return x1, x2, y1, y2
def pltRec(lines, xc, yc, fig=None, ax=None):
    """
    Plot a rectangle defined by the center and lines, illustrating the orientation and dimensions.

    This function creates a visualization of a rectangle that represents the orientation and dimensions
    of a set of lines with respect to a given center point. It optionally takes a matplotlib figure and
    axis for plotting or creates new ones if not provided. The rectangle is determined based on the
    aggregate properties (angles and lengths) of the lines DataFrame.

    Parameters
    ----------
    lines : pandas.DataFrame
        A DataFrame containing line data, with columns for angles ('theta') and lengths ('rho').
    xc : float
        The x-coordinate of the center point.
    yc : float
        The y-coordinate of the center point.
    fig : matplotlib.figure.Figure, optional
        The figure object for the plot. If None, a new figure is created.
    ax : matplotlib.axes._subplots.AxesSubplot, optional
        The axes object for the plot. If None, new axes are created on the provided or new figure.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object used for plotting.
    ax : matplotlib.axes._subplots.AxesSubplot
        The axes object used for plotting.
    length : float
        The calculated length of the rectangle.
    width : float
        The calculated width of the rectangle.

    Examples
    --------
    >>> # Example DataFrame of lines
    >>> import pandas as pd
    >>> lines_df = pd.DataFrame({'theta': [0, 90, 45], 'rho': [5, 5, 7]})
    >>> xc, yc = 0, 0  # Center point
    >>> fig, ax, length, width = pltRec(lines_df, xc, yc)
    >>> plt.show()

    Note
    ----
    The orientation of the rectangle is determined by the mean angle of the lines, and its dimensions are
    based on the mean and range of the lengths (rho) of the lines. The rectangle is plotted with its center
    at the given (xc, yc) point, rotated to match the average orientation of the lines.
    """
    col = lines.columns

    post = ['Theta', 'AvgTheta', 'theta']
    posr = ['Rho', 'AvgRho', 'rho']

    for p in post:
        if p in col:
            t = p

    for p in posr:
        if p in col:
            r = p

    if 'Average Rho (m)' in col:
        r = 'Average Rho (m)'
        t = r'Average Theta ($^\circ$)'

    if t == 'AvgTheta' or t == r'Average Theta ($^\circ$)':
        segl = 'R_Length'
    else:
        segl = 'seg_length'

    if fig is None or ax is None:
        fig, ax = plt.subplots()
    xi, yi = endpoints2(lines)
    x0 = xc
    y0 = yc
    size = len(lines)
    if abs(np.sum(np.sign(lines[t].values))) < size:
        crossZero = True
        ang = np.mean(abs(lines[t].values))
        tol = 6
        if np.isclose(ang, 0, atol=4):
            ang = np.mean((lines[t].values))
    else:
        crossZero = False
        ang = np.mean((lines[t].values))

    xp, yp = rotateXYShift(np.deg2rad(-1 * ang), xi, yi, x0, y0)

    width = np.ptp(xp.flatten())
    length = np.ptp(yp.flatten())

    xc = (max(xp) - min(xp)) / 2 + min(xp)
    yc = (max(yp) - min(yp)) / 2 + min(yp)

    xr = xc + width / 2
    xl = xc - width / 2
    yu = yc + length / 2
    yd = yc - length / 2
    xs = np.append(xr, xl)
    ys = np.append(yu, yd)

    xpi, ypi = unrotate(np.deg2rad(-1 * ang), xp, yp, x0, y0)

    Xedges = np.array([xs[0], xs[0], xs[1], xs[1], xs[0]])
    Yedges = np.array([ys[1], ys[0], ys[0], ys[1], ys[1]])

    Xmid = (np.max(xs) + np.min(xs)) / 2
    Ymid = (np.max(ys) + np.min(ys)) / 2

    xs, ys = unrotate(np.deg2rad(-1 * ang), Xedges, Yedges, x0, y0)
    Xmid, Ymid = unrotate(np.deg2rad(-1 * ang), Xmid, Ymid, x0, y0)

    ax.plot(xs, ys, 'k-.', alpha=0.7, linewidth=1)

    xstart, xend, ystart, yend = clustered_lines(xi, yi, ang, length, xmid=Xmid, ymid=Ymid)

    ax.plot([xstart, xend], [ystart, yend], 'g.-', linewidth=1)
    for i in range(0, len(lines)):
        ax.plot([xi[i], xi[i + len(lines)]], [yi[i], yi[i + len(lines)]], 'r-', linewidth=2)
    ax.plot(Xmid, Ymid, 'yp')
    ax.set_aspect('equal')

    return fig, ax, length, width


def labelcolors(labels, colormap):
    """
    Assigns colors to unique labels using a colormap.

    This function takes a list of labels and a colormap and assigns a unique color to each unique label based on the
    colormap. It returns a list of colors corresponding to the input labels.

    Parameters:
    ----------
    labels : list or pandas.Series
        A list of labels.
    colormap : matplotlib.colors.Colormap
        A colormap to assign colors from.

    Returns:
    -------
    colors : list
        A list of colors in hexadecimal format (#RRGGBB) corresponding to the input labels.
    colors_short : list
        A list of short color names (e.g., 'red', 'blue') corresponding to the input labels.

    Example:
    -------
    >>> # Assign colors to unique labels using a colormap
    >>> labels = ['A', 'B', 'A', 'C', 'B']
    >>> colormap = plt.get_cmap('viridis')
    >>> label_colors, short_colors = labelcolors(labels, colormap)
    """
    n = len(np.unique(labels))
    c = colormap(np.linspace(0, 1, n))
    colors = []
    colors_short = [RGBtoHex(c[i]) for i in range(n)]

    for i in range(len(labels)):
        c_loc = np.where(np.unique(labels) == labels.iloc[i])[0][0]
        c_val = RGBtoHex(c[c_loc])
        colors.append(c_val)

    return colors, colors_short


def plotlines(data, col, ax, alpha=1, myProj=None, maskar=None, linewidth=1,
              ColorBy=None, center=False, xc=None, yc=None, extend=False,
              cmap=cm.turbo, cbarStatus=False, SpeedUp=True, equal=True):
    """
    Plots line segments on a specified axis with optional attributes for color, transparency, and more.

    This function allows for the visualization of line segments from a DataFrame containing their start and end
    UTM coordinates. It offers various customization options including color coding by data attributes, transparency,
    line extension, and more.

    Parameters
    ----------
    data : pandas.DataFrame
        DataFrame containing 'Xstart', 'Xend', 'Yend', 'Ystart' columns representing UTM coordinates of line segments.
    col : str or RGB tuple
        Color for plotting the lines, specified as a color name string or an RGB tuple (e.g., (0.5, 0.5, 0.5)).
    ax : matplotlib.axes._axes.Axes
        Axes object for plotting.
    alpha : float, optional
        Transparency level of lines (0.0 transparent, 1.0 opaque). Default is 1.
    myProj : pyproj.Proj, optional
        PyProj projection object for converting UTM coordinates to lat/long. Default is None.
    maskar : array-like, optional
        Logical mask indicating which lines to plot. Default is None.
    linewidth : float, optional
        Width of the plotted lines. Default is 1.
    ColorBy : str, optional
        Column name from `data` to color lines based on its values. Default is None.
    center : bool, optional
        Whether to plot the center point. Default is False.
    xc : float, optional
        X-coordinate of the center point, required if `center` is True. Default is None.
    yc : float, optional
        Y-coordinate of the center point, required if `center` is True. Default is None.
    extend : bool, optional
        If True, extends lines beyond their endpoints. Default is False.
    cmap : matplotlib.colors.Colormap, optional
        Colormap for coloring lines based on `ColorBy`. Default is cm.turbo.
    cbarStatus : bool, optional
        If True, displays a colorbar for color-coded lines. Default is False.
    SpeedUp : bool, optional
        If True, downsamples data for faster plotting with large datasets. Default is True.
    equal : bool, optional
        If True, sets aspect ratio to equal, maintaining scale. Default is True.

    Returns
    -------
    None

    Example
    -------
    >>> plotlines(data, 'red', ax, alpha=0.6, ColorBy='Theta', cbarStatus=True)
    >>> # Plots line segments in red with alpha transparency, color lines based on 'Theta' column, and displays a colorbar.
    """

    if maskar is not None:
        temp=data.loc[maskar]
        if not isinstance(col,str):
            col=col[maskar]

    else :
        temp=data

    if SpeedUp and len(data)>2000:

        temp=data.sample(frac=0.25)
        print("Downsampled to", len(temp)," lines from", len(data))


    if len(data) > 15:
        alpha=0.8
    elif len(data) > 100:
        alpha=0.4

    if ColorBy is not None:
        C=data[ColorBy].values

        if type(C[0]) is str:
            C,cmap=StringColors(C)
            col=[RGBtoHex(a) for a in cmap(C)]
            print("in plotlines, colorby is a str")

        else:
            norm=Normalize(vmin=min(C), vmax=max(C))

            m = cm.ScalarMappable(norm=norm, cmap=cmap)
            col=m.to_rgba(C)




    for i in range(0,len(temp)):
        x1=temp['Xstart'].iloc[i]
        y1=temp['Ystart'].iloc[i]
        y2=temp['Yend'].iloc[i]
        x2=temp['Xend'].iloc[i]




        if myProj is not None:
            lon1, lat1 = myProj(x1, y1, inverse = True)
            lon2, lat2 = myProj(x2, y2, inverse = True)
        else:
            lon1=x1
            lat1=y1
            lon2=x2
            lat2=y2

        LAT = [lat1, lat2]
        LONG = [lon1, lon2]

        if extend:
            m=(lat1-lat2)/(lon1-lon2)
            b=lat1-lon1*m

            l=np.sqrt((lat1-lat2)**2+(lon1-lon2)**2)
            LONG=[lon1-l*4, lon2+l*4]
            LAT=np.multiply(LONG,m)+b

        #color handling

        if ColorBy is None:
            if isinstance(col,list):
                colo=col[i]
            elif isinstance(col, str):
                colo=col
            else:
                colo=col
        else:
            colo=col[i]


        ax.plot(LONG,LAT, c=colo, alpha=alpha, linewidth=linewidth)

    if center:
        if xc is None or yc is None:
            xc,yc=HT_center(data)
        ax.plot(xc,yc, "wp", markeredgecolor="black", markersize=10)

        if equal:

            FixCartesianLabels(ax)
            ax.set_aspect('equal')


def HThist(lines, rstep, tstep, weights=None, fig=None, ax=None, rbins=None, tbins=None, cmap=cm.Blues, gamma=0.3):
    """
    Plot a 2D histogram (Hough Transform space) of line data.

    This function creates a 2D histogram representing the Hough Transform space of a set of lines. It visualizes
    the distribution of lines' orientations (theta) and distances (rho) from the origin. The histogram is plotted
    with specified steps for rho and theta, and can optionally weight the lines. The plot can be customized with
    matplotlib's figure and axes objects, bin sizes, colormap, and a gamma correction for the colormap normalization.

    Parameters
    ----------
    lines : pandas.DataFrame or array-like
        Data containing line information. Must have columns or keys for 'rho' and 'theta', representing the distance
        and orientation of lines, respectively.
    rstep : float
        Step size for the rho (distance) dimension of the histogram.
    tstep : float
        Step size for the theta (orientation) dimension of the histogram.
    weights : array-like, optional
        Weights for each line, affecting how they contribute to the histogram. Default is None.
    fig : matplotlib.figure.Figure, optional
        Figure object on which to plot. If None, a new figure is created. Default is None.
    ax : matplotlib.axes.Axes, optional
        Axes object on which to plot. If None, new axes are created on the figure. Default is None.
    rbins : array-like, optional
        Bin edges for rho (distance) dimension. If None, bins are created based on rstep. Default is None.
    tbins : array-like, optional
        Bin edges for theta (orientation) dimension. If None, bins are created based on tstep. Default is None.
    cmap : matplotlib.colors.Colormap, optional
        Colormap for the histogram. Default is matplotlib.cm.Blues.
    gamma : float, optional
        Gamma correction for the colormap normalization. Default is 0.3.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object used for the plot.
    ax : matplotlib.axes.Axes
        The axes object used for the plot.
    histogram_data : list
        A list containing the histogram array (h), bin edges for x (xe), bin edges for y (ye), and the QuadMesh object (c).

    Example
    -------
    >>> lines = pd.DataFrame({'rho': np.random.rand(100) * 200, 'theta': np.random.rand(100) * 180 - 90})
    >>> fig, ax, hist_data = HThist(lines, 1, 10)
    >>> plt.show()

    """
    t,r=whichForm(lines)

    if fig is None or ax is None:
        fig,ax=plt.subplots()
    if rbins is None:
        rbins=np.arange(min(lines[r])/1000, max(lines[r])/1000, rstep)
    if tbins is None:
        tbins=np.arange(-90, 90, tstep)
    h,xe,ye, c=ax.hist2d(lines[t], lines[r]/1000, bins=[tbins, rbins], weights=weights, cmap=cmap, norm=mcolors.PowerNorm(gamma))
    fig.colorbar(c, label='Counts', ax=ax)
    ax.set_xlabel(r'Theta ($^\circ$)')
    ax.set_ylabel('Rho (km)')

    return fig, ax, [h, xe,ye,c]

def annotateWLines(ax, angles=None):
    """
    Annotate the given axis with lines at specified angles.

    This function adds lines to a matplotlib axis at specified angles to serve as annotations. It is designed to
    work when the axis aspect ratio is set to equal. The lines extend from a point just above the top of the axis,
    across the plotting area, at the specified angles.

    Parameters
    ----------
    ax : matplotlib.axes._axes.Axes
        The axes object to which the annotation lines will be added.
    angles : list of float, optional
        The angles in degrees at which to draw the lines. If None, default angles are used.

    Note
    ----
    The axis must have an equal aspect ratio for the lines to appear correctly oriented.

    Example
    -------
    >>> fig, ax = plt.subplots()
    >>> ax.plot([0, 1], [0, 1])
    >>> ax.set_aspect('equal')
    >>> annotateWLines(ax, angles=[-70, -30, 0, 30, 70])
    >>> plt.show()

    The lines are drawn from a point just above the current top of the axis, extending across the plotting area.
    """

    if angles is None:
        angles = [-70, -30, 0, 30, 70]
    
    bottom, top = ax.get_ylim()
    yrange = top - bottom
    newy = top + yrange * 0.02

    length = 1  # Adjust length as necessary
    
    for theta in angles:
        m = -1 / (np.tan(np.deg2rad(theta)) + 0.000000001)

        x0 = 0  # Starting x-coordinate, adjust as necessary
        y0 = newy
        dx = length**2 / (1 + m**2)
        dy = m * dx
        
        x1 = x0 + np.sqrt(dx)
        y1 = y0 + np.sqrt(dy) * np.sign(m)
        x2 = x0 - np.sqrt(dx)
        y2 = y0 - np.sqrt(dy) * np.sign(m)

        x = [x1, x2]
        y = [y1, y2]
        
        ax.plot(x, y, color='k')


def AngleHistograms(dikeset,lines, ax=None, fig=None, Trusted=True, Annotate=False):

    if ax is None:
        fig,ax=plt.subplots()

    ax.hist(dikeset['theta'], bins=np.arange(-90,100,10), density=True, facecolor='white', edgecolor='k', label='Segments')
    ax.hist(lines['AvgTheta'], bins=np.arange(-90,100,10), density=True, color='lightskyblue', alpha=0.5, label='All Clusters')

    if Trusted:
        ax.hist(lines[lines['TrustFilter']==1]['AvgTheta'], bins=np.arange(-90,100,10), density=True, color='mediumslateblue', alpha=0.5, label='Trusted Clusters')
    if Annotate:
        annotateWLines(ax)
    return ax

def AngleHistograms(dikeset, lines, ax=None, fig=None, Trusted=True, Annotate=False):
    """
    Plot histograms of angles for dike segments and clusters.

    This function creates histograms to compare the distribution of angles for dike segments and clusters, including an
    optional subset for 'trusted' clusters. It can also annotate the histogram with lines at specific angles if desired.

    Parameters
    ----------
    dikeset : pandas.DataFrame
        DataFrame containing the dike segments data, with a column 'theta' for their angles.
    lines : pandas.DataFrame
        DataFrame containing the clusters data, with columns 'AvgTheta' for their average angles and optionally
        'TrustFilter' to indicate trusted clusters.
    ax : matplotlib.axes._axes.Axes, optional
        Axes object on which to plot the histograms. If None, a new figure and axes are created. Default is None.
    fig : matplotlib.figure.Figure, optional
        Figure object for the plot. Only used if `ax` is None, in which case a new figure is created. Default is None.
    Trusted : bool, optional
        If True, includes a histogram for clusters marked as trusted based on the 'TrustFilter' column. Default is True.
    Annotate : bool, optional
        If True, adds annotation lines at specific angles to the histogram for visual reference. Default is False.

    Returns
    -------
    matplotlib.axes._axes.Axes
        The axes object used for the plot.

    Example
    -------
    >>> AngleHistograms(dikeset=df_segments, lines=df_clusters, Trusted=True, Annotate=True)
        Plots histograms of angles for dike segments and clusters, including 'trusted' clusters, with annotations.

    Note
    ----
    - The `dikeset` DataFrame must contain a 'theta' column representing the angles of dike segments.
    - The `lines` DataFrame must contain an 'AvgTheta' column for the average angles of clusters and can optionally
      contain a 'TrustFilter' boolean column to filter for trusted clusters.
    """
    import matplotlib.pyplot as plt
    import numpy as np

    if ax is None:
        fig, ax = plt.subplots()

    # Histogram for dike segments
    ax.hist(dikeset['theta'], bins=np.arange(-90, 100, 10), density=True, facecolor='white', edgecolor='k', label='Segments')
    
    # Histogram for all clusters
    ax.hist(lines['AvgTheta'], bins=np.arange(-90, 100, 10), density=True, color='lightskyblue', alpha=0.5, label='All Clusters')

    # Optional histogram for trusted clusters
    if Trusted:
        ax.hist(lines[lines['TrustFilter'] == 1]['AvgTheta'], bins=np.arange(-90, 100, 10), density=True, color='mediumslateblue', alpha=0.5, label='Trusted Clusters')
    
    # Optional annotation lines
    if Annotate:
        annotateWLines(ax)

    ax.legend()
    return ax


def DotsHT(fig, ax, lines, color=None, ColorBy="Dike Cluster Length (km)", label=None, cmap=cm.turbo, marker='o',
           rhoScale=True, Cbar=True, title=None, CbarLabels=True,
           axlabels=(True, True), StrOn=True, palette=None, alpha=0.4):
    """
    Create a scatter plot of rho and theta on a polar plot with customizable attributes.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure object to which the plot is added.
    ax : matplotlib.axes._subplots.PolarAxes
        The polar axes object on which to create the scatter plot.
    lines : pandas.DataFrame
        A DataFrame containing the data points to be plotted, with columns for rho and theta.
    color : str or None, optional
        The color of the data points. Default is None, which uses the default color.
    ColorBy : str, optional
        The name of the DataFrame column used to determine the color of the data points based on its values.
        Default is "Dike Cluster Length (km)".
    label : str or None, optional
        The label for the colorbar. Default is None, which uses the column name specified by `ColorBy`.
    cmap : matplotlib.colors.Colormap, optional
        The colormap used for coloring data points based on `ColorBy`. Default is `cm.turbo`.
    marker : str, optional
        The marker style for the data points. Default is 'o'.
    rhoScale : bool, optional
        If True, scales the rho values by dividing by 1000 to display in kilometers. Default is True.
    Cbar : bool, optional
        If True, displays a colorbar when coloring data points based on `ColorBy`. Default is True.
    title : str or None, optional
        The title of the scatter plot. Default is None.
    CbarLabels : bool, optional
        If True, displays tick labels on the colorbar. Default is True.
    axlabels : tuple of bool, optional
        Specifies whether to display labels for theta and rho axes. Default is (True, True).
    StrOn : bool, optional
        If True and `ColorBy` values are strings, enables string-based coloring. Default is True.
    palette : str or None, optional
        The name of the color palette to use for string-based coloring. Default is None.
    alpha : float, optional
        The transparency level of the data points. Default is 0.4.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object used for the plot.
    ax : matplotlib.axes._subplots.PolarAxes
        The polar axes object used for the plot.

    """

    t,r=whichForm(lines)

    if axlabels[0]:
        ax.set_xlabel(r'Theta ($^\circ$)')

    if not rhoScale:
        if axlabels[1]:
            ax.set_ylabel('Rho (m)')
            print("m scale label")
        rho=lines[r].values
    elif rhoScale:
        rho=lines[r].values/1000
        if axlabels[1]:
            ax.set_ylabel('Rho (km)')
            print("km scale label")


    if ColorBy==None and color==None:
        c='grey'
    if color is not None:
        c=color

    if ColorBy is not None:
        if type(lines[ColorBy].values[0]) is str:
            c,cmap=StringColors(lines[ColorBy].values, palette=palette)
            print("colorby value is type string")
        else:
            c=lines[ColorBy].values



    #ax[1], h2=HThist(lines['Average Rho (m)'], lines['Average Theta ($^\circ$)'], rstep, tstep, weights=lines['Dike Cluster Length (km)'], ax=ax[1],rbins=rbins)
    c2=ax.scatter(lines[t].values, rho, c=c, cmap=cmap, edgecolor='black', marker=marker, alpha=alpha)
    if title is not None:
        ax.set_title(title)


    if ColorBy is not None and Cbar:

        if type(lines[ColorBy].values[0]) is str:

            if StrOn and len(np.unique(c)) < 15:
                cbar=StringCbar(c2, fig, ax, lines[ColorBy].values.astype(str))
                cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation='vertical')
            else:
                cbar=fig.colorbar(c2, ax=ax)
        else:
            cbar=fig.colorbar(c2, ax=ax)

        if label is None:
            cbar.set_label(ColorBy)
        elif ColorBy is not None:
            cbar.set_label(label)

    if not CbarLabels and Cbar:
        cbar.ax.set_yticklabels([])
    ax.set_xlim([-90,90])
    plt.tight_layout()

    return fig,ax


def DotsLines(lines, ColorBy="seg_length", cmap=cm.turbo, linewidth=1, fig=None, ax=None, Cbar=True, CbarLabels=True, StrOn=False, color=None):
    """
    Plot line segments and their corresponding Hough Transform scatter plot side by side.

    Parameters
    ----------
    lines : pandas.DataFrame
        DataFrame containing data for line segments to be plotted.
    ColorBy : str, optional
        Column name in `lines` DataFrame to color the scatter plot points based on its values. Default is "seg_length".
    cmap : matplotlib.colors.Colormap, optional
        Colormap for the scatter plot points. Default is `cm.turbo`.
    linewidth : int, optional
        Line width for the line segments plot. Default is 1.
    fig : matplotlib.figure.Figure or None, optional
        The figure object for the plots. If None, a new figure is created. Default is None.
    ax : list of matplotlib.axes._subplots.PolarAxes or None, optional
        List of two polar axes objects for the plots. If None, new axes will be created. Default is None.
    Cbar : bool, optional
        If True, displays a colorbar for the scatter plot. Default is True.
    CbarLabels : bool, optional
        If True, displays labels on the colorbar. Default is True.
    StrOn : bool, optional
        If True and `ColorBy` is a string, uses string-based coloring for the scatter plot. Default is False.
    color : str or None, optional
        Color for the line segments. If None, uses default color. Default is None.

    Returns
    -------
    matplotlib.figure.Figure
        The figure object for the plots.
    list of matplotlib.axes._subplots.PolarAxes
        List of two polar axes objects for the plots.

    Example
    -------
    >>> fig, ax = DotsLines(lines_df)
    """


    if fig is None:
        fig,ax=plt.subplots(1,2)

    if color is None:
        c='k'
    else:
        c=color
    plotlines(lines, c, ax[0], ColorBy=ColorBy, cmap=cmap, linewidth=linewidth)

    DotsHT(fig, ax[1], lines, ColorBy=ColorBy, cmap=cmap,CbarLabels=CbarLabels, StrOn=StrOn, Cbar=Cbar, color=color)

    plt.tight_layout()
    FixAxisAspect(ax[1], ax[0])
    return fig,ax


def breakXaxis(xlim, numAxes=1):
    """
    function to break x axis into based in xlim
    based on matplotlib example
    https://matplotlib.org/stable/gallery/subplots_axes_and_figures/broken_axis.html

    num axes cannot be greater than 13

    input:
        xlim: tuple of x limits
        nAxes: number of axes you wish to make with the same breakpoints

    output:
        fig: figure object
        ax: list of axes objects
    """
    # f, axes = plt.subplots(numAxes,2)
    # ax=axes[:,0]
    # ax2=axes[:,1]
    try :
        numAxes>13
    except ValueError:
        print('You can not have a numAxes greater than 13 ')


    mosaic="""AAAB"""
    ax1labels=["A"]
    ax2labels=["B"]
    from string import ascii_uppercase
    j=2
    if numAxes>1:
        for i in range((numAxes-1)):

            letter1=ascii_uppercase[j]
            ax1labels.append(letter1)
            j=j+1
            letter2=ascii_uppercase[j]
            ax2labels.append(letter2)
            newline="\n"+letter1*3+letter2
            mosaic=mosaic+newline
            j=j+1

    print(mosaic)
    f = plt.figure(constrained_layout=True)
    ax_dict = f.subplot_mosaic(mosaic)
    #identify_axes(ax_dict)


    ax=[ax_dict[i] for i in ax1labels]
    ax2=[ax_dict[i] for i in ax2labels]


    d = .5  # proportion of vertical to horizontal extent of the slanted line
    kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
              linestyle="none", color='k', mec='k', mew=1, clip_on=False)

    for i in range(numAxes):
        if numAxes == 1:

            ax.set_xlim(xlim[0])
            ax2.set_xlim(xlim[1])

            ax.spines['right'].set_visible(False)
            ax2.spines['left'].set_visible(False)
            #ax.yaxis.tick_top()
            ax2.tick_params(labeleft=False, left=False)  # don't put tick labels at the top
            #ax2.yaxis.tick_bottom()

            #plots break symbols

            ax.plot([1, 1], [0, 1], transform=ax.transAxes, **kwargs)
            ax2.plot([0, 0], [0, 1], transform=ax2.transAxes, **kwargs)
            continue

        ax[i].set_xlim(xlim[0])
        ax2[i].set_xlim(xlim[1])

        ax[i].spines['right'].set_visible(False)
        ax2[i].spines['left'].set_visible(False)
        #ax[i].yaxis.tick_top()
        ax2[i].tick_params(labelleft=False, left=False)  # don't put tick labels at the top
        #ax2[i].yaxis.tick_bottom()

        #plots break symbols
        ax[i].plot([1, 1], [0, 1], transform=ax[i].transAxes, **kwargs)
        ax2[i].plot([0, 0], [0, 1], transform=ax2[i].transAxes, **kwargs)

    return f, ax, ax2

def splitData(xlim, x):
    """
    function to split data into two groups based on xlim
    assume only one breakpoint
    """
    x1=[]
    x2=[]
    for i in range(len(x)):
        if x[i]<max(xlim[0]):
            x1.append(x[i])
        else:
            x2.append(x[i])
    return x1, x2

def plotBreak(xlim, x, y, ax1, ax2,marker, **kwargs):
    """
    function to plot breakpoints
    """
    x1, x2=splitData(xlim, x)
    # ax1.plot(x1, y[0:len(x1)],marker,   **kwargs)
    # ax2.plot(x2, y[len(x1):], marker,  **kwargs)
    ax1.plot(x,y, marker, **kwargs)
    ax2.plot(x,y, marker, **kwargs)

def NumtoStringCoord(x,y):
    return "("+str(int(x))+","+str(int(y))+")"

def plotRadialOver(fig,ax1, ax2,xc,yc,Crange=50000,n=4, step=None, color='gray', color2="red", colorLines=False):
    ys=np.linspace(yc-Crange, yc+Crange,n)
    xs=np.linspace(xc-Crange, xc+Crange,n)
    xr,yr=np.meshgrid( xs, ys)
    cos=np.cos(np.linspace(-np.pi,np.pi))
    sin=np.sin(np.linspace(-np.pi,np.pi))
    r=[(x-xc)*cos+(y-yc)*sin for x,y in zip(xr.flatten(), yr.flatten())]
    i=0
    labelsAll=[str(x) for x in np.arange(0,len(xr.flatten()))] #[ NumtoStringCoord(x, y) for x,y in zip(xr.flatten(), yr.flatten())]
    colors = cm.rainbow(np.linspace(0, 1, len(ys)))
    # ax1.plot(xr, yr, '*', c=color2)
    # for i, label in zip(r, labels):
    #     ax2.plot( np.rad2deg(np.linspace(-np.pi,np.pi)), i, "-.",label=label, c=color, alpha=0.3)

    for y,c in zip(ys, colors):
        ax1.plot(xs, [y]*n, "*", c=c)

        r=[(x-xc)*cos+(y-yc)*sin for x in xs]
        labels=labelsAll[(i)*n :n*(i+1)]

        for rho, label in zip(r, labels):
            ax2.plot( np.rad2deg(np.linspace(-np.pi,np.pi)), rho, "-.",label=label, c=c, alpha=0.7)
        i=i+1

    #labellines.labelLines(ax2.get_lines())

    colors = cm.rainbow(np.linspace(0, 1, len(ys)))

    return r


def plotScatterHist(lines, x, y, hue=None, hue_norm=None, xlim=None, ylim=None, log_scale=(False, False), palette='Spectral', style=None, **kwargs):
    """
    Plot a scatter plot with marginal histograms for two variables and optionally color by a third variable.

    This function plots data from a DataFrame as a scatter plot along with histograms for the x and y variables
    on the top and right margins, respectively. The scatter plot points can be colored based on a third variable.
    The histograms help visualize the distribution of the x and y variables, providing a comprehensive view of the data.

    Parameters
    ----------
    lines : pandas.DataFrame
        The DataFrame containing the data to be plotted.
    x : str
        The column name in `lines` to be used for x-axis values in the scatter plot.
    y : str
        The column name in `lines` to be used for y-axis values in the scatter plot.
    hue : str, optional
        The column name in `lines` whose values are used to color the data points in the scatter plot. Default is None.
    hue_norm : tuple, optional
        The normalization range (min, max) for the hue variable. Applies only if `hue` is not None. Default is None.
    xlim : tuple, optional
        The limits for the x-axis as a tuple (min, max). Default is None, which autoscales the x-axis.
    ylim : tuple, optional
        The limits for the y-axis as a tuple (min, max). Default is None, which autoscales the y-axis.
    log_scale : tuple, optional
        A tuple of booleans specifying whether to apply a logarithmic scale to the x and y axes, respectively.
        Format is (x_log, y_log). Default is (False, False).
    palette : str or list, optional
        The color palette for the hue variable. Can be a string specifying a seaborn palette or a list of colors. Default is 'Spectral'.
    style : str, optional
        The marker style for the scatter plot. For example, 'o' for circles, 's' for squares. Default is None, which uses default markers.
    **kwargs : dict
        Additional keyword arguments passed to the scatter plot function.

    Returns
    -------
    matplotlib.figure.Figure
        The figure object containing the scatter plot and marginal histograms.
    list of matplotlib.axes._subplots.AxesSubplot
        A list containing the axes objects for the scatter plot, x-axis histogram, and y-axis histogram.

    Example
    -------
    >>> data = pd.DataFrame({'theta': np.random.rand(100), 'rho': np.random.rand(100), 'label': np.random.choice(['1', '2', '3'], 100)})
    >>> fig, axes = plotScatterHist(data, 'theta', 'rho', hue='label', palette='viridis')

    """
    sns.set_theme(style="ticks")

    fig = SetupAGUFig((115,190), 'landscape')
    gs = gridspec.GridSpec(4, 4)
    ax_main = plt.subplot(gs[1:4, :3])
    ax_xDist = plt.subplot(gs[0, :3])#,sharex=ax_main)
    ax_yDist = plt.subplot(gs[1:4, 3])#,sharey=ax_main)


    m=lines['Linked']==1
    if hue is not None:
        sns.scatterplot(lines, x=x, y=y, hue=hue,
                        palette=palette,
                        alpha=0.6, ax=ax_main,
                        edgecolor='k', hue_norm=hue_norm,
                        style=style)

        if len(np.unique(lines[hue])) < 10:
            h1=sns.histplot(lines, x=x, hue=hue,
                         multiple="stack",
                         palette=palette,
                         edgecolor=".3",
                         linewidth=.5,
                         ax=ax_xDist, stat='percent',
                         log_scale=(log_scale[0], False),
                         legend=False, hue_norm=hue_norm)

            h1.set(xticklabels=[])
            h1.set(xlabel=None)


            h2=sns.histplot(lines, y=y, hue=hue,
                         multiple="stack",
                         palette=palette,
                         edgecolor=".3",
                         linewidth=.5,
                         ax=ax_yDist, stat='percent',
                         log_scale=(False, log_scale[1]),
                         legend=False, hue_norm=hue_norm)
        else:
            h1=sns.histplot(lines, x=x,
                         edgecolor=".3",
                         linewidth=.5,
                         ax=ax_xDist, stat='percent',
                         log_scale=(log_scale[0], False),
                         legend=False)

            h1.set(xticklabels=[])
            h1.set(xlabel=None)

            h2=sns.histplot(lines, y=y,
                         edgecolor=".3",
                         linewidth=.5,
                         ax=ax_yDist, stat='percent',
                         log_scale=(False, log_scale[1]),
                         legend=False)


        lines=FilterLines(lines)
        h2.set(yticklabels=[])
        h2.set(ylabel=None)


    else:
        ax_main.scatter(lines[x].values, lines[y].values,  color='grey', alpha=0.4, edgecolor='k')

        ax_xDist.hist(lines[x].values, color='grey', alpha=0.4, weights=np.ones(len(lines[x].values)) / len(lines[x].values))

        ax_yDist.hist(lines[y].values, orientation='horizontal', color='grey', alpha=0.4, weights=np.ones(len(lines[y].values)) / len(lines[y].values))

        ax_yDist.set_xlabel('Percent')
        ax_yDist.xaxis.set_major_formatter(PercentFormatter(1))
        ax_xDist.set_ylabel('Percent')
        ax_xDist.yaxis.set_major_formatter(PercentFormatter(1))

    if log_scale[0]:
        #change x axis
        ax_main.set_xscale('log')
        ax_xDist.set_xscale('log')
    if log_scale[1]:
        #change x axis
        ax_main.set_yscale('log')
        ax_yDist.set_yscale('log')


    ax_main.set_xlabel(x)
    ax_main.set_ylabel(y)
    plt.tight_layout()


    return fig, [ax_main, ax_xDist, ax_yDist]

def plotRatioLine(ax, x, ratio, line_kw=None):
    """
    Plot a line with a specified ratio.

    This function plots a line on a given axis with a specified ratio (slope) by specifying the x values. You can customize the appearance of the line using the `line_kw` argument.

    Parameters:
    -----------
    ax : matplotlib.axes._subplots.AxesSubplot
        The axis object on which to plot the line.
    x : array-like
        The x values for the line.
    ratio : float
        The desired slope (ratio) of the line.
    line_kw : dict, optional 
        A dictionary of keyword arguments to customize the line's appearance (e.g., color, linestyle, label).

    Returns:
    --------
    ax : matplotlib.axes._subplots.AxesSubplot
        The modified axis object.
    l: list of matplotlib.lines.Line2D
        A list containing the line objects created.

    Example:
    --------
        # Plot a line with a 1:2 slope (y = 0.5 * x)
        fig, ax = plt.subplots()
        ax, line = plotRatioLine(ax, x=[0, 10], ratio=0.5, line_kw={'color': 'red', 'linestyle': '--', 'label': 'Line'})
        ax.legend()
    """


    xs=np.linspace(min(x), max(x))
    ys=ratio*xs
    l=ax.plot(xs,ys, **line_kw)
    return ax, l

def plotByLoc(lines, col, log_scale=(False, False)):
    """
    Plot histograms of a column against Xmid and Ymid, color-coded by the column values.

    This function creates a 2x2 subplot grid where the left side displays histograms of a specified column
    against 'Xmid' and 'Ymid', and the right side shows color bars corresponding to the histograms. The histograms
    are color-coded based on the values in the specified column, providing a visual distribution of the data. Options
    to apply a logarithmic scale to the x and y axes of the histograms are also available.

    Parameters
    ----------
    lines : pandas.DataFrame
        The DataFrame containing the data to be plotted. Must include 'Xmid', 'Ymid', and the specified column `col`.
    col : str
        The name of the column in `lines` DataFrame to plot and color-code against 'Xmid' and 'Ymid'.
    log_scale : tuple of bool, optional
        Specifies whether to apply a logarithmic scale on the x-axis and y-axis of the histograms, respectively.
        Format is (x_log, y_log). Default is (False, False).

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object containing the plots.
    axs : list of matplotlib.axes._subplots.AxesSubplot
        A list containing the two main axes objects for the histograms.
    axins : list of matplotlib.axes._axes.Axes
        A list containing the two inset axes objects for the color bars.

    """

    if "Xmid" not in lines.columns:
        lines=lines.assign( Xmid=(lines['Xstart']+lines['Xend'])/2, Ymid=(lines['Ystart']+lines['Yend'])/2)
    fig, ax=plt.subplots(1,2)


    axins1 = inset_axes(ax[0],
                    width="5%",  # width = 50% of parent_bbox width
                    height="25%",  # height : 5%
                    loc='upper right')

    axins2 = inset_axes(ax[1],
                    width="5%",  # width = 50% of parent_bbox width
                    height="25%",  # height : 5%
                    loc='upper right')


    cmap = sns.cubehelix_palette(start=0, light=1, as_cmap=True)
    cmap2 = sns.cubehelix_palette(start=2, light=1, as_cmap=True)
    g1=sns.histplot(data=lines, x='Xmid', y=col, ax=ax[0], cmap=cmap, cbar=True, stat='percent', cbar_ax=axins1, log_scale=log_scale)
    g2=sns.histplot(data=lines, x='Ymid', y=col, ax=ax[1], cmap=cmap2, cbar=True,  stat='percent', cbar_ax=axins2, log_scale=log_scale)

    axins1.yaxis.set_ticks_position("left")
    axins2.yaxis.set_ticks_position("left")

    return fig, [ax, axins1, axins2]
