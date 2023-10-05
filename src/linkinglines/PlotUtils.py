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
from HT import HT_center
from FitRectangle import *
from PrePostProcess import whichForm, FilterLines, getCartLimits
import scipy.cluster.hierarchy as sch
import labellines
import matplotlib.gridspec as gridspec
import string
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.colors import LogNorm

np.random.seed(5)

import matplotlib.ticker as mticker

class FixCartesianLabels():
    """
    Moves Offset axis tick label to axis label.

    This class is used to adjust the axis labels by moving the offset (scale factor) of the axis tick labels to the axis label itself. It's particularly helpful when dealing with plots where the offset notation is desired to be shown as part of the axis label.

    Attributes:
        None

    Methods:
        __init__(self, ax):
            Initializes an instance of the FixCartesianLabels class for a given axis.

        update(self, ax, lim):
            Updates the axis labels by moving the offset to the axis label.

    Example:
        import matplotlib.pyplot as plt

        # Create a sample plot
        fig, ax = plt.subplots()
        ax.plot([1, 2, 3], [1e6, 2e6, 3e6])

        # Initialize FixCartesianLabels for the y-axis
        y_label_fixer = FixCartesianLabels(ax.yaxis)

        plt.show()
    """
    def __init__(self, ax):
        # self.axis = {"y":ax.yaxis, "x":ax.xaxis}[axis]
        # self.labelx=""
        ax.callbacks.connect('Cartesian Plots Updated', self.update)
        ax.figure.canvas.draw()
        self.update(ax, None)

    def update(self, ax, lim):
        """
        Updates the axis labels by moving the offset to the axis label.

        Args:
            ax (matplotlib.axis.Axis): The axis for which the labels should be updated.
            lim: Unused parameter (needed for the callback).

        Returns:
            None
        """
        for i, l in zip([ax.yaxis, ax.xaxis], ['Y', 'X']):
            fmt = i.get_major_formatter()
            i.offsetText.set_visible(False)
            i.set_label_text(l + " (" + fmt.get_offset() + " m )")


def get_ax_size_inches(ax):
    """
    Get the size of a matplotlib axis in inches.

    This function calculates the size (width and height) of a matplotlib axis in inches, taking into account the current figure's DPI settings.

    Parameters:
        ax (matplotlib.axes._subplots.AxesSubplot): The axis for which the size should be calculated.

    Returns:
        tuple: A tuple containing the width and height of the axis in inches.

    Example:
        import matplotlib.pyplot as plt

        # Create a sample plot
        fig, ax = plt.subplots()

        # Get the size of the axis in inches
        width, height = get_ax_size_inches(ax)

        print(f"Width: {width} inches, Height: {height} inches")
    """
    fig = ax.get_figure()
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height

    return width, height
def FixAxisAspect(ax1, ax2):
    """
    Adjust the aspect ratios of two matplotlib axes to match each other.

    This function adjusts the aspect ratios of two axes to make them visually compatible. It ensures that the data displayed in both axes appears with the correct proportions.

    Parameters:
        ax1 (matplotlib.axes._subplots.AxesSubplot): The first axis to be adjusted.
        ax2 (matplotlib.axes._subplots.AxesSubplot): The second axis to be adjusted.

    Example:
        import matplotlib.pyplot as plt

        # Create two sample plots with different aspect ratios
        fig, ax1 = plt.subplots()
        ax2 = plt.axes([0.2, 0.2, 0.4, 0.4])

        # Adjust the aspect ratios to match
        FixAxisAspect(ax1, ax2)

        plt.show()
    """
    figW0, figH0 = ax1.get_figure().get_size_inches()
    # Axis size on figure
    _, _, w0, h0 = ax1.get_position().bounds
    # Ratio of display units
    disp_ratio0 = (figH0 * h0) / (figW0 * w0)

    w0i, h0i = get_ax_size_inches(ax1)

    figW1, figH1 = ax2.get_figure().get_size_inches()
    # Axis size on figure
    _, _, w1, h1 = ax2.get_position().bounds
    # Ratio of display units
    disp_ratio1 = (figH1 * h1) / (figW1 * w1)
    w1i, h1i = get_ax_size_inches(ax2)

    if h0i > h1i:
        y_min, y_max = ax2.get_ylim()
        yperinch = (y_max - y_min) / h1i
        deltah = h0i - h1i
        deltay = deltah * yperinch / 2
        ax2.set_ylim(y_min - deltay, y_max + deltay)
        print('reset Y')

    if w0i > w1i:
        x_min, x_max = ax2.get_xlim()
        xperinch = (x_max - x_min) / w1i
        deltaw = w0i - w1i
        deltax = deltaw * xperinch / 2
        ax2.set_xlim(x_min - deltax, x_max + deltax)
        print('reset x')


def labelSubplots(ax, labels=None, **kwargs):

    """
    Adds alphabet label to corner of each subplot

    Parameters:
        ax: list, dict, or array
        returns error for just one AxesSubplot object
        label: Default=None
        Defaults to "A", "B", "C " etc but
        could be any string

    Returns: none
    adds label to right top corner of plot

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
        How big the label should be.
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
        ax (matplotlib.axes._subplots.AxesSubplot): The axes for which to calculate the aspect ratio.

    Returns:
        float: The aspect ratio of the axes.

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
        vals     An RGB/RGBA tuple
        rgbtype  Valid valus are:
                          1 - Inputs are in the range 0 to 1
                        256 - Inputs are in the range 0 to 255

     Returns:
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
        c (list): A list of RGB or RGBA tuples.

    Returns:
        list: A list of hex strings in the form '#RRGGBB' or '#RRGGBBAA'.

    Example:
        # Convert an array of RGB (0-1 range) to Hex
        rgb_colors = [(0.2, 0.5, 0.8), (0.9, 0.1, 0.4)]
        hex_colors = RGBArraytoHexArray(rgb_colors)
        print("Hex Colors:", hex_colors)
    """
    return [RGBtoHex(i) for i in c]


def StringColors(values, palette="turbo"):
    """
    Map a list of strings to colors using a specified color palette.

    This function takes a list of strings and maps each unique string to a unique color from a specified color palette.

    Parameters:
        values (list): A list of strings to be mapped to colors.
        palette (str, optional): The name of the color palette to use. Defaults to "turbo".

    Returns:
        tuple: A tuple containing two elements:
            - color_idx (numpy.ndarray): An array of indices representing the colors for each string.
            - cm (matplotlib.colors.LinearSegmentedColormap): The colormap used for mapping the strings to colors.

    Example:
        # Map a list of categories to colors
        categories = ["Category A", "Category B", "Category C", "Category A"]
        color_indices, colormap = StringColors(categories, palette="viridis")
        print("Color Indices:", color_indices)
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
    """
    Create a colorbar for string-based categorical data.

    This function creates a colorbar for visualizing categorical data that has been mapped to colors using the StringColors function.

    Parameters:
        c (matplotlib.collections.Collection): The collection of colored elements (e.g., scatter points) in the plot.
        fig (matplotlib.figure.Figure): The figure object containing the plot.
        ax (matplotlib.axes.Axes): The axes object on which the plot is drawn.
        values (list): A list of strings representing the categorical data.

    Returns:
        matplotlib.colorbar.Colorbar: The colorbar object associated with the categorical data.

    Example:
        # Create a scatter plot with categorical colors
        values = ["Category A", "Category B", "Category C", "Category A"]
        color_indices, colormap = StringColors(values, palette="viridis")
        scatter = ax.scatter(x, y, c=color_indices, cmap=colormap)
        colorbar = StringCbar(scatter, fig, ax, values)
        plt.show()
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
    Generate list of items with fonts to change

    Parameters
    ----------
    fig : matplotlib figure object
        DESCRIPTION.
    ax : matplotlib axes object or list
        DESCRIPTION.

    Returns
    -------
    fontItems : list
        list of font items in fig and ax

    """

    if type(ax) is list:
        items=[]
        for i in ax:
            items.append([i.title, i.xaxis.label, i.yaxis.label] + i.get_xticklabels() + i.get_yticklabels())
        fontItems = [item for sublist in items for item in sublist]
    else:
        fontItems=[ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()

        fontItems.append(fig.get_legend().get_texts())
    return fontItems

def jgrSize(fig, ax, finalSize, units="mm"):
    '''
    Resize matplotlib figure to JGR appropriate size and set dpi to 600

    Parameters
    ----------
    fig : TYPE
        DESCRIPTION.
    finalSize : string or float
        "full", "half", "quarter" or float

    Returns
    -------
    fig : TYPE
        DESCRIPTION.

    '''

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
def SetupJGRFig(finalSize, orientation, units='mm'):
    """
    Set up a Matplotlib figure for creating a JGR (Journal of Geophysical Research) style plot.

    This function configures Matplotlib settings for creating a JGR-style plot, adjusting font sizes and figure size based on the desired final size and orientation.

    Parameters:
        finalSize (str, float, list, or tuple): The final size of the plot. Valid options are:
            - 'quarter': Quarter page size.
            - 'half': Half page size (landscape or portrait orientation).
            - 'full': Full page size.
            - A float value specifying the size as a fraction of full page size.
            - A list or tuple containing width and height values (in specified units).
        orientation (str): The orientation of the plot. Valid options are 'landscape' or 'portrait'.
        units (str): The units for specifying the final size if using a list or tuple. Valid options are 'mm', 'cm', or 'inches'.

    Returns:
        matplotlib.figure.Figure: The configured Matplotlib figure object for creating the plot.

    Example:
        # Set up a JGR-style figure with half page size in landscape orientation
        fig = SetupJGRFig('half', 'landscape')
        plt.plot(x, y)
        plt.xlabel('X Label')
        plt.ylabel('Y Label')
        plt.title('JGR-style Plot')
        plt.savefig('jgr_plot.png', dpi=300, bbox_inches='tight')
        plt.show()
    """
    SMALL_SIZE = 8
    MEDIUM_SIZE = 8
    BIGGER_SIZE = 8

    # Configure font sizes
    plt.rc('font', size=SMALL_SIZE)
    plt.rc('axes', titlesize=SMALL_SIZE)
    plt.rc('axes', labelsize=MEDIUM_SIZE)
    plt.rc('xtick', labelsize=SMALL_SIZE)
    plt.rc('ytick', labelsize=SMALL_SIZE)
    plt.rc('legend', fontsize=SMALL_SIZE)
    plt.rc('figure', titlesize=MEDIUM_SIZE)
    plt.rcParams['legend.title_fontsize'] = SMALL_SIZE

    # Create a Matplotlib figure
    fig = plt.figure()

    # Set the figure size based on finalSize and orientation
    if finalSize == 'quarter':
        w = 95 / 25.4  # cm/(cm/in)
        l = 115 / 25.4
    elif finalSize == 'half':
        if orientation == 'landscape':
            w = 190 / 25.4 / 2
            l = 230 / 25.4
        elif orientation == 'portrait':
            w = 190 / 25.4
            l = 230 / 25.4 / 2
    elif finalSize == 'full':
        w = 190 / 25.4
        l = 230 / 25.4
    elif isinstance(finalSize, float):
        w = 190 / 25.4 * finalSize
        l = 230 / 25.4 * finalSize
    elif isinstance(finalSize, (list, tuple)):
        if units == 'mm':
            w = finalSize[0] / 25.4
            l = finalSize[1] / 25.4
        elif units == 'cm':
            w = finalSize[0] / 2.54
            l = finalSize[1] / 2.54
        elif units == 'inches':
            w = finalSize[0]
            l = finalSize[1]
    else:
        w = 190 / 2 / 25.4
        l = 230 / 2 / 25.4

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

    Given a set of x and y coordinates, a central point (xmid, ymid), an angle (theta), and a length, this function
    calculates the coordinates of two points that represent a line segment with one end clustered around (xmid, ymid).

    Parameters:
        xs (array-like): Array of x-coordinates of data points.
        ys (array-like): Array of y-coordinates of data points.
        theta (float): Angle in degrees for the line segment.
        length (float): Length of the line segment.
        xmid (float, optional): X-coordinate of the central point. If None, it is calculated as the average of xs.
        ymid (float, optional): Y-coordinate of the central point. If None, it is calculated as the average of ys.

    Returns:
        tuple: A tuple containing four integers (x1, x2, y1, y2) representing the coordinates of two points that
        define the line segment.

    Example:
        # Calculate coordinates for a clustered line segment
        xs = [1, 2, 3, 4, 5]
        ys = [2, 3, 4, 5, 6]
        theta = 45  # Angle in degrees
        length = 3
        xmid, ymid = 3, 4  # Central point
        x1, x2, y1, y2 = clustered_lines(xs, ys, theta, length, xmid, ymid)
        print(f'Point 1: ({x1}, {y1})')
        print(f'Point 2: ({x2}, {y2})')
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
    Plot the rectangle defined by the center and lines, illustrating the orientation and dimensions.

    This function takes a set of lines, a center (xc, yc), and optionally a figure and axis to create a plot
    illustrating a rectangle that represents the orientation and dimensions of the lines with respect to the center.

    Parameters:
        lines (DataFrame): DataFrame containing line data, including angles (theta) and lengths (rho).
        xc (float): X-coordinate of the center.
        yc (float): Y-coordinate of the center.
        fig (matplotlib.figure.Figure, optional): The figure to use for plotting. If None, a new figure is created.
        ax (matplotlib.axes._subplots.AxesSubplot, optional): The axis to use for plotting. If None, a new axis is created.

    Returns:
        matplotlib.figure.Figure: The figure used for plotting.
        matplotlib.axes._subplots.AxesSubplot: The axis used for plotting.
        float: The length of the rectangle.
        float: The width of the rectangle.

    Example:
        # Plot a rectangle representing the orientation and dimensions of lines
        fig, ax, length, width = pltRec(lines_df, 0, 0)
        plt.show()
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
        t = 'Average Theta ($^\circ$)'

    if t == 'AvgTheta' or t == 'Average Theta ($^\circ$)':
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
        labels (list or pandas.Series): A list of labels.
        colormap (matplotlib.colors.Colormap): A colormap to assign colors from.

    Returns:
        list: A list of colors in hexadecimal format (#RRGGBB) corresponding to the input labels.
        list: A list of short color names (e.g., 'red', 'blue') corresponding to the input labels.

    Example:
        # Assign colors to unique labels using a colormap
        labels = ['A', 'B', 'A', 'C', 'B']
        colormap = plt.get_cmap('viridis')
        label_colors, short_colors = labelcolors(labels, colormap)
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
    Plots line segments on a specified axis.

    This function plots line segments based on input data, allowing customization of various plot attributes such as
    color, transparency, and more.

    Parameters:
        data (pandas.DataFrame): A DataFrame containing columns 'Xstart', 'Xend', 'Yend', and 'Ystart' representing
            UTM coordinates of line segments.
        col (str or RGB tuple): The color in which to plot the lines. Can be a string specifying a named color or
            an RGB tuple (e.g., (0.5, 0.5, 0.5)).
        ax (matplotlib.axes._axes.Axes): The axes object on which to plot the line segments.
        alpha (float, optional): The transparency level of the lines (0.0 for fully transparent, 1.0 for fully opaque).
            Default is 1.
        myProj (pyproj.Proj, optional): A PyProj projection object to convert UTM coordinates to lat/long. If None, UTM
            coordinates are assumed to be in lat/long format. Default is None.
        maskar (array-like, optional): A logical mask of the same length as the data, indicating which lines to plot.
            Default is None, which plots all lines.
        linewidth (float, optional): The width of the plotted lines. Default is 1.
        ColorBy (str, optional): A column name from the data DataFrame to color the lines based on a data attribute.
            Default is None, which results in a single color for all lines.
        center (bool, optional): If True, plots the center point of the data points on the map. Default is False.
        xc (float, optional): X-coordinate of the center point. Required if `center` is True and not provided in data.
        yc (float, optional): Y-coordinate of the center point. Required if `center` is True and not provided in data.
        extend (bool, optional): If True, extends the plotted lines beyond their endpoints. Default is False.
        cmap (matplotlib.colors.Colormap, optional): A colormap to use for coloring lines based on the `ColorBy`
            parameter. Default is the 'turbo' colormap.
        cbarStatus (bool, optional): If True, displays a colorbar when coloring lines based on `ColorBy`. Default is False.
        SpeedUp (bool, optional): If True and the data size is large, downsamples the data to speed up plotting.
            Default is True.
        equal (bool, optional): If True, ensures that the plot aspect ratio is equal, maintaining the correct scale.
            Default is True.

    Returns:
        None

    Example:
        # Plot line segments in red with transparency, color lines based on 'Value' column, and display a colorbar
        plotlines(data, 'red', ax, alpha=0.6, ColorBy='Theta', cbarStatus=True)
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


def HThist(lines, rstep, tstep, weights=None,fig=None, ax=None, rbins=None, tbins=None, cmap=cm.Blues, gamma=0.3):
    t,r=whichForm(lines)



    if fig is None or ax is None:
        fig,ax=plt.subplots()
    if rbins is None:
        rbins=np.arange(min(lines[r])/1000, max(lines[r])/1000, rstep)
    if tbins is None:
        tbins=np.arange(-90, 90, tstep)
    h,xe,ye, c=ax.hist2d(lines[t], lines[r]/1000, bins=[tbins, rbins], weights=weights, cmap=cmap, norm=mcolors.PowerNorm(gamma))
    fig.colorbar(c, label='Counts', ax=ax)
    ax.set_xlabel('Theta ($^\circ$)')
    ax.set_ylabel('Rho (km)')
    #ax.set_title("HT histogram")


    return fig, ax, [h, xe,ye,c]

def annotateWLines(ax, angles=None):
    #doesn't work unless axis are equal
    if angles is None:
        angles=[-70, -30, 1, 30, 70]
    bottom, top= ax.get_ylim()
    yrange=top-bottom
    newy=top+yrange*.02

    length=1 #0.5*yrange
    slopes=-1/(np.tan(np.deg2rad(angles))+0.000000001)
    for theta in angles:
        m=-1/(np.tan(np.deg2rad(theta))+0.000000001)

        x0 = theta
        y0 = newy
        dx=length**2/(2*(1+m**2))
        dy=dx*m
        x1 = (x0 + dx)
        y1 = (y0 + dy)
        x2 = (x0 - dx)
        y2 = (y0 - dy)

        x=[ x1, x2]
        y=[ y1, y2 ]
        print(x,y)
        ax.plot(x,y,color='k')


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


def BA_HT(dikeset,lines,rstep=5000):
    fig,ax=plt.subplots(1,3)
     #lines['StdRho'].mean()*2
    tstep=2
    rbins=np.arange(min(dikeset['rho']), max(dikeset['rho']), rstep)

    #ax[0],h1=HThist(dikeset['rho'], dikeset['theta'],rstep, tstep, weights=dikeset['seg_length'], ax=ax[0], rbins=rbins)
    ax[0],h1=HThist(dikeset['rho'], dikeset['theta'],rstep, tstep, ax=ax[0], rbins=rbins)
    ax[0].set_title('Raw Data')
    ax[0].set_xlabel('Theta (degrees)')
    ax[1].set_xlabel('Theta (degrees)')
    ax[2].set_xlabel('Theta (degrees)')
    ax[0].set_ylabel('Rho (m)')
    #ax[1], h2=HThist(lines['Average Rho (m)'], lines['Average Theta ($^\circ$)'], rstep, tstep, weights=lines['Dike Cluster Length (km)'], ax=ax[1],rbins=rbins)
    ax[1], h2=HThist(lines['Average Rho (m)'], lines['Average Theta ($^\circ$)'], rstep, tstep, ax=ax[1],rbins=rbins)
    ax[1].set_title('Clustered Data')
    fig.colorbar(h1[3], ax=ax[0])
    fig.colorbar(h2[3], ax=ax[1])
    hdiff= h2[0] - h1[0]
    x,y=np.meshgrid(h1[1], h1[2])
    divnorm = mcolors.TwoSlopeNorm(vcenter=0)
    c2=ax[2].pcolormesh(x,y,hdiff.T, cmap=cm.RdBu, norm=divnorm)
    ax[2].set_title('Change')
    fig.colorbar(c2, ax=ax[2])

    plt.tight_layout()

    return fig,ax, h1, h2



def DotsHT(fig, ax, lines, color=None, ColorBy="Dike Cluster Length (km)", label=None, cmap=cm.turbo, marker='o',
           rhoScale=True, Cbar=True, title=None, CbarLabels=True,
           axlabels=(True, True), StrOn=True, palette=None, alpha=0.4):
    """
    Create a scatter plot of rho and theta

    This function creates a scatter plot of data points on a polar plot, allowing for customization of various plot
    attributes such as colors, markers, scales, and more.

    Parameters:
        fig (matplotlib.figure.Figure): The Figure object to place the plot on.
        ax (matplotlib.axes._subplots.PolarAxes): The polar axes object on which to create the scatter plot.
        lines (pandas.DataFrame): A DataFrame containing data points to be plotted.
        color (str or None, optional): The color of the data points. Can be a string specifying a named color or
            None to use default color. Default is None.
        ColorBy (str, optional): The name of the DataFrame column to color the data points based on its values.
            Default is "Dike Cluster Length (km)".
        label (str or None, optional): The label to be shown in the colorbar if ColorBy is used. Default is None.
        cmap (matplotlib.colors.Colormap, optional): A colormap to use for coloring data points based on ColorBy.
            Default is the 'turbo' colormap.
        marker (str, optional): The marker style for data points. Default is 'o' (circle).
        rhoScale (bool, optional): If True, scales the rho values by dividing by 1000 to display in kilometers.
            Default is True.
        Cbar (bool, optional): If True, displays a colorbar when coloring data points based on ColorBy. Default is True.
        title (str or None, optional): The title of the scatter plot. Default is None (no title).
        CbarLabels (bool, optional): If True, displays tick labels on the colorbar; otherwise, hides them. Default is True.
        axlabels (tuple, optional): A tuple (ax_theta_label, ax_rho_label) specifying whether to display axis labels for
            theta and rho. Default is (True, True).
        StrOn (bool, optional): If True and ColorBy values are strings, enables string-based colorbar and labels.
            Default is True.
        palette (str or None, optional): The name of the color palette to use for string-based coloring when ColorBy
            values are strings. Default is None.
        alpha (float, optional): The transparency level of the data points (0.0 for fully transparent, 1.0 for fully
            opaque). Default is 0.4.

    Returns:
        matplotlib.figure.Figure: The modified Figure object.
        matplotlib.axes._subplots.PolarAxes: The modified polar axes object.

    Example:
        # Create a scatter plot of data points colored by the 'Value' column
        fig, ax = DotsHT(fig, ax, data, ColorBy='theta', cmap=cm.inferno)
    """

    t,r=whichForm(lines)

    if axlabels[0]:
        ax.set_xlabel('Theta ($^\circ$)')

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
    Create a side-by-side plot with line segments on the left and a scatter plot on the right.

    This function creates a side-by-side plot with line segments on the left panel and a scatter plot on the right panel.
    The scatter plot can be customized by specifying the ColorBy column, colormap, linewidth, and more.

    Parameters:
        lines (pandas.DataFrame): A DataFrame containing line segment data.
        ColorBy (str, optional): The name of the DataFrame column to color the scatter plot points based on its values.
            Default is "seg_length".
        cmap (matplotlib.colors.Colormap, optional): A colormap to use for coloring the scatter plot points based on ColorBy.
            Default is the 'turbo' colormap.
        linewidth (int, optional): The linewidth of the line segments in the left panel. Default is 1.
        fig (matplotlib.figure.Figure or None, optional): The Figure object to place the plot on. If None, a new Figure will be created.
        ax (list of matplotlib.axes._subplots.PolarAxes or None, optional): A list of two polar axes objects (left and right panels).
            If None, new axes will be created. Default is None.
        Cbar (bool, optional): If True, displays a colorbar when coloring scatter plot points based on ColorBy. Default is True.
        CbarLabels (bool, optional): If True, displays tick labels on the colorbar; otherwise, hides them. Default is True.
        StrOn (bool, optional): If True and ColorBy values are strings, enables string-based colorbar and labels.
            Default is False.
        color (str or None, optional): The color of the line segments. Can be a string specifying a named color or None to use default color.
            Default is None.

    Returns:
        matplotlib.figure.Figure: The modified Figure object.
        list of matplotlib.axes._subplots.PolarAxes: A list of two polar axes objects (left and right panels).

    Example:
        # Create a side-by-side plot of line segments and a scatter plot colored by the 'Value' column
        fig, ax = DotsLines(data, ColorBy='Length', cmap=cm.inferno)
    """

    if fig is None:
        fig,ax=plt.subplots(1,2)
        # fig = SetupJGRFig('quarter', 'landscape')
        # ax1 = fig.add_axes([0.5, 0.1, 0.9, 0.9])
        # ax2 = fig.add_axes([0.1, 0.1, 0.4, 0.9])
        # ax=[ax1,ax2]
    if color is None:
        c='k'
    else:
        c=color
    plotlines(lines, c, ax[0], ColorBy=ColorBy, cmap=cmap, linewidth=linewidth)


    #ax[1], h2=HThist(lines['Average Rho (m)'], lines['Average Theta ($^\circ$)'], rstep, tstep, weights=lines['Dike Cluster Length (km)'], ax=ax[1],rbins=rbins)
    DotsHT(fig, ax[1], lines, ColorBy=ColorBy, cmap=cmap,CbarLabels=CbarLabels, StrOn=StrOn, Cbar=Cbar, color=color)
    #ax[1].set_title('HT')



    plt.tight_layout()
    FixAxisAspect(ax[1], ax[0])
    return fig,ax

def DotsLinesHist(lines, rstep, tstep, cmap1=cm.turbo, cmap2=cm.gray, ColorBy=None):
    t,r=whichForm(lines)
    if ColorBy is None:
        ColorBy=t
    fig,ax=plt.subplots(1,3)    #lines['StdRho'].mean()*2
    plotlines(lines, 'k', ax[0], ColorBy=ColorBy, cmap=cmap1, center=True, alpha=0.4)
    ax[0].set_title('Cartesian')
    #ax[1], h2=HThist(lines['Average Rho (m)'], lines['Average Theta ($^\circ$)'], rstep, tstep, weights=lines['Dike Cluster Length (km)'], ax=ax[1],rbins=rbins)
    DotsHT(fig, ax[1], lines, ColorBy=ColorBy, cmap=cmap1)
    ax[1].set_title('HT')


    fig,ax, h=HThist(lines, rstep, tstep, cmap=cmap2, fig=fig, ax=ax[2])

    plt.tight_layout()

    return fig, ax


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
    Create a scatter plot with histograms for two variables.

    This function creates a scatter plot with histograms for two variables (x and y) and the option to color points by a third variable (hue).

    Parameters:
        lines (pandas.DataFrame): A DataFrame containing the data to plot.
        x (str): The name of the column to use for the x-axis.
        y (str): The name of the column to use for the y-axis.
        hue (str, optional): The name of the column to use for coloring points.
        hue_norm (tuple, optional): A tuple specifying the normalization range for the hue variable (min, max).
        xlim (tuple, optional): A tuple specifying the x-axis limits (min, max).
        ylim (tuple, optional): A tuple specifying the y-axis limits (min, max).
        log_scale (tuple, optional): A tuple specifying whether to use a log scale for the x and y axes (x_log, y_log).
        palette (str or list, optional): The color palette to use for hue values.
        style (str, optional): The style of the scatter plot (e.g., 'o', 's', 'D').
        **kwargs: Additional keyword arguments to pass to the scatterplot function.

    Returns:
        matplotlib.figure.Figure: The modified Figure object.
        list of matplotlib.axes._subplots.AxesSubplot: A list of three axes objects (scatter plot, x-axis histogram, and y-axis histogram).

    Example:
        # Create a scatter plot with histograms
        fig, axes = plotScatterHist(data_df, x='X', y='Y', hue='Z', xlim=(0, 100), log_scale=(False, True))
    """
    # Function code goes here...

    sns.set_theme(style="ticks")

    fig = SetupJGRFig((115,190), 'landscape')
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
        #xbins=np.arange(0,1,0.1)
        #ax_xDist.hist(lines['Overlap'].values, bins=xbins)

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


            #ybins=np.arange(0,90,5)
            #ax_yDist.hist(lines['EnEchelonAngleDiff'].values, orientation='horizontal', bins=ybins)
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


            #ybins=np.arange(0,90,5)
            #ax_yDist.hist(lines['EnEchelonAngleDiff'].values, orientation='horizontal', bins=ybins)
            h2=sns.histplot(lines, y=y,
                         edgecolor=".3",
                         linewidth=.5,
                         ax=ax_yDist, stat='percent',
                         log_scale=(False, log_scale[1]),
                         legend=False)


        lines=FilterLines(lines)
        h2.set(yticklabels=[])
        h2.set(ylabel=None)

        #ax_yDist.set(yticklabels=[], ylabel=None)
        #ax_xDist.set(xticklabels=[], xlabel=None)
        #ax_main.scatter(lines[x].values, lines[y].values,  color='white', alpha=0.4, edgecolor='k', marker='*')

    else:
        ax_main.scatter(lines[x].values, lines[y].values,  color='grey', alpha=0.4, edgecolor='k')

        ax_xDist.hist(lines[x].values, color='grey', alpha=0.4, weights=np.ones(len(lines[x].values)) / len(lines[x].values))

        ax_yDist.hist(lines[y].values, orientation='horizontal', color='grey', alpha=0.4, weights=np.ones(len(lines[y].values)) / len(lines[y].values))

        ax_yDist.set_xlabel('Percent')
        ax_yDist.xaxis.set_major_formatter(PercentFormatter(1))
        #ax_yDist.set(yticklabels=[])
        #ax_xDist.set(xticklabels=[])
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
        ax (matplotlib.axes._subplots.AxesSubplot): The axis object on which to plot the line.
        x (array-like): The x values for the line.
        ratio (float): The desired slope (ratio) of the line.
        line_kw (dict, optional): A dictionary of keyword arguments to customize the line's appearance (e.g., color, linestyle, label).

    Returns:
        matplotlib.axes._subplots.AxesSubplot: The modified axis object.
        list of matplotlib.lines.Line2D: A list containing the line objects created.

    Example:
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
    Plot histograms of a column against Xmid and Ymid, with color-coded distributions.

    This function creates a 2x2 subplot grid with two main plots on the left and two color bars (inset axes) on the right. It uses seaborn's histplot to create histograms of a specified column col, mapping the color to the values in col. The function also provides options to use a logarithmic scale on the x-axis and y-axis of the main plots.

    Parameters:
        lines (pandas.DataFrame): The DataFrame containing the data.
        col (str): The column name to plot against Xmid and Ymid.
        log_scale (tuple of bool, optional): A tuple of two boolean values specifying whether to use a logarithmic scale for the x-axis and y-axis of the main plots, respectively. Default is (False, False).

    Returns:
        matplotlib.figure.Figure: The created figure object.
        list of matplotlib.axes._subplots.AxesSubplot: A list containing the main axis objects.
        list of matplotlib.axes._axes.Axes: A list containing the inset (color bar) axis objects.

    Example:
        # Plot histograms of a column 'Value' against Xmid and Ymid with logarithmic y-axis scale
        fig, [ax, axins1, axins2] = plotByLoc(lines, col='Value', log_scale=(False, True))
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
