import matplotlib.pyplot as plt
import numpy as np
import cartopy
from matplotlib import colors
from matplotlib import gridspec

import pandas as pd
import h3
from h3.unstable import vect  # used for hexgrid vectorization

# Claudio: renamed the class to hexGrid to avoid confusion with other functions
# I renamed plot_hex_hist to pcolorhex and made it a method of the hexGrid class
# I also renamed previosly names pcolorhex to plot_hexagons.
class hexGrid:
    """
    'Reference object' for creating transition matrices. This acts as a sort of template
    to form transition matrices with, as well as a dictionary to convert between integer
    indices and hex positions.

    To do [priority]:
     - Allow for shrinking the matrix based on empty elements [low]
     - Differentiate between 2D and 3D [low]
    """

    def __init__(self, hexagons, level_bounds=[0, np.inf], h3_res=3):
        """
        Initialize the grid.

        Parameters
        ----------
        hexagons : set
            The hexagons of the region of interest
        level_bounds : list
            List of floats defining the vertical level boundaries of the grid.
        h3_res : int
            Resolution of the H3 grid.
        """

        self.level_bounds = level_bounds
        self.h3_res = h3_res
        self.hexagons = list(hexagons)
        self.hexint = np.array([int(a, 16) for a in self.hexagons])
        self.n_levels = len(level_bounds) - 1  # Number of vertical levels
        self.levels = np.arange(self.n_levels)  # Levels of the grid
        self.n_hex = len(self.hexagons)  # Number of hexagons
        self.n = self.n_levels * self.n_hex  # Total number of bins
        self.edges = dict()  # Edges of the grid
        self.matidx = np.arange(0, self.n, dtype=np.uint(64))  # Matrix index

        # These series can be used to quickly convert between hexagon and matrix index
        self.map_to_mat_nodepth = pd.Series(
            index=self.hexagons, data=self.matidx[: self.n_hex]
        )
        self.map_from_mat_depth = pd.Series(
            index=self.matidx, data=np.repeat(self.levels, self.n_hex)
        )
        self.map_from_mat_hex = pd.Series(
            index=self.matidx, data=np.tile(self.hexagons, self.n_levels)
        )

    def compute_centroids(self):
        self.centroids = [h3.h3_to_geo(hex) for hex in self.hexagons]
        self.cen_lats = [c[0] for c in self.centroids]
        self.cen_lons = [c[1] for c in self.centroids]

    def info(self):
        print(
            f"Number of hexagons in the region with grid resolution {self.h3_res}:",
            len(self.hexagons),
        )
        print(f"Levels: {self.n_levels}")
        print(f"Total number of bins:", self.n)

    def count_2d(self, lon, lat, normalize=False):
        """
        Count particles in a 2D grid.

        Parameters
        ----------
        lon : np.array
            Array of longitudes
        lat : np.array
            Array of latitudes
        normalize : bool, optional
            If True, normalize the counts. Default is False.

        Returns
        -------
        np.array
            Array of counts per bin, normalized if specified.
        """

        interpolated_hexes_as_int = vect.geo_to_h3(
            lat, lon, self.h3_res
        )  # Interpolated hexes as integers
        count = np.zeros(self.n_hex, dtype=np.uint(64))  # Initialize count array

        counted_unique_int, counted_unique_int_occurences = np.unique(
            interpolated_hexes_as_int, return_counts=True
        )  # Count unique hexes
        counted_unique_hexes = [
            hex(hexagon)[2:] for hexagon in counted_unique_int
        ]  # Convert to hex strings
        self.miscount = 0  # Initialize miscount

        for hexagon_idx, hexagon in enumerate(counted_unique_hexes):  # Count particles
            try:
                count_idx = self.hexagons.index(hexagon)  # Get index of hexagon
                count[count_idx] = counted_unique_int_occurences[
                    hexagon_idx
                ]  # Count particles
            except ValueError:  # If the hexagon is not in the region of interest
                self.miscount += 1  # Increment miscount

        if normalize:
            total_count = np.sum(count)
            if total_count > 0:  # Avoid division by zero
                count = count / total_count  # Normalize counts

        # if self.miscount > 0:  # Print miscount
            # print(
            #     f"{self.miscount} particles were not counted because they were outside the region of interest."
            # )
        return count
    
    def plot_my_centroids(self, domain):
        centroids = [h3.h3_to_geo(hex) for hex in self.hexagons]
        cen_lats = [c[0] for c in centroids]
        cen_lons = [c[1] for c in centroids]

        fig = plt.figure()
        ax = plt.axes(projection=cartopy.crs.PlateCarree())
        ax.scatter(cen_lons, cen_lats, s=0.3, c="r", transform=cartopy.crs.PlateCarree())
        ax.coastlines()
        ax.gridlines(draw_labels=True, zorder=0, linestyle="--", linewidth=0.5)
        plt.show()

    def pcolorhex(self, counts, cmap="viridis", maxnorm=None, ax=None, draw_edges=False, alpha=1.0):
        """
        Plot a histogram of particle counts in a hexagonal grid

        Parameters
        ----------
        counts : array-like
            Array of particle counts. Should be the same length as the grid
        cmap : str, optional
            Colormap, by default 'viridis'
        maxnorm : int, optional
            Maximum value of the colorbar, by default None
        ax : matplotlib.axes.Axes, optional
            Axes to plot to, by default None
        """
        if not maxnorm:
            maxnorm = counts.max() if counts.size > 0 else 1

        if ax is None:
            fig, ax = plt.subplots(
                1, 1, subplot_kw={"projection": cartopy.crs.PlateCarree()}
            )

        cmap_instance = plt.cm.get_cmap(cmap)
        cmap_instance.set_bad("w")
        color_values = get_colors(counts, cmap_instance, 0, maxnorm)
        plot_hexagons(
            ax,
            self.hexagons,
            color_values,
            draw_edges=draw_edges,
            alpha=alpha,
            label="concentration",
        )

        sm = plt.cm.ScalarMappable(
            cmap=cmap_instance, norm=plt.Normalize(vmin=0, vmax=maxnorm)
        )

        return sm


# Aditional Functions
def plot_hexagons(
    ax,
    hexagons,
    colors=None,
    draw_edges=True,
    fill_polygons=True,
    transform=cartopy.crs.PlateCarree(),
    **kwargs,
):
    """
    Draw a collection of hexagons colored by a value. Based on a script from Mikael.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to draw to
    hexes : list
        A list of `h3` hexagons (i.e. related to tools.countGrid.hexagons)
    colors : list, optional
        A list of colors to use for each hexagon, by default None
    draw_edges : bool, optional
        Whether to draw the edges of the hexagons, by default True
    fill_polygons : bool, optional
        Whether to fill the hexagons with color, by default True
    **kwargs
        Additional keyword arguments passed to `matplotlib.pyplot.fill`

    Returns
    -------
    None
        This function does not return anything, but it does draw to the axes
        that are passed in.
    """
    for i1, hex_ in enumerate(hexagons):
        # get latitude and longitude coordinates of hexagon
        lat_long_coords = np.array(h3.h3_to_geo_boundary(str(hex_)))
        x = lat_long_coords[:, 1]
        y = lat_long_coords[:, 0]

        # ensure that all longitude values are between 0 and 360
        x_hexagon = np.append(x, x[0])
        y_hexagon = np.append(y, y[0])
        if x_hexagon.max() - x_hexagon.min() > 25:
            x_hexagon[x_hexagon < 0] += 360

        # draw edges
        if draw_edges:
            ax.plot(x_hexagon, y_hexagon, "k-", transform=transform, linewidth=0.2)

        # fill polygons
        if fill_polygons:
            ax.fill(
                x_hexagon, y_hexagon, color=colors[i1], transform=transform, **kwargs
            )


def get_colornorm(vmin=None, vmax=None, center=None, linthresh=None, base=None):
    """ "
    Return a normalizer

    Parameters
    ----------
    vmin : float (default=None)
        Minimum value of the data range
    vmax : float (default=None)
        Maximum value of the data range
    center : float (default=None)
        Center value for a two-slope normalization
    linthresh : float (default=None)
        Threshold for a symmetrical log normalization
    base : float (default=None)
        Base for a symmetrical log normalization

    Returns
    -------
    norm : matplotlib.colors.Normalize object
        A normalizer object
    """
    if type(base) is not type(None) and type(linthresh) is type(None):
        norm = colors.LogNorm(vmin=None, vmax=None)
    elif type(linthresh) is not type(None) and type(base) is not type(None):
        norm = colors.SymLogNorm(linthresh=linthresh, base=base, vmin=None, vmax=None)
    elif center:
        norm = colors.TwoSlopeNorm(vmin=vmin, vcenter=center, vmax=vmax)
    else:
        norm = colors.Normalize(vmin, vmax)
    return norm


def get_colors(
    inp, colormap, vmin=None, vmax=None, center=None, linthresh=None, base=0
):
    """ "
    Based on input data, minimum and maximum values, and a colormap, return color values

    Parameters
    ----------
    inp : array-like
        Input data
    colormap : matplotlib.colors.Colormap object
        Colormap to use
    vmin : float (default=None)
        Minimum value of the data range
    vmax : float (default=None)
        Maximum value of the data range
    center : float (default=None)
        Center value for a two-slope normalization
    linthresh : float (default=None)
        Threshold for a symmetrical log normalization
    base : float (default=None)
        Base for a symmetrical log normalization

    Returns
    -------
    colors : array-like
        Array of color values
    """
    norm = get_colornorm(vmin, vmax, center, linthresh, base)
    return colormap(norm(inp))


def pcolorhex(counts, grid, cmap="viridis", maxnorm=None, ax=None):
    """
    Creates a hexagonal binning plot of particle counts on a specified grid, with customizable colormap and normalization.

    Parameters
    ----------
    counts : array-like
        The particle counts for each hexagon in the grid. Length must match that of the grid.
    grid : tool.countGrid object
        An object representing the hexagonal grid layout.
    cmap : str, optional
        The colormap for the plot. Defaults to 'viridis'.
    maxnorm : int, optional
        The maximum value for color normalization. If None, uses the maximum count. Defaults to None.
    ax : matplotlib.axes.Axes, optional
        The matplotlib axes to plot on. If None, a new figure and axes are created. Defaults to None.

    Notes
    -----
    - The function automatically adjusts the color intensity of each hexagon based on the provided counts, using the specified colormap.
    - If `ax` is not provided, the function creates a new figure and axes with a PlateCarree projection.
    - The function returns a ScalarMappable object which can be used to create a colorbar.

    Returns
    -------
    matplotlib.cm.ScalarMappable
        A ScalarMappable instance created with the specified colormap and normalization, useful for creating a colorbar.
    """
    if not maxnorm:
        maxnorm = counts.max()

    if ax is None:
        fig, ax = plt.subplots(
            1, 1, subplot_kw={"projection": cartopy.crs.PlateCarree()}, cmap=plt.cm.viridis
        )

    plot_hexagons(
        ax,
        grid.hexagons,
        get_colors(counts, plt.cm.viridis, 0, maxnorm),
        draw_edges=False,
        alpha=1.0,
        label=" concentration",
    )

    cmap = plt.cm.get_cmap(cmap)
    cmap.set_bad("w")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=maxnorm))

    return sm


# function to convert hexint to hex
def int_to_hex(hexint):
    """
    Convert a list of integers to a list of hex strings.
    
    Parameters
    ----------
    hexint : list
        A list of integers.
        
    Returns
    -------
    list
        A list of hex strings.
    """
    return [hex(hexagon)[2:] for hexagon in hexint]

# function to convert hex to hexint
def hex_to_int(hex):
    """
    Convert a list of hex strings to a list of integers.
    
    Parameters
    ----------
    hex : list
        A list of hex strings.
        
    Returns
    -------
    list
        A list of integers.

    """
    return int(hex, 16)