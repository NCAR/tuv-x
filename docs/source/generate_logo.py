import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import matplotlib.textpath
import matplotlib.patches
import numpy as np
from PIL import Image
import os


# https://matplotlib.org/stable/gallery/lines_bars_and_markers/gradient_bar.html
def gradient_image(ax, extent, direction=0.3, cmap_range=(0, 1), **kwargs):
    """
    Draw a gradient image based on a colormap.

    Parameters
    ----------
    ax : Axes
        The axes to draw on.
    extent
        The extent of the image as (xmin, xmax, ymin, ymax).
        By default, this is in Axes coordinates but may be
        changed using the *transform* keyword argument.
    direction : float
        The direction of the gradient. This is a number in
        range 0 (=vertical) to 1 (=horizontal).
    cmap_range : float, float
        The fraction (cmin, cmax) of the colormap that should be
        used for the gradient, where the complete colormap is (0, 1).
    **kwargs
        Other parameters are passed on to `.Axes.imshow()`.
        In particular useful is *cmap*.
    """
    phi = direction * np.pi / 2
    v = np.array([np.cos(phi), np.sin(phi)])
    X = np.array([[v @ [1, 0], v @ [1, 1]],
                  [v @ [0, 0], v @ [0, 1]]])
    a, b = cmap_range
    X = a + (b - a) / X.max() * X
    im = ax.imshow(X, extent=extent, interpolation='bicubic',
                   vmin=0, vmax=1, **kwargs)
    return im

def make_logo(tuvx_static_path = '_static', 
              dpi = 192, 
              width = 1200, 
              height = 288, 
              alpha = 1.0):
  """ 
  Create the TUV-x logo

    Parameters
    ----------
    tuvx_static_path : str, default 'static_'
        The directory that the logo and favicon are saved to
    dpi: int, default 192
      The dots per inch. This is used to calculate the figsize.
    width : int, default 1200
      The width of the figure, in pixels
    height : int, default 288
      The height of the figure, in pixels
    alpha: float, defuault 1.0
      The transparency applied to the background gradient and the edge of the 
      'TUV-x' text
  """

  # https://stackoverflow.com/questions/13714454/specifying-and-saving-a-figure-with-exact-size-in-pixels
  figsize = (width/dpi, height/dpi)
  fig, ax = plt.subplots(figsize=figsize)

  # define text before gradient to get extent
  fp = FontProperties(family='DejaVu Sans', weight='bold')
  text = matplotlib.textpath.TextPath((0.0, 0.0), 'TUV-x',
                                              size=1, prop=fp)

  extent = text.get_extents().extents[[0, 2, 1, 3]]
  im = gradient_image(ax, direction=1, extent=extent,
                cmap=plt.cm.winter, cmap_range=(0.2, 0.8), alpha=alpha)

  im.set_clip_path(text, transform=ax.transData)

  tuv_patch = matplotlib.patches.PathPatch(text, facecolor='none',
     edgecolor='violet', alpha=alpha, linewidth=2)
  ax.add_patch(tuv_patch)

  ax.spines[['top', 'bottom', 'left', 'right']].set_visible(False)
  ax.set_xticks(())
  ax.set_yticks(())

  # the default ylim cuts off the very bottom of the U for som reason
  # ax.set_ylim((-.2, ax.get_ylim()[1]))

  ax.tick_params(width=0, which='both')

  logo = r'logo.png'

  fig.savefig(logo, format='png', transparent=True, dpi=dpi)
  fig.savefig(f'{tuvx_static_path}/logo.svg', format='svg', transparent=True)

  img = Image.open(logo)
  img.save(f'{tuvx_static_path}/favicon.ico')

  os.remove(logo)
  
if __name__ == '__main__':
  make_logo()