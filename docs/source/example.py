# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Example Usage

# %% nbsphinx="hidden"
# %matplotlib inline
# %config InlineBackend.figure_formats = ['svg']
from IPython.core.display import display, Markdown
def printmd(object):
    display(Markdown(str(object).replace('\n', '<br>')))

# %%
import chalc as ch
import numpy as np
import matplotlib.pyplot as plt

# %% nbsphinx="hidden"
plt.rcParams["animation.html"] = "jshtml"
import h5py, logging, sys

# Disable some deprecation warnings from the seaborn plotting library
logging.basicConfig(stream=sys.stdout, level=logging.WARNING)
logging.captureWarnings(True)

# %% [markdown]
# For our data we sample 100 points on a circle with some noise and 100 points from inside the unit disk.

# %%
rng = np.random.default_rng(40)
num_points_circle = 100
num_points_disk = 100
mean = [0, 0]
cov = np.eye(2)*0.01
circle = np.array([[np.sin(2*np.pi*t), np.cos(2*np.pi*t)] for t in rng.random(num_points_circle)]).T +\
    rng.multivariate_normal(mean, cov, num_points_circle).T # points as columns
disk = np.sqrt(rng.random(num_points_disk)) * np.array([[np.sin(2*np.pi*t), np.cos(2*np.pi*t)] for t in rng.random(num_points_disk)]).T
points = np.concatenate((circle, disk), axis=1)
plt.scatter(circle[0, :], circle[1, :], s=10)
plt.scatter(disk[0, :], disk[1, :], s=10)
plt.gca().set_aspect('equal')
colours = [0]*num_points_circle + [1]*num_points_disk
plt.title('Point cloud')
plt.show()


# %% [raw] raw_mimetype="text/restructuredtext"
# We can compute the chromatic alpha/Delaunay--Cech/Delaunay--Rips complex $K$ of the point cloud using the functions :func:`chromatic.alpha <chalc.chromatic.alpha>`, :func:`chromatic.delcech <chalc.chromatic.delcech>`, and :func:`chromatic.delrips <chalc.chromatic.delrips>` respectively.
# In each of these cases, :math:`K` has far fewer simplices than either the Čech or Vietoris--Rips complex, which each have :math:`\displaystyle \binom{200}{2} = 19900` edges and :math:`\displaystyle \binom{200}{3} = 1313400` 2-simplices.

# %%
K, _ = ch.chromatic.alpha(points, colours)
printmd(f'$K$ has {len(K.simplices[1])} 1-simplices and {len(K.simplices[2])} 2-simplices')

# %% [raw] raw_mimetype="text/restructuredtext"
# :math:`K` is a :class:`FilteredComplex <chalc.filtration.FilteredComplex>` object, and simplices in the filtration are objects of type :class:`Simplex <chalc.filtration.Simplex>` (although you should not generally need to interact with the API of these objects).
# We can visualise the filtration at any given time using :func:`plotting.draw_filtration <chalc.plotting.draw_filtration>`. For example, at time 0.3 we have:

# %%
fig, ax = plt.subplots(1, 3, figsize=(12, 4))
ch.plotting.draw_filtration(K, points, time=0.3, include_colours=[0], ax=ax[0])
ax[0].set_title('Colour 0')
ch.plotting.draw_filtration(K, points, time=0.3, include_colours=[1], ax=ax[1])
ax[1].set_title('Colour 1')
ch.plotting.draw_filtration(K, points, time=0.3, ax=ax[2])
ax[2].set_title('Both colours')
for i in range(3):
    ax[i].set_xlim(-1.2, 1.2)
    ax[i].set_ylim(-1.2, 1.2)
fig.suptitle('Complexes at $t=0.3$')
plt.tight_layout()
plt.show()

# %% [raw] raw_mimetype="text/restructuredtext"
# The functions :func:`sixpack.from_filtration <chalc.sixpack.from_filtration>` and :func:`sixpack.compute <chalc.sixpack.compute>` compute the 6-pack of persistent homology diagrams of the chromatic alpha/Delaunay--Čech complex/Delaunay--Rips complexes.

# %%
# using the previously computed filtration
dgms_alpha = ch.sixpack.from_filtration(K, dom=[0])
# or directly from the point cloud
dgms_delcech = ch.sixpack.compute(points, colours, dom=[0], method="delcech")
dgms_delrips = ch.sixpack.compute(points, colours, dom=[0], method="delrips")

# %% [raw] raw_mimetype="text/restructuredtext"
# Each 6-pack is an object of type :class:`DiagramEnsemble <chalc.sixpack.DiagramEnsemble>`, with attributes ``ker``, ``cok``, ``im``, ``dom``, ``cod``, ``rel``, ``entrance_times``, and ``dimensions``.
# Each of the first 6 of those attributes is an object of type :class:`SimplexPairings <chalc.sixpack.SimplexPairings>` that has sets of paired and unpaired simplices, while the last two attributes are lists indicate the filtration time and dimension of those simplices. 
# For example-

# %%
print(dgms_alpha.ker)

# %% [raw] raw_mimetype="text/restructuredtext"
# A convenience function :func:`get <chalc.sixpack.DiagramEnsemble.get>` is provided to extract the diagrams as a matrix of birth/death times:

# %%
np.set_printoptions(threshold=10)
print(dgms_alpha.get('ker', 1)) # get the kernel in dimension one
# dgms_alpha.get('ker', [0, 1]) to get the kernel in dimensions zero and one
# dgms_alpha.get('ker') to get the kernel in all dimensions from zero to max(dgms_alpha.dimensions)

# %% [raw] raw_mimetype="text/restructuredtext"
# The convenience functions :func:`plotting.plot_sixpack <chalc.plotting.plot_sixpack>` and :func:`plotting.plot_diagram <chalc.plotting.plot_diagram>` are provided for plotting either a 6-pack or a specific diagram from a 6-pack.

# %%
fig1, ax1 = ch.plotting.plot_sixpack(dgms_alpha, tolerance = 1e-3)
fig1.suptitle('Chromatic Alpha')
fig1.set_figwidth(15)
fig1.set_figheight(9)
fig1.subplots_adjust(top=0.92)
plt.show()
fig2, ax2 = plt.subplots(1, 3)
# You can specify the dimensions of features to include in the plot.
ch.plotting.plot_diagram(dgms_alpha, 'rel', dimensions = {0, 2}, truncation = 0.3, ax = ax2[0], tolerance = 1e-3)
ax2[0].set_title('Chromatic Alpha')
# You can also specify a single dimension as an integer.
ch.plotting.plot_diagram(dgms_delcech, 'rel', dimensions = 1, truncation = 0.3, ax = ax2[1], tolerance = 1e-3)
ax2[1].set_title('Chromatic Delaunay-Cech')
# If no dimensions are specified, all features will be included.
ch.plotting.plot_diagram(dgms_delrips, 'rel', truncation = 0.3, ax = ax2[2], tolerance = 1e-3)
ax2[2].set_title('Chromatic Delaunay-Rips')
fig2.suptitle('Individual relative diagrams')
fig2.set_figwidth(15)
plt.show()

# %% [markdown]
# We can save the diagrams to file or load a diagram from file:

# %%
with h5py.File('test.h5', 'w') as f:
    ch.sixpack.save_diagrams(dgms_alpha, f)

with h5py.File('test.h5', 'r') as f:
    dgms_alpha = ch.sixpack.load_diagrams(f)

# %% [markdown]
# We can visualise the 2-skeleton of the filtration for points in 2D:

# %%
animation = ch.plotting.animate_filtration(
    K, points, filtration_times=list(np.linspace(0, 1.0, 45)), animation_length=5)
animation
