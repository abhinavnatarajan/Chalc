import chalc as ch
import numpy as np

# Standard 2-simplex
K = ch.standard_simplex(2)
# Flatten out K
K.flat_representation()
# Change the filtration value of the vertex labelled [0, 1]
K.simplices[1][0].filtration_value = 2.0
# Check that the flat representation is sorted by filtration values
K.flat_representation()
# Change the vertex colours
K.simplices[0][0].colours = 1
K.simplices[0][1].colours = 1
K.simplices[0][2].colours = 2
K.propagate_colours()
K.flat_representation()

# Delaunay triangulation of N points on the circle
N = 999
x = np.array([(np.cos(theta), np.sin(theta)) for theta in np.array(range(N))/(2*np.pi)]).transpose()
K = ch.delaunay_complex(x)

# Split the points into two classes randomly
colours = np.random.randint(2, size=N)
# Stratify the points
y = ch.stratify(x, colours)
# y is now 3-dimensional
y.shape

# Weak alpha chromatic complex of x with previously chosen random colours
L = ch.weak_chromatic_alpha_complex(x, colours)
# Get the number of simplices in L
print(L.num_simplices)
# Check if L has the simplex [0, 4, 1, 5]
L.has_simplex([0, 4, 1, 5])
# Get the list of 2-simplices in L
L.simplices[2]
# Get handle to an arbitrary 2-simplex
temp = np.random.randint(len(L.simplices[2]))
s = L.simplices[2][list(L.simplices[2].keys())[temp]]
# Get the vertex labels of s
s.vertices
# Get the facets of s
s.facets
# Check that the [0]th face of s has the same vertices as the last 2 vertices of s
s.facets[0].vertices == s.vertices[1:]
# Check that the facets appear before s
all([f.filtration_value <= s.filtration_value for f in s.facets])
# Get the 2-skeleton of L
K = ch.FilteredComplex(L, 2)
# Check that s is present in K
K.has_simplex(s.vertices)
# Get the filtration value of s
s_filt_value = s.filtration_value
# s is a handle, so if we change its filtration value then its value changes in both K and L
s.filtration_value = -1 # only for demonstration!
K.simplices[2][s.label].filtration_value == L.simplices[2][s.label].filtration_value
L.simplices[2][s.label].filtration_value == -1.0
# Rectify the changed filtration value by propagating the filtration values upwards
L.propagate_filt_values(start_dim=1, upwards=True)
# Check that the filtration value of s is back to its original value
# this is exact, not just an approximate floating point equality
s.filtration_value == s_filt_value 