import matplotlib.pyplot as plt
import numpy as np
import meshio

# Assumptions:
# 1) The center of the connecting rod head is at point (0,0,0).
# 2) Axes are assigned as follows: X (s1), Y (s2), Z (s3).
# 3) Plane (x,y,0) is a symmetry plane of the connecting rod.

# Load the mesh from STL file
mesh = meshio.read("conrod_suzuki_si.msh")
#mesh = meshio.read("conrod_rabaman_si.msh")

# Extract points and cells (tetrahedrons)
points = np.array(mesh.points)
cells = mesh.cells_dict["tetra"]

# Connecting distance between centers of small and large hole [m]
L = 0.100
# Crank radius [m]
R = 0.0279
# Coefficient lambda = R/L
l = R / L
# Shaft rotational speed
f = 9500
# Shaft angular speed
omega = (2*np.pi/60) * f
# Density rho = 7800 [kg/m^3]
rho = 7722
# Starting volume [m^3]
V = 0
# Starting mass [kg]
mc = 0
# Starting values of moments m_{c,1} and m_{c,2} [m*kg]
mc1 = 0
mc2 = 0

# Numerical integration - summation by tetrahedrons
for i in range(len(cells)):
    # Take tetrahedron vertices indices
    ti = cells[i]
    # Convert vertices indices to (x,y,z) coordinates:
    pi = points[ti, :].T
    # Convert vertices A,B,C,D to 3 vectors with origin at A
    wi = pi[:, 1:4] - pi[:, 0:1]
    # Volume of tetrahedron ABCD as 1/6 of vectors volume
    Vi = abs(np.linalg.det(wi)) / 6
    V += Vi
    # Compute tetrahedron mass
    mci = rho * Vi
    mc += mci
    # Tetrahedron ABCD centroid
    d = np.sum(pi, axis=1) / 4
    # Integration by s1*dm and s2*dm
    mc1 += d[0] * mci
    mc2 += d[1] * mci

# Print volume in [m^3]
print("Volume: {} m^3".format(V))
# Print mass in [kg]
print("Mass: {} kg".format(mc))
# Print moment m_{c,1} in [m*kg]
print("Moment m_c1: {} m*kg".format(mc1))
# Print moment m_{c,2} in [m*kg]
print("Moment m_c2: {} m*kg".format(mc2))

# Prepare plots of inertia force

# Space of angles
a = np.linspace(0, 2 * np.pi, 360*2)
# X coordinate of inertia force
def inertia_x(R, L, omega, mc, mc1, mc2, a):
    l = R/L
    return R*omega*omega*(mc1/L*np.cos(a) + (mc-mc1/L)*(np.cos(a)+l*np.cos(2*a)) + mc2/L*np.sin(a))
# Y coordinate of inertia force
def inertia_y(R, L, omega, mc1, mc2, a):
    l = R/L
    return R*omega*omega*(-mc1/L*np.sin(a) + mc2/L*np.cos(a) - mc2/L*(np.cos(a)+l*np.cos(2*a)))
# Wector of Fx and Fy coordinates for all angles
Fx = inertia_x(R, L, omega, mc, mc1, mc2, a)
Fy = inertia_y(R, L, omega, mc1, mc2, a)
# Magniture of inertia force
absF = np.sqrt(Fx**2 + Fy**2)
# Red points on inertia force vector trajectory
red_angles = np.array([0, np.pi/2, np.pi, 3*np.pi/2])
red_angles_txt = ['$0$', '$\pi/2$', '$\pi$', '$3\pi/2$']
red_points_x = inertia_x(R, L, omega, mc, mc1, mc2, red_angles)
red_points_y = inertia_y(R, L, omega, mc1, mc2, red_angles)
# Increase fons size
plt.rcParams.update({
    'font.size': 14,
    'axes.titlesize': 20,
    'axes.labelsize': 18,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 16,
    'figure.titlesize': 20
})
# Set up figure
plt.figure(figsize=(18, 6))
# Inertia force rotation plot
plt.subplot(1, 3, 1)
plt.plot(Fx, Fy, label='$F(\\alpha)$', linewidth=3)
# Set origin of polar system (0, 0)
plt.plot(0, 0, 'ro')
plt.text(300, 300, '$O$', fontsize=16, verticalalignment='bottom')
# Set red points on trajectory for alpha = 0, pi/2, pi, 3pi/2
for i, angle in enumerate(red_angles):
    plt.plot(red_points_x[i], red_points_y[i], 'ro')
    plt.text(red_points_x[i] + 100, red_points_y[i] + 300, red_angles_txt[i], fontsize=16)
plt.title('Inertia force rotation')
plt.xlabel('$F_x$ [N]')
plt.ylabel('$F_y$ [N]')
plt.axis('equal')
plt.grid(True)
plt.legend(loc='lower right')

# Red points on inertia force magnituda plot
red_angles = np.array([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
red_angles_txt = ['$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$']
absFx = inertia_x(R, L, omega, mc, mc1, mc2, red_angles)
absFy = inertia_y(R, L, omega, mc1, mc2, red_angles)
red_points_x = red_angles
red_points_y = np.sqrt(absFx*absFx + absFy*absFy)
# Magnitude of inertia force plot
plt.subplot(1, 3, (2,3))
plt.plot(a, absF, label='$|F(\\alpha)|$', linewidth=3)
# Set red points on plot for alpha = 0, pi/2, pi, 3pi/2, 2pi
for i, angle in enumerate(red_angles):
    plt.plot(red_points_x[i], red_points_y[i], 'ro')
plt.title('Magnitude of inertia force')
plt.xlabel('$\\alpha$ [rad]')
plt.xticks(ticks=[0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi], labels=['$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'])
plt.ylabel('\n\n$|F|$ [N]')
plt.grid(True)
plt.legend(loc='lower right')

plt.tight_layout()
plt.savefig('num_example_1.eps', format='eps')
plt.show()
