import matplotlib.pyplot as plt
import numpy as np
import meshio

# Connecting distance between centers of small and large hole [m]
L = 0.2667
# Crank radius [m]
R = 0.073
# Coefficient lambda = R/L
l = R / L
# Shaft rotational speed
f = 1000
# Shaft angular speed
omega = (2*np.pi/60) * f
# Density rho = 7800 [kg/m^3]
rho = 7722
# Starting mass [kg]
mc = 8.221
# Starting values of moments m_{c,1} and m_{c,2} [m*kg]
mc1 = L*5.775
mc2 = L*1.6442

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
plt.savefig('num_example_2.eps', format='eps')
plt.show()
