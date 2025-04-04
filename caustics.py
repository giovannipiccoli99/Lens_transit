import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

# Define optical depth function in cylindrical coordinates
def l_OD_test(h, phi):
    # Example function with height (h) dependence and some phi dependence
    # h is the height above the reference radius
    # phi is the azimuthal angle
    z = np.exp(-h/20) * (1 + 0.1 * np.cos(2*phi))
    return z

# def l_OD_test(h, phi):
#     # Base exponential profile with height
#     base = np.exp(-h/20)
    
#     # Add latitude-dependent features (φ=0 is equator in this example)
#     # For example: a band at latitude 30° (φ=π/6)
#     feature = 0.2 * np.exp(-((phi - np.pi/6)**2)/(0.1**2))
    
#     return base * (1 + feature)

# Create grid for interpolation
h_data = np.arange(0, 200, 1)
phi_data = np.arange(-np.pi, np.pi, 0.1)
l_OD_data = l_OD_test(h_data[:, None], phi_data[None, :])

# Create interpolation function
l_OD = interpolate.RectBivariateSpline(h_data, phi_data, l_OD_data)

# Plot the optical depth profile at phi=0
plt.figure(figsize=(10, 6))
plt.plot(h_data, l_OD.ev(h_data, np.zeros_like(h_data)), label='Optical Depth at φ=0')
plt.xlabel('Height (h)')
plt.ylabel('Optical Depth')
plt.title('Optical Depth Profile')
plt.grid(True)
plt.legend()
plt.savefig('optical_depth_profile.png')
plt.close()

# Define magnification matrix in cylindrical coordinates
def inv_M(h, phi, args):
    R, D = args  # R is planet radius, D is distance to source
    r = h + R  # Total radius from center

    M_hh = 1/D + l_OD.partial_derivative(2, 0)(h, phi)
    M_pp = r[:, None]/D + l_OD.partial_derivative(1, 0)(h, phi)/r[:, None] + l_OD.partial_derivative(0, 2)(h, phi)/r[:, None]**2
    M_hp = l_OD.partial_derivative(0, 1)(h, phi)/r[:, None]**2 - l_OD.partial_derivative(1, 1)(h, phi)/r[:, None]
    
    res = np.array([[M_hh, M_hp], [M_hp, M_pp]])
    res = np.transpose(res, (2, 3, 0, 1))
    
    return res

# Calculate magnification determinant
def magnification_det(h, phi, args):
    A = np.zeros([2, 2, np.size(h), np.size(phi)])
    det_inv_M = np.zeros([np.size(h), np.size(phi)])
    A = inv_M(h, phi, args)
    det_inv_M = np.linalg.det(A)
    return (det_inv_M)

# Parameters
R = 10          # Planet radius
D = 100000      # Distance to source
args = (R, D) # Tuple containing parameters
# Create a grid of h and phi values
h = np.arange(0, 200, 0.1)
phi = np.arange(-np.pi, np.pi, 0.1)
h_mesh, phi_mesh = np.meshgrid(h, phi)  # Create a meshgrid for h and phi
# Calculate determinant
det_inv_M = magnification_det(h, phi, args)

# Plot the determinant as a function of phi at a fixed height
plt.figure(figsize=(10, 6))
h_fixed_idx = 10  # Some index for a fixed height
plt.plot(phi, det_inv_M[h_fixed_idx, :])
plt.xlabel('φ (radians)')
plt.ylabel('Determinant')
plt.title(f'Determinant at h = {h[h_fixed_idx]:.1f}')
plt.grid(True)
plt.savefig('determinant_fixed_height.png')
plt.close()

# Plot the determinant as a 2D heatmap
plt.figure(figsize=(12, 8))
plt.imshow(det_inv_M.T, origin='lower', 
           extent=[h.min(), h.max(), phi.min(), phi.max()],
           aspect='auto', cmap='viridis')
plt.colorbar(label='Determinant')
plt.xlabel('Height (h)')
plt.ylabel('φ (radians)')
plt.title('Magnification Determinant in Cylindrical Coordinates')

# Find critical curves (where determinant = 0)
contour = plt.contour(h_mesh, phi_mesh, np.transpose(det_inv_M), levels=[0], colors='red')
plt.clabel(contour, inline=True, fontsize=8)
plt.savefig('determinant_contours.png')
plt.close()

# Extract critical curves
critical_curves = []
# Check if contour has paths and handle differently
if len(contour.allsegs) > 0:
    for segment in contour.allsegs[0]:  # allsegs[0] contains segments for the first level
        h_critical = segment[:, 0]
        phi_critical = segment[:, 1]
        critical_curves.append((h_critical, phi_critical))

# Plot critical curves separately if needed
if critical_curves:
    plt.figure(figsize=(10, 8))
    for h_cr, phi_cr in critical_curves:
        plt.plot(h_cr, phi_cr, 'r-')
    plt.xlabel('Height (h)')
    plt.ylabel('φ (radians)')
    plt.title('Critical Curves in (h, φ) Space')
    plt.grid(True)
    plt.savefig('critical_curves.png')
    plt.close()

# Function to map critical curves to caustics in the source plane
def find_caustics(h_cr, phi_cr, args):
    R, D = args
    r_cr = h_cr + R
    
    # Evaluate derivatives with grid=False to treat inputs as points
    dtau_dh = np.array([l_OD.partial_derivative(1, 0)(h, p, grid=False) for h, p in zip(h_cr, phi_cr)])
    dtau_dphi = np.array([l_OD.partial_derivative(0, 1)(h, p, grid=False) for h, p in zip(h_cr, phi_cr)])
    
    # Calculate caustic coordinates
    h_caustic = h_cr - D * dtau_dh
    phi_caustic = phi_cr - (D * dtau_dphi) / r_cr**2
    
    return h_caustic, phi_caustic

# Calculate and plot caustics if critical curves exist
if critical_curves:
    plt.figure(figsize=(10, 8))
    for h_cr, phi_cr in critical_curves:
        h_caustic, phi_caustic = find_caustics(h_cr, phi_cr, args)
        plt.plot(h_caustic, phi_caustic, 'b-')
    plt.xlabel('Height in Source Plane')
    plt.ylabel('φ in Source Plane')
    plt.title('Caustics Mapped to Source Plane')
    plt.grid(True)
    plt.savefig('caustics_source_plane.png')
    plt.close()