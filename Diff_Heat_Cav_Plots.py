import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.interpolate import make_interp_spline

# Physical data
Ra = 1e-6
L = 1.0
lambd = 1.0
rho = 1.0
cp = 0.71
alpha = lambd / (rho * cp)

# Choose plots
print("Choose the graphs you want to visualize:")
print("- To display colormaps of T, u, and v, type: 'COLOR'")
print("- To plot velocity profiles at L/2, type: 'VELOCITY'")
print("- To plot Grid Independence Study graph, type: 'GRID'")
graph = input("Plot: ").strip().upper()

if graph == 'COLOR':

    ################### COLORMAP: TEMPERATURE T ##################

    # Read temperature data file (formatted with 3 columns)
    df = pd.read_csv("TemperatureDistribution.txt", sep=r"\s+", engine="python", header=None, names=["x", "y", "T"])
    x = df["x"].values
    y = df["y"].values
    T = df["T"].values

    # Original grid resolution
    nx = len(np.unique(x))
    ny = len(np.unique(y))
    scale = 3  # Upscaling for finer interpolation

    # Create interpolated grid
    xi = np.linspace(np.min(x), np.max(x), nx * scale)
    yi = np.linspace(np.min(y), np.max(y), ny * scale)
    Xi, Yi = np.meshgrid(xi, yi)
    Ti = griddata((x, y), T, (Xi, Yi), method="cubic")

    # Plot contourf map
    fig, ax = plt.subplots(figsize=(6, 6))
    cf = ax.contourf(Xi, Yi, Ti, levels=100, cmap='jet', vmin=0, vmax=1)

    # Black isolines
    contours = ax.contour(Xi, Yi, Ti, levels=np.arange(0, 1.01, 0.1), colors='black', linewidths=0.7)
    ax.clabel(contours, inline=True, fontsize=7, fmt="%.1f")

    # Colorbar
    #cbar = plt.colorbar(cf, ax=ax)
    #cbar.set_label("Temperature")

    # Labels and title
    ax.set_xlabel("X (m)")
    ax.set_ylabel("Y (m)")

    #ax.set_title("Temperature Map with Ra = {:.0e}".format(Ra))
    #ax.set_aspect("equal")

    # Layout and save
    plt.tight_layout()
    plt.savefig("TemperatureMap_Ra_{:.0e}.png".format(Ra), dpi=100)
    plt.show()


    ################### COLORMAP: VELOCITY U ##################

    # Read velocity-u data file
    df_u = pd.read_csv("VelocityUDistribution.txt", sep=r"\s+", engine="python", header=None, names=["x", "y", "u"])
    x_u = df_u["x"].values
    y_u = df_u["y"].values
    u_vals = df_u["u"].values

    # Original grid resolution
    nx_u = len(np.unique(x_u))
    ny_u = len(np.unique(y_u))
    scale_u = 3  # Upscaling for finer interpolation

    # Create interpolated grid
    xi_u = np.linspace(np.min(x_u), np.max(x_u), nx_u * scale_u)
    yi_u = np.linspace(np.min(y_u), np.max(y_u), ny_u * scale_u)
    Xi_u, Yi_u = np.meshgrid(xi_u, yi_u)
    Ui = griddata((x_u, y_u), u_vals, (Xi_u, Yi_u), method="cubic")

    # Plot contourf map
    fig_u, ax_u = plt.subplots(figsize=(6, 6))
    cf_u = ax_u.contourf(Xi_u, Yi_u, Ui, levels=100, cmap='jet')

    # Black isolines
    vmin_u = np.percentile(Ui, 1)
    vmax_u = np.percentile(Ui, 99)
    levels_u = np.linspace(vmin_u, vmax_u, 11)  # 10 intervals

    contours_u = ax_u.contour(Xi_u, Yi_u, Ui, levels=levels_u, colors='black', linewidths=0.8)
    ax_u.clabel(contours_u, inline=True, fontsize=7, fmt="%.1f")

    # Colorbar
    #cbar_u = plt.colorbar(cf_u, ax=ax_u)
    #cbar_u.set_label("Velocity u (m/s)")

    # Labels and title
    ax_u.set_xlabel("X (m)")
    ax_u.set_ylabel("Y (m)")

    #ax_u.set_title("Velocity u(y) Map with Ra = {:.0e}".format(Ra))
    #ax_u.set_aspect("equal")
    ax_u.set_xlim(0.01, 0.99)
    ax_u.set_ylim(0, 1)
    ax_u.set_aspect('equal')

    # Layout and save
    plt.tight_layout()
    plt.savefig("VelocityUMap_Ra_{:.0e}.png".format(Ra), dpi=100)
    plt.show()



    ################### COLORMAP: VELOCITY V ##################

    # Read velocity-v data file
    df_v = pd.read_csv("VelocityVDistribution.txt", sep=r"\s+", engine="python", header=None, names=["x", "y", "v"])
    x_v = df_v["x"].values
    y_v = df_v["y"].values
    v_vals = df_v["v"].values

    # Original grid resolution
    nx_v = len(np.unique(x_v))
    ny_v = len(np.unique(y_v))
    scale_v = 3  # Upscaling for finer interpolation

    # Create interpolated grid
    xi_v = np.linspace(np.min(x_v), np.max(x_v), nx_v * scale_v)
    yi_v = np.linspace(np.min(y_v), np.max(y_v), ny_v * scale_v)
    Xi_v, Yi_v = np.meshgrid(xi_v, yi_v)
    Vi = griddata((x_v, y_v), v_vals, (Xi_v, Yi_v), method="cubic")

    # Plot contourf map
    fig_v, ax_v = plt.subplots(figsize=(6, 6))
    cf_v = ax_v.contourf(Xi_v, Yi_v, Vi, levels=100, cmap='jet')

    # Black isolines
    vmin_v = np.percentile(Vi, 1)
    vmax_v = np.percentile(Vi, 99)
    levels_v = np.linspace(vmin_v, vmax_v, 11)
    contours_v = ax_v.contour(Xi_v, Yi_v, Vi, levels=levels_v, colors='black', linewidths=0.8)
    ax_v.clabel(contours_v, inline=True, fontsize=7, fmt="%.1f")

    # Colorbar
    #cbar_v = plt.colorbar(cf_v, ax=ax_v)
    #cbar_v.set_label("Velocity v (m/s)")

    # Labels and title
    ax_v.set_xlabel("X (m)")
    ax_v.set_ylabel("Y (m)")
    #ax_v.set_title("Velocity v(x) Map with Ra = {:.0e}".format(Ra))
    #ax_v.set_aspect("equal")
    ax_v.set_xlim(0, 1)
    ax_v.set_ylim(0.01, 0.99)
    ax_v.set_aspect('equal')

    # Layout and save
    plt.tight_layout()
    plt.savefig("VelocityVMap_Ra_{:.0e}.png".format(Ra), dpi=100)
    plt.show()



elif graph == 'VELOCITY':

    ################### VELOCITY PROFILE u(y) at x = L/2 ##################

    # Read data from 3-column file
    df = pd.read_csv("VelocityUDistribution.txt", delim_whitespace=True, header=None, names=["x", "y", "u"])

    # Filter rows where x = L/2 (= 0.5)
    filtered = df[abs(df["x"] - 0.5) < 1e-6]

    # Extract data for plotting
    y_vals_u = filtered["y"].values
    u_vals = filtered["u"].values

    # Sort for smooth plotting
    sorted_indices = y_vals_u.argsort()
    y_vals_u = y_vals_u[sorted_indices]
    u_vals = u_vals[sorted_indices]

    # Spline interpolation for smoother curve
    y_smooth = np.linspace(y_vals_u.min(), y_vals_u.max(), 300)
    spline = make_interp_spline(y_vals_u, u_vals, k=3)
    u_smooth = spline(y_smooth)

    # Plot interpolated curve
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(u_smooth, y_smooth, color='red', linewidth=2)
    ax.set_xlabel("u (m/s)")
    ax.set_ylabel("y (m)")
    # ax.set_title(f"Velocity Profile u(y) at x = L/2 with Ra = {Ra:.0e}")
    ax.grid(True)
    ax.set_ylim([0, 1])
    plt.tight_layout()
    plt.savefig(f"Profile_Uy_Ra_{Ra:.0e}.png", dpi=100)
    plt.show()


    # Print maximum u* velocity and y position at x = L/2
    filename = 'VelocityUDistribution.txt'

    # Load data without header
    df = pd.read_csv(filename, delim_whitespace=True, header=None, names=['x', 'y', 'velocity'])

    # Filter rows where x = 0.5
    filtered = df[abs(df['x'] - 0.5) < 1e-6]

    # Find row with max velocity
    if not filtered.empty:
        max_row = filtered.loc[filtered['velocity'].idxmax()]
        print(f"Maximum u* velocity at x = L/2 is {max_row['velocity'] * L / alpha}, occurring at y = {max_row['y']}")
    else:
        print("No data found with x = 0.5")



    ################### VELOCITY PROFILE v(x) at y = L/2 ##################

    # Read data from 3-column file
    df = pd.read_csv("VelocityVDistribution.txt", delim_whitespace=True, header=None, names=["x", "y", "v"])

    # Filter rows where y = L/2 (= 0.5)
    filtered = df[abs(df["y"] - 0.5) < 1e-6]

    # Extract data for plotting
    x_vals_v = filtered["x"].values
    v_vals = filtered["v"].values

    # Sort for smooth plotting
    sorted_indices = x_vals_v.argsort()
    x_vals_v = x_vals_v[sorted_indices]
    v_vals = v_vals[sorted_indices]

    # Spline interpolation for smoother curve
    x_smooth = np.linspace(x_vals_v.min(), x_vals_v.max(), 300)
    spline = make_interp_spline(x_vals_v, v_vals, k=3)
    v_smooth = spline(x_smooth)

    # Plot interpolated curve
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(x_smooth, v_smooth, color='blue', linewidth=2)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("v (m/s)")
    # ax.set_title(f"Velocity Profile v(x) at y = L/2 with Ra = {Ra:.0e}")
    ax.grid(True)
    ax.set_xlim([0, 1])
    plt.tight_layout()
    plt.savefig(f"Profile_Vx_Ra_{Ra:.0e}.png", dpi=100)
    plt.show()


    # Print maximum v* velocity and x position at y = L/2
    filename = 'VelocityVDistribution.txt'

    # Load data without header
    df = pd.read_csv(filename, delim_whitespace=True, header=None, names=['x', 'y', 'velocity'])

    # Filter rows where y = 0.5
    filtered = df[abs(df['y'] - 0.5) < 1e-6]

    # Find row with max velocity
    if not filtered.empty:
        max_row = filtered.loc[filtered['velocity'].idxmax()]
        print(f"Maximum v* velocity at y = L/2 is {max_row['velocity'] * L / alpha}, occurring at x = {max_row['x']}")
    else:
        print("No data found with y = 0.5")

elif graph == 'GRID':

    # Grid sizes used (expressed as N for N×N)
    grid_points = [21, 31, 41, 51]

    # Relative errors [%] for each Rayleigh number and each quantity
    errors_data = {
        "Ra=1e3": {
            "Nu": [1.237, 0.575, 0.224, 0],
            "umax": [4.77, 1.89, 0.61, 0],
            "vmax": [5.10, 2.20, 0.73, 0],
        },
        "Ra=1e4": {
            "Nu": [2.67, 1.14, 0.38, 0],
            "umax": [3.46, 1.80, 0.68, 0],
            "vmax": [4.59, 1.47, 0.40, 0],
        },
        "Ra=1e5": {
            "Nu": [6.55, 2.56, 0.89, 0],
            "umax": [6.63, 2.92, 0.94, 0],
            "vmax": [6.51, 0.16, 0.40, 0],
        },
        "Ra=1e6": {
            "Nu": [17.19, 6.73, 2.30, 0],
            "umax": [10.03, 5.61, 1.62, 0],
            "vmax": [5.07, 2.11, 5.73, 0],
        }
    }

    # Create one subplot for each quantity
    fig, axs = plt.subplots(3, 1, figsize=(8, 12), sharex=True)

    quantities = ["Nu", "umax", "vmax"]
    titles = [
        "Relative Error on Average Nusselt Number",
        "Relative Error on Max Horizontal Velocity $u_{max}$",
        "Relative Error on Max Vertical Velocity $v_{max}$"
    ]

    # Plot each quantity with all Rayleigh numbers
    for i, quantity in enumerate(quantities):
        for ra, data in errors_data.items():
            axs[i].plot(grid_points, data[quantity], marker='o', label=ra)
        axs[i].set_ylabel("Relative Error [%]")
        axs[i].set_title(titles[i])
        axs[i].grid(True)
        axs[i].legend()

    axs[2].set_xlabel("Grid Size (N × N)")
    plt.tight_layout()
    plt.savefig("grid_independence_errors.png")
    plt.show()


else:
    print("Invalid input. Please choose 'COLOR' or 'VELOCITY'.")
