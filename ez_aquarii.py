import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

def read_data(filename):
    """
    Reads the simulation data from file
    """
    data = {
        'time': [],
        'A_x': [],
        'A_y': [],
        'A_z': [],
        'A_vx': [],
        'A_vy': [],
        'A_vz': [],
        'B_x': [],
        'B_y': [],
        'B_z': [],
        'B_vx': [],
        'B_vy': [],
        'B_vz': [],
        'C_x': [],
        'C_y': [],
        'C_z': [],
        'C_vx': [],
        'C_vy': [],
        'C_vz': []
    }
    with open(filename, 'r') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split()
            if len(parts) < 19:
                continue  # Skip incomplete lines
            try:
                data['time'].append(float(parts[0]))
                data['A_x'].append(float(parts[1]))
                data['A_y'].append(float(parts[2]))
                data['A_z'].append(float(parts[3]))
                data['A_vx'].append(float(parts[4]))
                data['A_vy'].append(float(parts[5]))
                data['A_vz'].append(float(parts[6]))
                data['B_x'].append(float(parts[7]))
                data['B_y'].append(float(parts[8]))
                data['B_z'].append(float(parts[9]))
                data['B_vx'].append(float(parts[10]))
                data['B_vy'].append(float(parts[11]))
                data['B_vz'].append(float(parts[12]))
                data['C_x'].append(float(parts[13]))
                data['C_y'].append(float(parts[14]))
                data['C_z'].append(float(parts[15]))
                data['C_vx'].append(float(parts[16]))
                data['C_vy'].append(float(parts[17]))
                data['C_vz'].append(float(parts[18]))
            except ValueError:
                # Handle any conversion errors
                continue
    return data

def create_animation(data, method_name, output_filename):
    """
    Creates a 3D animation of the simulation and saves it as a GIF.

    """
    from matplotlib.animation import FuncAnimation
    import matplotlib

    # Ensure we can save as GIF
    matplotlib.use('Agg')

    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Define colors for each star
    colors = {
        'A': 'blue',
        'B': 'green',
        'C': 'red'
    }

    # Initialize the plot elements we want to animate
    line_A, = ax.plot([], [], [], label='Star A', color=colors['A'])
    line_B, = ax.plot([], [], [], label='Star B', color=colors['B'])
    line_C, = ax.plot([], [], [], label='Star C', color=colors['C'])
    point_A, = ax.plot([], [], [], marker='o', color=colors['A'])
    point_B, = ax.plot([], [], [], marker='o', color=colors['B'])
    point_C, = ax.plot([], [], [], marker='o', color=colors['C'])

    # Set plot limits
    # Compute min and max for all positions
    all_x = data['A_x'] + data['B_x'] + data['C_x']
    all_y = data['A_y'] + data['B_y'] + data['C_y']
    all_z = data['A_z'] + data['B_z'] + data['C_z']

    max_range = np.max([np.max(all_x) - np.min(all_x),
                        np.max(all_y) - np.min(all_y),
                        np.max(all_z) - np.min(all_z)]) / 2.0

    mid_x = (np.max(all_x) + np.min(all_x)) / 2.0
    mid_y = (np.max(all_y) + np.min(all_y)) / 2.0
    mid_z = (np.max(all_z) + np.min(all_z)) / 2.0

    ax.set_xlim((mid_x - max_range)/6,( mid_x + max_range)/6)
    ax.set_ylim((mid_y - max_range)/6, (mid_y + max_range)/6)
    ax.set_zlim((mid_z - max_range)/6, (mid_z + max_range)/6)

    ax.set_xlabel('X Position (m)')
    ax.set_ylabel('Y Position (m)')
    ax.set_zlabel('Z Position (m)')
    ax.set_title(f'{method_name} Method 3D Trajectories Animation')
    ax.legend()

    # Number of frames in the animation
    num_frames = 500
    skip_frames = max(1, num_frames // 1000)  # Adjust this number to control the number of frames
    frames = range(0, num_frames, skip_frames)

    def update(num):
        idx = frames[num]
        # Update the data for the lines
        line_A.set_data(data['A_x'][:idx], data['A_y'][:idx])
        line_A.set_3d_properties(data['A_z'][:idx])

        line_B.set_data(data['B_x'][:idx], data['B_y'][:idx])
        line_B.set_3d_properties(data['B_z'][:idx])

        line_C.set_data(data['C_x'][:idx], data['C_y'][:idx])
        line_C.set_3d_properties(data['C_z'][:idx])

        # Update the data for the points
        point_A.set_data(data['A_x'][idx - 1:idx], data['A_y'][idx - 1:idx])
        point_A.set_3d_properties(data['A_z'][idx - 1:idx])

        point_B.set_data(data['B_x'][idx - 1:idx], data['B_y'][idx - 1:idx])
        point_B.set_3d_properties(data['B_z'][idx - 1:idx])

        point_C.set_data(data['C_x'][idx - 1:idx], data['C_y'][idx - 1:idx])
        point_C.set_3d_properties(data['C_z'][idx - 1:idx])

        return line_A, line_B, line_C, point_A, point_B, point_C

    # Create the animation
    ani = FuncAnimation(fig, update, frames=len(frames), interval=30)

    # Save the animation as a GIF
    ani.save(output_filename, writer='pillow', fps=30)

    print(f'Animation saved as {output_filename}')
    
def create_2d_animation(data, method_name, output_filename):
    """
    Creates a 2D animation of the simulation and saves it as a GIF.

    """
    from matplotlib.animation import FuncAnimation
    import matplotlib

    # Ensure we can save as GIF
    matplotlib.use('Agg')

    fig, ax = plt.subplots(figsize=(10, 8))

    # Define colors for each star
    colors = {
        'A': 'blue',
        'B': 'green',
        'C': 'red'
    }

    # Initialize the plot elements we want to animate
    line_A, = ax.plot([], [], label='Star A', color=colors['A'])
    line_B, = ax.plot([], [], label='Star B', color=colors['B'])
    line_C, = ax.plot([], [], label='Star C', color=colors['C'])
    point_A, = ax.plot([], [], marker='o', color=colors['A'])
    point_B, = ax.plot([], [], marker='o', color=colors['B'])
    point_C, = ax.plot([], [], marker='o', color=colors['C'])

    # Set plot limits
    # Compute min and max for all positions
    all_x = data['A_x'] + data['B_x'] + data['C_x']
    all_y = data['A_y'] + data['B_y'] + data['C_y']

    max_range = np.max([np.max(all_x) - np.min(all_x),
                        np.max(all_y) - np.min(all_y)]) / 2.0

    mid_x = (np.max(all_x) + np.min(all_x)) / 2.0
    mid_y = (np.max(all_y) + np.min(all_y)) / 2.0

    # Adjust the scaling factor as needed
    scaling_factor = 5  # Adjust to zoom in or out

    ax.set_xlim(mid_x - max_range / scaling_factor, mid_x + max_range / scaling_factor)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)

    ax.set_xlabel('X Position (m)')
    ax.set_ylabel('Y Position (m)')
    ax.set_title(f'{method_name} Method 2D Trajectories Animation')
    ax.legend()
    ax.grid(True)

    # Number of frames in the animation
    num_frames = len(data['time'])
    skip_frames = max(1, num_frames // 1000)  # Adjust this number to control the number of frames
    frames = range(0, num_frames, skip_frames)

    def update(num):
        idx = frames[num]
        # Update the data for the lines
        line_A.set_data(data['A_x'][:idx], data['A_y'][:idx])
        line_B.set_data(data['B_x'][:idx], data['B_y'][:idx])
        line_C.set_data(data['C_x'][:idx], data['C_y'][:idx])

        # Update the data for the points
        point_A.set_data(data['A_x'][idx - 1:idx], data['A_y'][idx - 1:idx])
        point_B.set_data(data['B_x'][idx - 1:idx], data['B_y'][idx - 1:idx])
        point_C.set_data(data['C_x'][idx - 1:idx], data['C_y'][idx - 1:idx])

        return line_A, line_B, line_C, point_A, point_B, point_C

    # Create the animation
    ani = FuncAnimation(fig, update, frames=len(frames), interval=30)

    # Save the animation as a GIF
    ani.save(output_filename, writer='pillow', fps=30)

    print(f'2D animation saved as {output_filename}')

def main():
    # Path to the simulation output file
    rk4_file = 'ez_aquarii_output.txt'

    # Read the simulation data
    rk4_data = read_data(rk4_file)

    # Plot static trajectories (optional)
    plt.figure()
    plt.plot(rk4_data['A_x'], rk4_data['A_y'], label='Star A')
    plt.plot(rk4_data['C_x'], rk4_data['C_y'], label='Star C')
    plt.xlabel('X Position (m)')
    plt.ylabel('Y Position (m)')
    plt.title('Star Trajectories')
    plt.legend()
    plt.axis('equal')
    plt.show()

    # Create animation
    create_animation(rk4_data, 'RK4', 'ez_aquarii_animation.gif')
    create_2d_animation(rk4_data, 'RK4', 'ez_aquarii_2d_animation.gif')

if __name__ == "__main__":
    main()