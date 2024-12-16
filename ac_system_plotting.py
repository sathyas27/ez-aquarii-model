import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

def read_data(filename):
    """
    Reads the simulation data from a file and organizes it into dictionaries.
    """
    data = {
        'time': [],
        'A_x': [],
        'A_y': [],
        'A_z': [],
        'A_vx': [],
        'A_vy': [],
        'A_vz': [],
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
            parts = line.strip().split('\t')
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
        'C': 'red'
    }

    # Initialize the plot elements we want to animate
    line_A, = ax.plot([], [], [], label='Star A', color=colors['A'])
    line_C, = ax.plot([], [], [], label='Star C', color=colors['C'])
    point_A, = ax.plot([], [], [], marker='o', color=colors['A'])
    point_C, = ax.plot([], [], [], marker='o', color=colors['C'])

    # Set axis limits for the animation plot
    ax.set_xlim(-1e10, 1e10)
    ax.set_ylim(-1e10, 1e10)
    ax.set_zlim(-1e10, 1e10) 
    
    ax.set_xlabel('X Position (m)')
    ax.set_ylabel('Y Position (m)')
    ax.set_zlabel('Z Position (m)')
    ax.set_title(f'{method_name} Method 3D Trajectories Animation')
    ax.legend()

    # Number of frames in the animation
    num_frames = len(data['time'])
    skip_frames = 100  
    frames = range(0, num_frames, skip_frames)

    def update(num):
        idx = frames[num]
        # Update the data for the lines
        line_A.set_data(data['A_x'][:idx], data['A_y'][:idx])
        line_A.set_3d_properties(data['A_z'][:idx])


        line_C.set_data(data['C_x'][:idx], data['C_y'][:idx])
        line_C.set_3d_properties(data['C_z'][:idx])

        # Update the data for the points
        point_A.set_data(data['A_x'][idx - 1:idx], data['A_y'][idx - 1:idx])
        point_A.set_3d_properties(data['A_z'][idx - 1:idx])


        point_C.set_data(data['C_x'][idx - 1:idx], data['C_y'][idx - 1:idx])
        point_C.set_3d_properties(data['C_z'][idx - 1:idx])

        return line_A, line_C, point_A, point_C

    # Create the animation
    ani = FuncAnimation(fig, update, frames=len(frames), interval=30, blit=True)

    # Save the animation as a GIF
    ani.save(output_filename, writer='pillow', fps=30)

    print(f'Animation saved as {output_filename}')

def main():
    # Paths to the simulation output files
    rk4_file = 'ac_output.txt'
    
    # Read the simulation data
    rk4_data = read_data(rk4_file)
    
    # Create animations
    create_animation(rk4_data, 'RK4', 'ac_system_animate.gif')

if __name__ == "__main__":
    main()
    
    
    #bbbbbbbbbbbbdjshfkljasdhfasipdjfhpaisdhfphasdkj;fhh;jkasdf;jhsad;lfh;kjasdhfkljahsdfkjaskdhflkjasdhflkjashdfljkasdhfhlkjasdfjhlaksjdhfkjlasdhfkjalsdfhjlkasdfhkalsdf