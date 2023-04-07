import numpy as np
import matplotlib.pyplot as plt

def cell_evolution(timesteps, initial_cells, phase_durations, max_divisions):
    """
    Simulate the evolution of the number of cells over time.

    Parameters:
    timesteps (int): The number of timesteps to simulate.
    initial_cells (int): The initial number of cells.
    phase_durations (dict): A dictionary containing the duration of each cell cycle phase (G1, S, G2, M).
    max_divisions (int): The maximum number of divisions a cell can undergo before dying.

    Returns:
    np.array: The number of cells at each timestep.
    """
    # Calculate the total cell cycle duration
    cell_cycle_duration = sum(phase_durations.values())

    # Initialize an array to store the number of cells at each timestep
    cells = np.zeros(timesteps)

    # Set the initial number of cells
    cells[0] = initial_cells

    # Create an array to store the remaining divisions for each cell
    remaining_divisions = np.full(initial_cells, max_divisions)

    for t in range(1, timesteps):
        print(t)
        print(cells[t-1])
        # Get the indices of cells that are dividing at this timestep
        dividing_cells = np.where(remaining_divisions > 0)

        # Update the remaining divisions for dividing cells
        remaining_divisions[dividing_cells] -= 1

        # Calculate the number of new cells produced by dividing cells
        new_cells = len(dividing_cells[0])

        # Update the number of cells
        cells[t] = cells[t - 1] + new_cells

        # Update the remaining_divisions array to include new cells
        remaining_divisions = np.concatenate((remaining_divisions, np.full(new_cells, max_divisions)))

    return cells

# Example usage
timesteps = 25
initial_cells = 10
phase_durations = {"G1": 10, "S": 5, "G2": 10, "M": 2}
max_divisions = 5

cell_count = cell_evolution(timesteps, initial_cells, phase_durations, max_divisions)

fig, ax = plt.subplots()
ax.plot(cell_count)
ax.set_xlabel("Time (timesteps)")
ax.set_ylabel("Number of cells")
plt.show()

