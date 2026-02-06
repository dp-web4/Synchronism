"""
Intent_Transfer_Models.py
Synchronism Framework Core Module - Planck-scale intent transfer simulations
"""

import numpy as np
from typing import Tuple

class PlanckGrid3D:
    def __init__(self, dimensions: Tuple[int, int, int] = (32, 32, 32)):
        """
        Initialize 3D Planck-scale grid with quantum intent values
        :param dimensions: Tuple (x, y, z) dimensions of the grid
        """
        self.dimensions = dimensions
        self.grid = np.random.randint(0, 4, dimensions, dtype=np.uint8)
        self.tension_field = np.zeros(dimensions, dtype=np.float32)
        
    def calculate_tension(self) -> None:
        """Calculate tension tensor based on neighboring cell intent differences"""
        padded_grid = np.pad(self.grid, 1, mode='wrap')
        
        # Calculate differences in all 6 directions (3D)
        diffs = [
            padded_grid[2:, 1:-1, 1:-1] - padded_grid[:-2, 1:-1, 1:-1],  # X-axis
            padded_grid[1:-1, 2:, 1:-1] - padded_grid[1:-1, :-2, 1:-1],  # Y-axis
            padded_grid[1:-1, 1:-1, 2:] - padded_grid[1:-1, 1:-1, :-2],  # Z-axis
        ]
        
        # Sum absolute differences for tension magnitude
        self.tension_field = np.sum(np.abs(diffs), axis=0) / 6

    def transfer_intent(self) -> None:
        """Perform intent transfer between adjacent cells based on tension"""
        new_grid = self.grid.copy()
        
        # Iterate through all cells
        for x in range(self.dimensions[0]):
            for y in range(self.dimensions[1]):
                for z in range(self.dimensions[2]):
                    current_intent = self.grid[x, y, z]
                    
                    # Get neighbors (wrapping at boundaries)
                    neighbors = [
                        self.grid[(x-1)%self.dimensions[0], y, z],
                        self.grid[(x+1)%self.dimensions[0], y, z],
                        self.grid[x, (y-1)%self.dimensions[1], z],
                        self.grid[x, (y+1)%self.dimensions[1], z],
                        self.grid[x, y, (z-1)%self.dimensions[2]],
                        self.grid[x, y, (z+1)%self.dimensions[2]],
                    ]
                    
                    # Calculate intent flow
                    for n in neighbors:
                        if current_intent > n:
                            transfer = min((current_intent - n) // 4, 1)
                            new_grid[x, y, z] -= transfer
                            new_grid[x, y, z] = max(new_grid[x, y, z], 0)
        
        self.grid = np.clip(new_grid, 0, 3)

    def calculate_coherence(self) -> float:
        """Calculate grid coherence score (0-1)"""
        unique, counts = np.unique(self.grid, return_counts=True)
        return np.max(counts) / np.sum(counts) if counts.size > 0 else 0

    def tick(self) -> None:
        """Advance simulation by one Planck time unit"""
        self.calculate_tension()
        self.transfer_intent()

if __name__ == "__main__":
    # Example usage
    grid = PlanckGrid3D((16, 16, 16))
    
    print("Initial coherence:", grid.calculate_coherence())
    for _ in range(10):
        grid.tick()
    print("Post-evolution coherence:", grid.calculate_coherence())
    
    # Basic visualization (requires matplotlib)
    try:
        import matplotlib.pyplot as plt
        plt.imshow(grid.grid[:, :, 8])  # Show middle Z-layer
        plt.title("Intent Distribution Slice")
        plt.colorbar(label="Intent Level")
        plt.show()
    except ImportError:
        print("Matplotlib not installed - visualization skipped")
