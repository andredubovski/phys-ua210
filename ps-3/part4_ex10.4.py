import numpy as np
import matplotlib.pyplot as plt

# Initial number of Tl atoms and decay constant
initial_atoms_Tl = 1000
tau_Tl = 3.053 * 60
decay_constant = np.log(2) / tau_Tl

# Generate random decay times for Tl atoms
decay_times = -1 / decay_constant * np.log(1 - np.random.random(initial_atoms_Tl))
decay_times = np.sort(decay_times)

# Calculate the number of surviving atoms and decayed atoms at each time point
decay = np.arange(1, initial_atoms_Tl + 1)
survive = initial_atoms_Tl - decay

# Plot the decay curve
plt.plot(decay_times, survive)
plt.plot(decay_times, decay)
plt.xlabel("Time (s)")
plt.ylabel("N")
plt.legend(["Tl", "Pb"])
plt.savefig("part4_Decay")
plt.show()
