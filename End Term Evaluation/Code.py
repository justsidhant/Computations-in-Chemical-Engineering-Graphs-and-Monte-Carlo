# End term evaluation by Sidhant Thalor (221055) and Alankrit Gupta (220110)
import numpy as np
import matplotlib.pyplot as plt

volume = np.array([2.4, 1, 0.21])
total_molecules = 1e5
steps = 220
iterations = 38

mole_frac = np.zeros((steps + 1, 3))
vol_prod = np.zeros((steps, 3))
mole_frac[0] = np.array([0.85, 0.12, 0.03])

for step in range(1, steps + 1):
    molecules = total_molecules * mole_frac[step - 1]

    for _ in range(iterations):
        for component in range(3):
            vol_prod[step - 1, component] = volume[component] * molecules[component]

        max_product = np.max(vol_prod[step - 1])
        ratio = vol_prod[step - 1] / max_product

        random = np.random.random()  #random float between 0 and 1
        molecules = np.where(ratio > random, molecules - 1, molecules)

    mole_frac[step] = molecules / np.sum(molecules)

# Plotting the ternary diagram
plt.figure(figsize=(8, 6))
plt.scatter(mole_frac[:, 0], mole_frac[:, 1], c=mole_frac[:, 2], cmap='viridis', s=0.5)
plt.colorbar(label='Cumene Mole Fraction')
plt.xlabel('Benzene Mole Fraction')
plt.ylabel('Toluene Mole Fraction')
plt.title('Mole Fraction Distribution in Ternary System')
plt.grid(True)
plt.show()

# Plotting the mole fraction changes over steps
plt.figure(figsize=(8, 6))
plt.plot(range(steps + 1), mole_frac[:, 0], label='Benzene')
plt.plot(range(steps + 1), mole_frac[:, 1], label='Toluene')
plt.plot(range(steps + 1), mole_frac[:, 2], label='Cumene')
plt.xlabel('Simulation Steps')
plt.ylabel('Mole Fraction')
plt.title('Evolution of Mole Fractions Over Time')
plt.legend()
plt.grid(True)
plt.show()
