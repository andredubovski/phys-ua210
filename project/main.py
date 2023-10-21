import matplotlib.pyplot as plt
import numpy as np

from cells import Cell

cells = []

for i in range(10):
    rho_L = 10
    e_L = 4
    cells.append(Cell(rho_L, 0.1, e_L))

for i in range(10):
    rho_R = 1
    e_R = 5
    cells.append(Cell(rho_R, 0.1, e_R))

times = []
rho_vs_time = []
for i in range(len(cells)-2):
    rho_vs_time.append([])


for t in range(10):
    times.append(t)
    for i, cell in enumerate(cells):
        if(i > 0 and i < len(cells)-1):
            left_cell = cells[i-1]
            right_cell = cells[i+1]

            print("left cell F_half from main:", left_cell.F_half[0])
            print("cell F_half from main:", cell.F_half[0])
            cell.find_F_half(left_cell)
            cell.evolve(1, 1, left_cell)
            rho_vs_time[i-1].append(cell.U[0])
            # print("alpha plus:", cell.alpha_plus)
            # print("alpha minus:", cell.alpha_minus)

print(rho_vs_time)
rho_vs_x = np.transpose(rho_vs_time)

for line in rho_vs_x:
    plt.plot(range(len(line)), line)

plt.show()