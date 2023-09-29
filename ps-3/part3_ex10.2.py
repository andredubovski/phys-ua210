import numpy as np
import matplotlib.pyplot as plt
import random

t_max = 20000
N_213Bi = 10000
step_s = 1
N_209Bi, N_Pb, N_Tl = 0, 0, 0
_213Bi_points, Pb_points, Tl_points, _209Bi_points = [], [], [], []

for t in np.arange(0, t_max, step_s):
    tau_Pb = 3.3 * 60
    p_Pb = 1 - 2 ** (-step_s / tau_Pb)
    tau_Tl = 2.2 * 60
    p_Tl = 1 - 2 ** (-step_s / tau_Tl)
    tau_213Bi = 46 * 60
    p_Bi = 1 - 2 ** (-step_s / tau_213Bi)
    p_Pb_Bi = p_Bi * 0.9791
    p_Tl_Bi = p_Bi * 0.0209

    Pb_points.append(N_Pb)
    Tl_points.append(N_Tl)
    _209Bi_points.append(N_209Bi)
    _213Bi_points.append(N_213Bi)

    decay_pb = 0
    for i in range(N_Pb):
        if random.random() < p_Pb:
            decay_pb += 1
    N_Pb -= decay_pb
    N_209Bi += decay_pb

    decay_Tl = 0
    for i in range(N_Tl):
        if random.random() < p_Tl:
            decay_Tl += 1
    N_Tl -= decay_Tl
    N_Pb += decay_Tl

    decay_to_Pb = 0
    decay_to_Tl = 0
    for i in range(N_213Bi):
        if random.random() < p_Tl_Bi:
            decay_to_Tl += 1
        elif p_Tl_Bi <= random.random() < (p_Pb_Bi + p_Tl_Bi):
            decay_to_Pb += 1
    N_Pb += decay_to_Pb
    N_Tl += decay_to_Tl
    N_213Bi -= (decay_to_Pb + decay_to_Tl)

plt.plot(_213Bi_points)
plt.plot(Pb_points)
plt.plot(Tl_points)
plt.plot(_209Bi_points)
plt.xlabel("Time (s)")
plt.ylabel("N")
plt.legend(["213Bi", "Pb", "Tl", "209Bi"])
plt.savefig("part3_radioactive_decay")
plt.show()
