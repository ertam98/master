import csv
import matplotlib.pyplot as plt
import numpy as np
from helpfunctions import geometric_mean

def main():
    print('---GEOMETRIC MEAN CONE---')
    gm_times = print_stats('output_gm.stat')

    print('---EXPONENTIAL CONE---')
    exp_times = print_stats('output_exp.stat')

    max_gm = max(gm_times.values())
    max_exp = max(exp_times.values())
    maxtime = max(max_gm, max_exp)
    
    gm, exp = list(), list()
    for key in gm_times.keys():
        if key in exp_times.keys():
            gm.append(gm_times[key])
            exp.append(exp_times[key])
        else:
            gm.append(gm_times[key])
            exp.append(10)
    for key in exp_times.keys():
        if key not in gm_times.keys():
            gm.append(10)
            exp.append(exp_times[key])

    fig, ax = plt.subplots()
    ax.scatter(gm, exp)
    ax.set_xlabel('Geometric mean cone')
    ax.set_ylabel('Exponential cone')
    line = np.arange(0, maxtime, 0.01)
    ax.plot(line, line)
    plt.savefig('scatter.png')

def print_stats(filename):
    nSolved = 0 # number of problems solved correctly
    nPrimIll = 0
    nDualIll = 0
    nPriminf = 0
    nDualinf = 0
    times = dict()

    with open(filename, 'r') as csvfile:
        reader = csv.DictReader(csvfile)

        for row in reader:
            # prim and dual feasible and optimal
            if (int(row['Problem status']) == 1 and
                int(row['Solution status']) == 1):
                nSolved += 1
                times[row['Instance']] = float(row['Time'])
                # times.append(float(row['Time']))
            elif int(row['Solution status']) == 5:
                nPriminf += 1
            elif int(row['Solution status']) == 6:
                nDualinf += 1
            elif int(row['Solution status']) == 7:
                nPrimIll += 1
            elif int(row['Solution status']) == 8:
                nDualIll += 1

    print('Problems solved correctly: %d' %(nSolved))
    print('Geometic mean of time: %.6f' %(geometric_mean(times.values())))
    print('Primal infeasible: %d' %(nPriminf))
    print('Dual infeasible: %d' %(nDualinf))
    print('Primal illposed: %d' %(nPrimIll))
    print('Dual illposed: %d' %(nDualIll))

    return times

main()