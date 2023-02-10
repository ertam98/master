import csv
import matplotlib.pyplot as plt
import numpy as np
from helpfunctions import geometric_mean

def main():
    print('---GEOMETRIC MEAN CONE---')
    gm_times = print_stats('output_gm.stat')

    print('---EXPONENTIAL CONE---')
    exp_times = print_stats('output_exp.stat')

    x, y = list(), list()
    for key in gm_times.keys():
        if key in exp_times.keys():
            x.append(gm_times[key])
            y.append(exp_times[key])

    fig, ax = plt.subplots()
    ax.scatter(x, y)
    ax.set_xlabel('Geometric mean cone')
    ax.set_ylabel('Exponential cone')
    maxtime = max(max(x), max(y))
    line = np.arange(0, maxtime, 0.01)
    ax.plot(line, line)
    plt.savefig('scatter.png')

def print_stats(filename):
    nSolved = 0 # number of problems solved correctly
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

    print('Problems solved correctly: %d' %(nSolved))
    print('Geometic mean of time: %.6f' %(geometric_mean(times.values())))

    return times

main()