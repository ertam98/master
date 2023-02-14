import csv
import matplotlib.pyplot as plt
import numpy as np
from helpfunctions import geometric_mean

def main():
    print('---GEOMETRIC MEAN CONE---')
    GMplot = MyPlotter()
    GMplot.getdata('output_gm.stat')
    print(GMplot)

    print('---EXPONENTIAL CONE---')
    Expplot = MyPlotter()
    Expplot.getdata('output_exp.stat')
    print(Expplot)

    GMplot.exportTimePlot(Expplot, 'Geometric Mean cone',
                         'Exponential cone', 'scatter.png')

class MyPlotter:
    def __init__(self):
        self.filename = str()
        self.nSolved = 0 # number of problems solved correctly
        self.nMaxtime = 0 # number of problems timed out
        self.nPrimIll = 0
        self.nDualIll = 0
        self.nPriminf = 0
        self.nDualinf = 0
        self.nUnknown = 0
        self.times = dict()
        self.iterations = dict()

    def getdata(self, filename):
        with open(filename, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                # prim and dual feasible and optimal
                if (int(row['Problem status']) == 1 and
                    int(row['Solution status']) == 1):
                    self.nSolved += 1
                    self.times[row['Instance']] = float(row['Time'])
                    self.iterations[row['Instance']] = int(row['Iterations'])
                    # times.append(float(row['Time']))
                elif int(row['Solution status']) == 5:
                    self.nPriminf += 1
                elif int(row['Solution status']) == 6:
                    self.nDualinf += 1
                elif int(row['Solution status']) == 7:
                    self.nPrimIll += 1
                elif int(row['Solution status']) == 8:
                    self.nDualIll += 1
                elif row['Response code'] == 'rescode.trm_max_time':
                    self.nMaxtime += 1
                elif (int(row['Problem status']) == 0 and
                    int(row['Solution status']) == 0):
                    self.nUnknown += 1

    def __str__(self):
        string = 'Problems solved correctly: %d\n' %(self.nSolved)
        string += 'Geometic mean of time: %.6f\n' %(geometric_mean(self.times.values()))
        string += 'Problems timed out: %d\n' %(self.nMaxtime)
        string += 'Primal infeasible: %d\n' %(self.nPriminf)
        string += 'Dual infeasible: %d\n' %(self.nDualinf)
        string += 'Primal illposed: %d\n' %(self.nPrimIll)
        string += 'Dual illposed: %d\n' %(self.nDualIll)
        string += 'Unknown status: %d' %(self.nUnknown)
        return string

    def exportTimePlot(self, other, selflabel, otherlabel, filename):
        max_self = max(self.times.values())
        max_other = max(other.times.values())
        max_time = max(max_self, max_other)

        solved_both = {'self': list(), 'other': list()}
        solved_self = {'self': list(), 'other': list()}
        solved_other = {'self': list(), 'other': list()}
        for key in self.times.keys():
            if key in other.times.keys():
                solved_both['self'].append(self.times[key])
                solved_both['other'].append(other.times[key])
            else:
                solved_self['self'].append(self.times[key])
                solved_self['other'].append(max_time)
        for key in other.times.keys():
            if key not in self.times.keys():
                solved_other['self'].append(max_time)
                solved_other['other'].append(other.times[key])

        _, ax = plt.subplots()
        ax.scatter(solved_both['self'], solved_both['other'], c = 'tab:blue')
        ax.scatter(solved_self['self'], solved_self['other'], c = 'tab:red')
        ax.scatter(solved_other['self'], solved_other['other'], c = 'tab:green')
        ax.set_xlabel(selflabel)
        ax.set_ylabel(otherlabel)
        line = np.arange(0, max_time, 0.01)
        ax.plot(line, line)
        plt.savefig(filename)

if __name__ == '__main__':
    main()