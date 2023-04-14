import mosek
from MyTask import MyTask
from helpfunctions import streamprinter
from optimizer import createcsv
import csv
inf = 0.0

def main():
    importfolder = 'benchmark_presolved'
    exportfolder = 'benchmark_optimalface_org_presolve'
    statsfile = 'timing_optimalface_org_presolve.stat'

    header = ['Instance',
              'Response code',
              'Problem status',
              'Solution status',
              'Time',
              'Iterations',
              'Primal norm',
              'Dual norm']

    createcsv(statsfile, header)

    with open('benchmark-v2.txt', 'r') as file:
        i = 0
        for line in file:
            i += 1
            instance = line.replace('.mps.gz\n', '')

            with MyTask() as task:
                task.readdata(importfolder + '/' + instance + '.task.gz')

                # LP-relaxtion
                n = task.getnumvar()
                task.putvartypelist(range(n),[mosek.variabletype.type_cont,]*n)

                # Find basic solution
                task.putdouparam(mosek.dparam.optimizer_max_time, 100.0)
                rescode = task.optimize()

                # Restrict to optimal face
                task.restrictToOptimalFace()

                # Export
                task.writedata(exportfolder + '/' + instance + '.task.gz')
                print('%d: Exported the optimal face of %s' %(i, instance))

                task.updatesolutioninfo(mosek.soltype.bas)
                data = [instance,
                        rescode,
                        task.getintinf(mosek.iinfitem.sol_bas_prosta),
                        task.getintinf(mosek.iinfitem.sol_bas_solsta),
                        task.getdouinf(mosek.dinfitem.optimizer_time),
                        task.getintinf(mosek.iinfitem.intpnt_iter),
                        task.getdouinf(mosek.dinfitem.sol_bas_nrm_xx),
                        task.getdouinf(mosek.dinfitem.sol_bas_nrm_y)]

                with open(statsfile, 'a') as csvfile:
                    writer = csv.writer(csvfile)
                    writer.writerow(data)

if __name__ == '__main__':
    main()