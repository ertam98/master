import mosek
from helpfunctions import streamprinter
import csv

def main():
    instances = list()
    with open('benchmark-v2.txt', 'r') as file:
        for line in file:
            instances.append(line.replace('.mps.gz\n', '.task.gz'))
            #instances.append(line[0:-8] + '.task.gz')

    header = ['Instance',
              'Response code',
              'Problem status',
              'Solution status',
              'Time',
              'Iterations',
              'Primal norm',
              'Dual norm']

    createcsv('output_exp_presolve.stat', header)
    createcsv('output_gm_presolve.stat', header)

    for instance in instances:
        solve_export('benchmark_exp_presolve_1', instance, 'exp')
        solve_export('benchmark_gm_presolve_1', instance, 'gm')
        
def createcsv(filename, header):
    with open(filename, 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)

def solve_export(folder, instance, model):
    with mosek.Task() as task:
        # task.set_Stream(mosek.streamtype.log, streamprinter)

        task.readdataformat(folder + '/' + instance,
                            mosek.dataformat.task,
                            mosek.compresstype.gzip)

        task.putdouparam(mosek.dparam.optimizer_max_time, 100.0)
        task.putintparam(mosek.iparam.intpnt_solve_form,
                         mosek.solveform.primal) # solve on primal form
        rescode = task.optimize()
        task.updatesolutioninfo(mosek.soltype.itr)

        print('%s for %s terminated with response code: %s' %(instance, model, rescode))

        data = [instance,
                rescode,
                task.getintinf(mosek.iinfitem.sol_itr_prosta),
                task.getintinf(mosek.iinfitem.sol_itr_solsta),
                task.getdouinf(mosek.dinfitem.optimizer_time),
                task.getintinf(mosek.iinfitem.intpnt_iter),
                task.getdouinf(mosek.dinfitem.sol_itr_nrm_xx),
                task.getdouinf(mosek.dinfitem.sol_itr_nrm_y)]

        with open('output_%s_presolve.stat' %(model), 'a') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(data)

if __name__ == '__main__':     
    main()