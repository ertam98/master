import mosek
from helpfunctions import streamprinter
import csv

def main():
    instances = list()
    with open('benchmark-v2-error-2.txt', 'r') as file:
        for line in file:
            instances.append(line[0:-8] + '.ptf.gz')

    header = ['Instance', 
              'Problem status',
              'Solution status',
              'Time',
              'Iterations',
              'Primal norm',
              'Dual norm']

    #createcsv('output_exp.stat', header)
    #createcsv('output_gm.stat', header)

    for instance in instances:
        solve_export('benchmark_exp', instance, 'exp')
        solve_export('benchmark_gm', instance, 'gm')
        
def createcsv(filename, header):
    with open(filename, 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)

def solve_export(folder, instance, model):
    with mosek.Task() as task:
        # task.set_Stream(mosek.streamtype.log, streamprinter)

        task.readdataformat(folder + '/' + instance,
                            mosek.dataformat.ptf,
                            mosek.compresstype.gzip)

        task.putdouparam(mosek.dparam.optimizer_max_time, 10.0)
        # turn off dualiser?
        rescode = task.optimize()
        task.updatesolutioninfo(mosek.soltype.itr)

        if rescode == mosek.rescode.ok:
            print('Solved %s for %s' %(instance, model))
            data = [instance,
                    task.getintinf(mosek.iinfitem.sol_itr_prosta),
                    # task.getprosta(mosek.soltype.itr),
                    task.getintinf(mosek.iinfitem.sol_itr_solsta),
                    # task.getsolsta(mosek.soltype.itr),
                    task.getdouinf(mosek.dinfitem.optimizer_time),
                    task.getintinf(mosek.iinfitem.intpnt_iter),
                    task.getdouinf(mosek.dinfitem.sol_itr_nrm_xx),
                    task.getdouinf(mosek.dinfitem.sol_itr_nrm_y)]

            with open('output_%s.stat' %(model), 'a') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(data)
        elif rescode == mosek.rescode.trm_max_time:
            print('%s for %s terminated due to max time' %(instance, model))
        else:
            print('Unexpected error in %s for %s: %s' %(instance, model, rescode))
            
main()