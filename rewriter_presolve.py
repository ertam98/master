import mosek
from MyTask import MyTask
from timeit import default_timer as timer

def main():
    with open('benchmark-v2-small.txt', 'r') as file:
        for line in file:
            instance = line.replace('.mps.gz\n', '')

            rewriter_presolve('benchmark_optimalface_org_presolve', instance, '.task.gz', 'benchmark_optimalface_org_presolve_presolved')

def rewriter_presolve(importfolder, instance, filetype, exportfolder):
    with MyTask() as task:
        task.readdata(importfolder + "/" + instance + filetype)

        try:
            print('Starting presolve')
            start = timer()
            task.presolve()
            end = timer()
            with open('timing_org_presolve_optimalface.stat', 'a') as file:
                file.write('%s, %f\n' %(instance, end-start))
            print('Finished presolve for %s in %f sek' %(instance, end-start))

            task.writedata(exportfolder + "/" + instance + '.task.gz')
        except Exception:
            with open('timing_org_presolve_optimalface.stat', 'a') as file:
                file.write('%s, %s\n' %(instance, 'infeasible'))
            print(instance, 'infeasible by presolve')

if __name__ == '__main__':
    main()