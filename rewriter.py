import mosek
from helpfunctions import combine_lists
from MyTask import MyTask
from timeit import default_timer as timer


def main():
    with open('benchmark-v2.txt', 'r') as file:
        for line in file:
            instance = line.replace('.mps.gz\n', '')

            rewriter('exp', 'benchmark_presolved', instance, '.task.gz', 'benchmark_exp_presolve_1', False)
            rewriter('gm', 'benchmark_presolved', instance, '.task.gz', 'benchmark_gm_presolve_1', False)


def rewriter(model, importfolder, instance, fileformat, exportfolder, withpresolve):
    if not (model == 'gm' or model == 'exp'):
        raise ValueError("Parameter 'model' needs to be 'gm' or 'exp'.")

    inf = 0.0

    with MyTask() as task:
        task.readdata(importfolder + "/" + instance + fileformat)
        #task.readdataformat(importfolder + "/" + problemfile, 
        #                    mosek.dataformat.mps, 
        #                    mosek.compresstype.gzip)

        if withpresolve:
            try:
                print('Starting presolve')
                start = timer()
                task.presolve()
                # task.removeemptyrows()
                # task.presolve_domain(1e8)
                # task.remove_redundant(1e-6)
                # task.presolve_lindep()
                end = timer()
                with open('timing_presolve_%s.stat' %(model), 'a') as file:
                    file.write('%s, %f\n' %(instance, end-start))
                print('Finished presolve in %f sek' %(end-start))
            except Exception:
                with open('timing_presolve_%s.stat' %(model), 'a') as file:
                    file.write('%s, %s\n' %(instance, 'infeasible'))
                print(instance, 'infeasible by presolve')

        n = task.getnumvar() # vars in original problem

        task.getindexset()
        task.updatep()

        if model == 'gm':
            task.buildgmcone()
            
            # change objective function
            for j in range(n):
                task.putcj(j, 0.0)
            task.putcj(n, 1.0)
            task.putobjsense(mosek.objsense.maximize)

            # make t a free variable
            task.putvarbound(n, mosek.boundkey.fr, -inf, +inf)
        elif model == 'exp':
            task.buildexpcone()

            # change objective function
            for j in range(n):
                task.putcj(j, 0.0)
            for j in range(n, n+task.p):
                task.putcj(j, 1.0)
            task.putobjsense(mosek.objsense.maximize)

            # make all t's free variables
            task.putvarboundsliceconst(n, n+task.p,
                                       mosek.boundkey.fr, -inf, +inf)

        # remove redundant constraints
        # combine lists to avoid duplicates
        # I[6] are rows which are 0
        task.removecons(combine_lists(task.I[0], task.I[1], task.I[6]))

        # remove redundant variable bounds
        # combine lists to avoid duplicates
        task.putvarboundlistconst(combine_lists(task.I[2], task.I[3]),
                                  mosek.boundkey.fr, -inf, +inf)

        # remove integrality to get LP-relaxation
        task.putvartypelist(range(n), [mosek.variabletype.type_cont,]*n)

        # export (remove .mps.gz from filename)
        task.writedata(exportfolder + "/" + instance + '.task.gz')
        #               + problemfile.replace('.mps.gz', '.task.gz')) 
        #task.writedata(exportfolder + "/" + problemfile[0:-7] + ".task.gz")      

    print('Sucessfully exported %s to %s for %s' %(instance, exportfolder,
                                                   model))

if __name__ == '__main__':
    main()