import mosek
from helpfunctions import combine_lists
from MyTask import MyTask

def main():
    with open('benchmark-v2.txt', 'r') as file:
        for line in file:
            instance = line[0:-1] # remove '\n'

            rewriter('exp', 'benchmark', instance, 'benchmark_exp')
            rewriter('gm', 'benchmark', instance, 'benchmark_gm')

    # for instance in instances:
    #     exponential_rewriter("benchmark", instance, "benchmark_exp")
    #     geometric_mean_rewriter("benchmark", instance, "benchmark_gm")

    # rewriter('exp', 'benchmark', "markshare_4_0.mps.gz", 'debug_exp')
    # rewriter('gm', 'benchmark', "markshare_4_0.mps.gz", 'debug_gm')
    
    #problemfile = "running_original.mps.gz"
    #importfolder = "benchmark"
    #problemfile = "markshare_4_0.mps.gz"

def rewriter(model, importfolder, problemfile, exportfolder):
    if not (model == 'gm' or model == 'exp'):
        raise ValueError("Parameter 'model' needs to be 'gm' or 'exp'.")

    inf = 0.0

    with MyTask() as task:
        task.readdataformat(importfolder + "/" + problemfile, 
                            mosek.dataformat.mps, 
                            mosek.compresstype.gzip)

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
        task.writedata(exportfolder + "/" + problemfile[0:-7] + ".ptf.gz")

    print('Sucessfully exported %s to %s for %s' %(problemfile, exportfolder, model))

main()