import mosek
from MyTask import MyTask
from helpfunctions import streamprinter
inf = 0.0

def main():
    with MyTask() as task:
        task.set_Stream(mosek.streamtype.log, streamprinter)

        # import instance
        #task.readdata('benchmark/30n20b8.mps.gz')
        task.appendvars(2)
        task.appendcons(2)
        c = [3.0, 1.5]
        #A = [[1.0, 2.0], [2.0, 1.0]]
        sub = [0, 1]
        ptrb = [0, 2]
        ptre = [2, 4]
        asub = [0, 1, 0, 1]
        aval = [1.0, 2.0, 2.0, 1.0]
        
        bkx = [mosek.boundkey.lo, mosek.boundkey.lo]
        blx = [0.0, 0.0]
        bux = [inf, inf]
        bkc = [mosek.boundkey.up, mosek.boundkey.up]
        blc = [-inf, -inf]
        buc = [8.0, 7.0]

        task.putclist([0, 1], c)
        task.putarowlist(sub, ptrb, ptre, asub, aval)
        task.putvarboundslice(0, 2, bkx, blx, bux)
        task.putconboundslice(0, 2, bkc, blc, buc)

        task.putobjsense(mosek.objsense.maximize)
        
        # LP-relaxtion
        n = task.getnumvar()
        task.putvartypelist(range(n), [mosek.variabletype.type_cont,]*n)

        # find basic solution
        task.putdouparam(mosek.dparam.optimizer_max_time, 100.0)
        task.optimize()
        task.solutionsummary(mosek.streamtype.log)
        
        # restric to optimal face
        task.restrictToOptimalFace()
        task.writedata('debug/test_optface.ptf')


if __name__ == '__main__':
    main()