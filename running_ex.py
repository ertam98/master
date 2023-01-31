# Analytical centre of:
# x_1 >= -2
# x_2 >= 1
# x_1 + 2*x_2 <= 6
# Ans: (0, 2)

import mosek
import sys

inf = 0.0

def streamprinter(text):
    sys.stdout.write(text)
    sys.stdout.flush()

def main():
    with mosek.Task() as task:
        task.set_Stream(mosek.streamtype.log, streamprinter)

        numvar = 2
        numcon = 1
        p = 2*numvar + numcon

        task.appendvars(p) # num of x + num of t
        task.appendcons(numcon)

        c = [0.0, 0.0, 1.0, 1.0, 1.0]
        task.putcslice(0, 5, c)

        for i in range(p):
            task.putvarbound(i, mosek.boundkey.fr, -inf, inf)

        # A = [1, 2]
        # task.putaij(0, 0, 1)
        # task.putaij(0, 1, 2)
        
        # Ax <= 6
        # task.putconbound(0, mosek.boundkey.up, -inf, 6.0)

        numafe = 9
        numcone = 3
        task.appendafes(numafe)
        
        # Build F
        Fsub = [[0, 6],
                [3, 6],
                [2,],
                [5,],
                [8,]]
        Fval = [[1.0, -1.0],
                [1.0, -2.0],
                [1.0,],
                [1.0,],
                [1.0,]]

        for j in range(p):
            task.putafefcol(j, Fsub[j], Fval[j])

        # Build g
        gsub = [0, 1, 3, 4, 6, 7]
        gval = [2.0, 1.0, -1.0, 1.0, 6.0, 1.0]
        task.putafeglist(gsub, gval)

        expdomain = task.appendprimalexpconedomain()
        for i in range(0, 9, 3):
            task.appendacc(expdomain, list(range(i, i+3)), None)
            
        task.putobjsense(mosek.objsense.maximize)

        task.optimize()
        task.solutionsummary(mosek.streamtype.msg)
        solsta = task.getsolsta(mosek.soltype.itr)

        if solsta == mosek.solsta.optimal:
            xx = task.getxx(mosek.soltype.itr)
            print("Optimal solution: ")
            for i in range(p):
                print("x[%d]=%.2f" %(i, xx[i]))

        elif (solsta == mosek.solsta.dual_infeas_cer or
             solsta == mosek.solsta.prim_infeas_cer):
             print("Primal or dual infeasibility certificate found.\n")
        elif solsta == mosek.solsta.unknown:
            print("Unkown solution status")
        else:
            print("Other solution status")
main()