# max x_2
# s.t. x_1 + 2*x_2 <= 6
#      x_1 >= -2
#      x_2 >= 1

import mosek
from helpfunctions import streamprinter

inf = 0.0

def main():
    with mosek.Task() as task:
        task.set_Stream(mosek.streamtype.log, streamprinter)

        numvar = 2
        numcon = 1

        task.appendvars(numvar)
        task.appendcons(numcon)

        task.putcj(1, 1.0)

        bkx = [mosek.boundkey.lo, mosek.boundkey.lo]
        blx = [-2.0, 1.0]
        bux = [+inf, +inf]
        task.putvarboundlist([0,1], bkx, blx, bux)

        task.putaijlist([0,0], [0,1], [1.0, 2.0])
        task.putconbound(0, mosek.boundkey.up, -inf, 6.0)

        task.putobjsense(mosek.objsense.maximize)
        task.optimize()
        task.optimizersummary(mosek.streamtype.msg)

        solsta = task.getsolsta(mosek.soltype.itr)

        if solsta == mosek.solsta.optimal:
            xx = task.getxx(mosek.soltype.itr)
            print("Optimal solution: ")
            for i in range(numvar):
                print("x[%d]=%.2f" %(i, xx[i]))

        elif (solsta == mosek.solsta.dual_infeas_cer or
             solsta == mosek.solsta.prim_infeas_cer):
             print("Primal or dual infeasibility certificate found.\n")
        elif solsta == mosek.solsta.unknown:
            print("Unkown solution status")
        else:
            print("Other solution status")


        #task.writedata("running_original.mps.gz")
main()