import mosek
from helpfunctions import streamprinter
from rewriter_functions import getindexset

inf = 0.0

def main():
    problemfile = "running_original.mps.gz"

    with mosek.Task() as task:
        # read problemfile
        task.readdataformat(problemfile,
                            mosek.dataformat.mps,
                            mosek.compresstype.gzip)

        # Get dimensions of original problem
        n = task.getnumvar()
        m = task.getnumcon()

        bkc, blc, buc = task.getconboundslice(0, m)
        bkx, blx, bux = task.getvarboundslice(0, n)

        # Calculate index sets
        I = getindexset(n, m, bkc, blc, buc, bkx, blx, bux)

        # parameter for number of t's to add
        p = sum(len(I[i]) for i in range(4))

        # Add t's
        task.appendvars(p)

        # Build F for cones
        task.appendafes(3*p)
        # First build rows corresponding to first element in each exp cone
        currentrow = 0
        for i in I[0]:
            nzi, subi, vali = task.getarow(i)
            for k in range(nzi):
                # -rows in I[0] in F
                task.putafefentry(currentrow, subi[k], -vali[k])
            # u_c in g
            task.putafeg(currentrow, buc[i])
            currentrow += 3
        for i in I[1]:
            nzi, subi, vali = task.getarow(i)
            # +rows in I[1]
            task.putafefrow(currentrow, subi, vali)
            # -l_c in g
            task.putafeg(currentrow, -blc[i])
            currentrow += 3
        for j in I[2]:
            # -x_j in cone
            task.putafefentry(currentrow, j, -1.0)
            # u_x in g
            task.putafeg(currentrow, bux[j])
            currentrow += 3
        for j in I[3]:
            # x_j in cone
            task.putafefentry(currentrow, j, 1.0)
            # -l_x in g
            task.putafeg(currentrow, -blx[j])
            currentrow += 3
            
        # 2nd element in each cone is 1
        task.putafeglist(range(1, 3*p, 3), [1.0,]*p)

        # 3rd element in each cone is t_i. t_i has position n+i in vars
        task.putafefentrylist(range(2, 3*p, 3), range(n, n+p), [1.0]*p)
           
        expdomain = task.appendprimalexpconedomain()
        for cone in range(p):
            task.appendacc(expdomain, range(3*cone, 3*cone+3), None)
            
        # remove redudant constraints
        task.removecons(I[0]+I[1])

        # remove redundant variable bounds
        task.putvarboundlistconst(I[2]+I[3], mosek.boundkey.fr, -inf, +inf)

        # make all t's free variables
        task.putvarboundsliceconst(n, n+p, mosek.boundkey.fr, -inf, +inf)
            
        # change objective function
        for j in range(n):
            task.putcj(j, 0.0)
        for j in range(n, n+p):
            task.putcj(j, 1.0)
        task.putobjsense(mosek.objsense.maximize)

        # export
        task.writedata(problemfile[0:-7] + "_exp.ptf")

main()