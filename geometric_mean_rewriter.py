import mosek
from rewriter_functions import getindexset
from helpfunctions import combine_lists

inf = 0.0

def geometric_mean_rewriter(importfolder, problemfile, exportfolder):
    #problemfile = "running_original.mps.gz"
    #folder = "benchmark"
    #problemfile = "markshare_4_0.mps.gz"

    with mosek.Task() as task:
        # read problemfile
        task.readdataformat(importfolder + "/" + problemfile, 
                            mosek.dataformat.mps, 
                            mosek.compresstype.gzip)

        # Get dimensions of original problem
        n = task.getnumvar()
        m = task.getnumcon()

        bkc, blc, buc = task.getconboundslice(0, m)
        bkx, blx, bux = task.getvarboundslice(0, n)

        I = getindexset(n, m, bkc, blc, buc, bkx, blx, bux)

        # parameter for size of cone
        p = sum(len(I[i]) for i in range(4))

        # add t as variable
        task.appendvars(1)

        # Build F for cone
        task.appendafes(p+1)
        currentrow = 0
        for i in I[0]:
            nzi, subi, vali = task.getarow(i)
            for k in range(nzi):
                # -rows in I[0] in F
                task.putafefentry(currentrow, subi[k], -vali[k])
            # u_c in g
            task.putafeg(currentrow, buc[i])
            currentrow += 1
        for i in I[1]:
            nzi, subi, vali = task.getarow(i)
            # +rows in I[1]
            task.putafefrow(currentrow, subi, vali)
            # -l_c in g
            task.putafeg(currentrow, -blc[i])
            currentrow += 1
        for j in I[2]:
            # -x_j in cone
            task.putafefentry(currentrow, j, -1.0)
            # u_x in g
            task.putafeg(currentrow, bux[j])
            currentrow += 1
        for j in I[3]:
            # x_j in cone
            task.putafefentry(currentrow, j, 1.0)
            # -l_x in g
            task.putafeg(currentrow, -blx[j])
            currentrow += 1
        # t in cone as last element. t has position n in vars
        task.putafefentry(p, n, 1.0)

        # Set <F, (x,t)> + g in GMcone
        gmdomain = task.appendprimalgeomeanconedomain(p+1)
        task.appendacc(gmdomain, range(p+1), None)

        # remove redundant constraints
        # combine lists to avoid duplicates
        task.removecons(combine_lists(I[0], I[1]))

        # remove redundant variable bounds
        # combine lists to avoid duplicates
        task.putvarboundlistconst(combine_lists(I[2], I[3]), mosek.boundkey.fr, -inf, +inf)

        # remove integrality to get LP-relaxation
        task.putvartypelist(range(n),
                            [mosek.variabletype.type_cont,]*n)

        # make t a free variable
        task.putvarbound(n, mosek.boundkey.fr, -inf, +inf)

        # change objective function
        for j in range(n):
            task.putcj(j, 0.0)
        task.putcj(n, 1.0)
        task.putobjsense(mosek.objsense.maximize)

        # export (remove .mps.gz from filename)
        task.writedata(exportfolder + "/" + problemfile[0:-7] + ".ptf.gz")