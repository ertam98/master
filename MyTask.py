import mosek
from helpfunctions import myhash
from math import isclose

class MyTask(mosek.Task):
    def __init__(self):
        super().__init__()
        self.I = [list() for i in range(7)] # index sets
        self.p = int() # parameter for cone
        self.neginf = float('-inf')
        self.posinf = float('inf')

    def getarownumnzlist(self, start, stop):
        return [self.getarownumnz(i) for i in range(start, stop)]

    def updatep(self):
        self.p = sum(len(self.I[i]) for i in range(4))

    def getindexset(self):
        n = self.getnumvar()
        m = self.getnumcon()

        bkc, _, _ = self.getconboundslice(0, m)
        bkx, _, _ = self.getvarboundslice(0, n)
        nz = self.getarownumnzlist(0, m)

        for i in range(m):
            if nz[i] == 0: # rows of A that are all zero
                self.I[6].append(i)
            else:
                boundtype = bkc[i]
                if boundtype == mosek.boundkey.ra:
                    self.I[0].append(i)
                    self.I[1].append(i)
                if boundtype == mosek.boundkey.up:
                    self.I[0].append(i)
                elif boundtype == mosek.boundkey.lo:
                    self.I[1].append(i)
                elif boundtype == mosek.boundkey.fx:
                    self.I[4].append(i)
        
        for j in range(n):
            boundtype = bkx[j]
            if boundtype == mosek.boundkey.ra:
                self.I[2].append(j)
                self.I[3].append(j)
            elif boundtype == mosek.boundkey.up:
                self.I[2].append(j)
            elif boundtype == mosek.boundkey.lo:
                self.I[3].append(j)
            elif boundtype == mosek.boundkey.fx:
                self.I[5].append(j)

    def presolve_domain(self):
        for j in range(self.getnumvar()): # for each variable
            nzj, subj, valj = self.getacol(j) # get all non-zero coeffs
            
            for count in range(nzj): # for each non-zero coeffs
                i = subj[count] # constraint index
                
                bkc, blc, buc = self.getconbound(i) # get bounds
                bkx, blx, bux = self.getvarbound(j) # updates each iteration

                if bkx == mosek.boundkey.fx: # already fixed
                    break

                inf, sup = self.getinfsuprow(i, j) # calculate inf and sup

                # if boundupdate is well-defined
                if self.__isufinite(bkc) and inf != self.neginf:
                    newbound = (buc - inf)/valj[count]
                    if valj[count] > 0.0: # if coeff is positive update upper
                        if self.__isufinite(bkx):
                            bux = min(newbound, bux)
                        else:
                            bux = newbound
                        self.chgvarbound(j, 0, 1, bux)
                    else: # if coeff is negative update lower
                        if self.__islfinite(bkx):
                            blx = max(newbound, blx)
                        else:
                            blx = newbound
                        self.chgvarbound(j, 1, 1, blx)

                bkx, blx, bux = self.getvarbound(j)

                # if boundupdate is well-defined
                if self.__islfinite(bkc) and sup != self.posinf:
                    newbound = (blc - sup)/valj[count]
                    if valj[count] > 0: # if coeff is positive update lower
                        if self.__islfinite(bkx):
                            blx = max(newbound, blx)
                        else:
                            blx = newbound
                        self.chgvarbound(j, 1, 1, blx)
                    else: # if coeff is negative update upper
                        if self.__isufinite(bkx):
                            bux = min(newbound, bux)
                        else:
                            bux = newbound
                        self.chgvarbound(j, 0, 1, bux)

    def __isufinite(self, bk):
        return bk in [mosek.boundkey.ra, mosek.boundkey.up, mosek.boundkey.fx]

    def __islfinite(self, bk):
        return bk in [mosek.boundkey.ra, mosek.boundkey.lo, mosek.boundkey.fx]

    def getinfsuprow(self, i, k):
        nzi, subi, vali = self.getarow(i)
        bkx, blx, bux = self.getvarboundlist(subi)

        if mosek.boundkey.fr in bkx:
            return self.neginf, self.posinf
        else:
            inf = self.__getinfrow(k, nzi, subi, vali, bkx, blx, bux)
            sup = self.__getsuprow(k, nzi, subi, vali, bkx, blx, bux)

        return inf, sup
    def __getinfrow(self, k, nzi, subi, vali, bkx, blx, bux):
        inf = 0.0

        for index in range(nzi):
            if subi[index] == k:
                pass
            elif vali[index] > 0.0:
                if bkx[index] == mosek.boundkey.up:
                    return self.neginf # x_j not bounded from below
                else:
                    inf += vali[index] * blx[index]
            else:
                if bkx[index] == mosek.boundkey.lo:
                    return self.neginf
                else:
                    inf += vali[index] * bux[index]

        return inf

    def __getsuprow(self, k, nzi, subi, vali, bkx, blx, bux):
        sup = 0.0

        for index in range(nzi):
            if subi[index] == k:
                pass
            elif vali[index] > 0.0:
                if bkx[index] == mosek.boundkey.lo:
                    return self.posinf
                else:
                    sup += vali[index] * bux[index]
            else:
                if bkx[index] == mosek.boundkey.up:
                    return self.posinf
                else:
                    sup += vali[index] * blx[index]

        return sup

    def getvarboundlist(self, subj):
        bkx, blx, bux = list(), list(), list()
        for j in subj:
            bk, bl, bu = self.getvarbound(j)
            bkx.append(bk)
            blx.append(bl)
            bux.append(bu)

        return bkx, blx, bux

    def getconboundlist(self, subi):
        bkc, blc, buc = list(), list(), list()
        for i in subi:
            bk, bl, bu = self.getconbound(i)
            bkc.append(bk)
            blc.append(bl)
            buc.append(bu)

        return bkc, blc, buc

    def presolve_lindep(self):
        m = self.getnumcon()
        ptrb, ptre, sub, val = self.getarowslice(0, m)
        rowpattern = dict() # dictionary for rows with same non-zero elements
        for i in range(m):
            subi = sub[ptrb[i]:ptre[i]]
            vali = val[ptrb[i]:ptre[i]]
            normvali = [element / vali[0] for element in vali]
            key = myhash(subi, normvali)
            if key in rowpattern.keys():
                rowpattern[key].append(i)
            else:
                rowpattern[key] = [i,]

        deletedrows = list()
        for rows in rowpattern.values():
            i = 0
            while i < len(rows)-1:
                r = rows[i]
                subr = sub[ptrb[r]:ptre[r]]
                valr = val[ptrb[r]:ptre[r]]
                valr = self.__sortsparselist(subr, valr)
                subr = sorted(subr)

                j = i + 1
                while j < len(rows):
                    q = rows[j]
                    subq = sub[ptrb[q]:ptre[q]]
                    valq = val[ptrb[q]:ptre[q]]
                    valq = self.__sortsparselist(subq, valq)
                    subq = sorted(subq)
                    
                    # check that the rows are in fact parallell since hash
                    # might not be perfect
                    if not subr == subq:
                        break
                    s = valq[0] / valr[0] # Aq = s*Ar
                    for k in range(len(subr)):
                        # use isclose due to floating number comparison
                        if not isclose(valq[k], s*valr[k]):
                            break

                    bkr, blr, bur = self.getconbound(r)
                    bkq, blq, buq = self.getconbound(q)

                    if s > 0:
                        newlr = max(blq/s, blr)
                        newur = min(buq/s, bur)

                        lfinite = self.__islfinite(bkq) or self.__islfinite(bkr)
                        ufinite = self.__isufinite(bkq) or self.__isufinite(bkr)

                    elif s < 0:
                        newlr = max(buq/s, blr)
                        newur = min(blq/s, bur)

                        lfinite = self.__isufinite(bkq) or self.__islfinite(bkr)
                        ufinite = self.__islfinite(bkq) or self.__isufinite(bkr)

                    self.chgconbound(r, 1, lfinite, newlr)
                    self.chgconbound(r, 0, ufinite, newur)

                    deletedrows.append(q)
                    rows.remove(q)
                    j += 1
                i += 1

        self.removecons(deletedrows)

    def __sortsparselist(self, subi, vali):
        # Adapted from https://stackoverflow.com/a/6618543
        return [val for _, val in sorted(zip(subi, vali))]

    def buildgmcone(self):
        n = self.getnumvar()
        m = self.getnumcon()

        _, blc, buc = self.getconboundslice(0, m)
        _, blx, bux = self.getvarboundslice(0, n)

        # add t as variable
        self.appendvars(1)

        # Build F for cone
        self.appendafes(self.p+1)
        currentrow = 0
        for i in self.I[0]:
            nzi, subi, vali = self.getarow(i)
            for k in range(nzi):
                # -rows in I[0] in F
                self.putafefentry(currentrow, subi[k], -vali[k])
            # u_c in g
            self.putafeg(currentrow, buc[i])
            currentrow += 1
        for i in self.I[1]:
            nzi, subi, vali = self.getarow(i)
            # +rows in I[1]
            self.putafefrow(currentrow, subi, vali)
            # -l_c in g
            self.putafeg(currentrow, -blc[i])
            currentrow += 1
        for j in self.I[2]:
            # -x_j in cone
            self.putafefentry(currentrow, j, -1.0)
            # u_x in g
            self.putafeg(currentrow, bux[j])
            currentrow += 1
        for j in self.I[3]:
            # x_j in cone
            self.putafefentry(currentrow, j, 1.0)
            # -l_x in g
            self.putafeg(currentrow, -blx[j])
            currentrow += 1
        # t in cone as last element. t has position n in vars
        self.putafefentry(self.p, n, 1.0)

        # Set <F, (x,t)> + g in GMcone
        gmdomain = self.appendprimalgeomeanconedomain(self.p+1)
        self.appendacc(gmdomain, range(self.p+1), None)

    def buildexpcone(self):
        n = self.getnumvar()
        m = self.getnumcon()

        _, blc, buc = self.getconboundslice(0, m)
        _, blx, bux = self.getvarboundslice(0, n)

        # Add t's
        self.appendvars(self.p)

        # Build F for cones
        self.appendafes(3*self.p)
        # First build rows corresponding to first element in each exp cone
        currentrow = 0
        for i in self.I[0]:
            nzi, subi, vali = self.getarow(i)
            for k in range(nzi):
                # -rows in I[0] in F
                self.putafefentry(currentrow, subi[k], -vali[k])
            # u_c in g
            self.putafeg(currentrow, buc[i])
            currentrow += 3
        for i in self.I[1]:
            nzi, subi, vali = self.getarow(i)
            # +rows in I[1]
            self.putafefrow(currentrow, subi, vali)
            # -l_c in g
            self.putafeg(currentrow, -blc[i])
            currentrow += 3
        for j in self.I[2]:
            # -x_j in cone
            self.putafefentry(currentrow, j, -1.0)
            # u_x in g
            self.putafeg(currentrow, bux[j])
            currentrow += 3
        for j in self.I[3]:
            # x_j in cone
            self.putafefentry(currentrow, j, 1.0)
            # -l_x in g
            self.putafeg(currentrow, -blx[j])
            currentrow += 3
            
        # 2nd element in each cone is 1
        self.putafeglist(range(1, 3*self.p, 3), [1.0,]*self.p)

        # 3rd element in each cone is t_i. t_i has position n+i in vars
        self.putafefentrylist(range(2, 3*self.p, 3),
                              range(n, n+self.p),
                              [1.0,]*self.p)
           
        expdomain = self.appendprimalexpconedomain()
        for cone in range(self.p):
            self.appendacc(expdomain, range(3*cone, 3*cone+3), None)