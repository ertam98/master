import mosek
from helpfunctions import myhash
from math import isclose, floor, ceil

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

    def removeemptyrows(self):
        m = self.getnumcon()
        nz = self.getarownumnzlist(0, m)

        self.removecons([i for i in range(m) if nz[i] == 0])

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

    def presolve_domain(self, maxiters):
        print('Domain propagation started')
        m = self.getnumcon()
        propcons = set(range(m)) # we don't want duplicate elements
        noboundchanged = True

        self.__calculate_minmaxact()

        iters = 0
        while len(propcons) > 0 and iters < maxiters:
            #print(iters, len(propcons))
            iters += 1
            con = propcons.pop()

            nzi, subi, vali = self.getarow(con)
            bkc, blc, buc = self.getconbound(con) # get bounds
            for k in range(nzi):
                var = subi[k]
                aij = vali[k]
                
                bkx, blx, bux = self.getvarbound(var) # get bounds
                upperboundchanged, lowerboundchanged = \
                    self.__domain_propagation(con, var, aij,
                    bkc, blc, buc, bkx, blx, bux)

                # if bound changed, add all cons with that variable
                if upperboundchanged or lowerboundchanged:
                    noboundchanged = False

                    _, subj, _ = self.getacol(var)
                    propcons.update(subj)

                    nzj, subj, valj = self.getacol(var)
                    newbkx, newblx, newbux = self.getvarbound(var)
                    if not newblx <= newbux:
                        raise Exception('Infeasible')

                if upperboundchanged:
                    self.__update_minmaxact_upper(var, bkx, bux,
                                                    newbkx, newbux,
                                                    nzj, subj, valj)
                if lowerboundchanged:
                    self.__update_minmaxact_lower(var, bkx, blx, 
                                                    newbkx, newblx, 
                                                    nzj, subj, valj)
        
        #print(iters)

        #self.__remove_redundant(1e-6)

        #self.__del_minmax()
        print('Domain propagation finished')
        return noboundchanged

    def remove_redundant(self, tol):
        m = self.getnumcon()
        bkc, blc, buc = self.getconboundslice(0, m)

        self.__calculate_minmaxact()

        deletedrows = list()
        for con in range(m):
            # Check feasibility:
            if self.maxactinfcnt[con] == 0 and self.__islfinite(bkc[con]):
                if self.maxact[con] < blc[con] - tol:
                    raise Exception('Infeasible')

            if self.minactinfcnt[con] == 0 and self.__isufinite(bkc[con]):
                if self.minact[con] > buc[con] + tol:
                    raise Exception('Infeasible')

            # Check redundancy:
            upper_redundant, lower_redundant = False, False

            if self.maxactinfcnt[con] == 0 and self.__isufinite(bkc[con]):
                if self.maxact[con] <= buc[con] + tol:
                    upper_redundant = True

            if self.minactinfcnt[con] == 0 and self.__islfinite(bkc[con]):
                if self.minact[con] >= blc[con]  - tol:
                    lower_redundant = True

            if (bkc[con] in {mosek.boundkey.fx, mosek.boundkey.ra} and
                upper_redundant and lower_redundant):
                deletedrows.append(con)
            elif bkc[con] == mosek.boundkey.up and upper_redundant:
                deletedrows.append(con)
            elif bkc[con] == mosek.boundkey.lo and lower_redundant:
                deletedrows.append(con)

        self.removecons(deletedrows)
        return len(deletedrows) == 0 # True if no removed

    def presolve_singleton(self):
        m = self.getnumcon()
        singlecons = tuple(con for con in range(m) if self.getarownumnz(con) == 1)
        
        # vari = [0,] * len(singlecons)
        # aijs = [0.0,] * len(singlecons)
        # for k in range(len(singlecons)):
        #     _, subi, vali = self.getarow(singlecons[k])
        #     vari[k] = subi[0]
        #     aijs[k] = vali[0]

        for con in singlecons:
            _, subi, vali = self.getarow(con)
            bkc, blc, buc = self.getconbound(con)
            bkx, blx, bux = self.getvarbound(subi[0])

            if vali[0] > 0:
                newu = min(buc/vali[0], bux)
                self.chgvarbound(subi[0], False, newu < 1e8, newu)

                newl = max(blc/vali[0], blx)
                self.chgvarbound(subi[0], True, newl > -1e8, newl)

            else: # aij < 0
                newu = min(blc/vali[0], bux)
                self.chgvarbound(subi[0], False, newu < 1e8, newu)

                newl = max(buc/vali[0], blx)
                self.chgvarbound(subi[0], False, newl > -1e8, newl)


    def __del_minmax(self):
        del self.minact
        del self.maxact
        del self.minactinfcnt
        del self.maxactinfcnt

    def __calculate_minmaxact(self):
        m = self.getnumcon()
        n = self.getnumvar()

        self.minact = [0.0,]*m
        self.maxact = [0.0,]*m
        self.minactinfcnt = [0,]*m
        self.maxactinfcnt = [0,]*m

        ptrb, ptre, sub, val = self.getarowslice(0,m)
        bkx, blx, bux = self.getvarboundslice(0,n)

        for i in range(m):
            nzi = self.getarownumnz(i)
            subi = sub[ptrb[i]:ptre[i]]
            vali = val[ptrb[i]:ptre[i]]
            for index in range(nzi):
                j = subi[index]
                aij = vali[index]

                if self.__isufinite(bkx[j]):
                    if aij > 0.0:
                        self.maxact[i] += aij * bux[j]
                    else:
                        self.minact[i] += aij * bux[j]
                else:
                    if aij > 0.0:
                        self.maxactinfcnt[i] += 1
                    else:
                        self.minactinfcnt[i] += 1

                if self.__islfinite(bkx[j]):
                    if aij > 0.0:
                        self.minact[i] += aij * blx[j]
                    else:
                        self.maxact[i] += aij * blx[j]
                else:
                    if aij > 0.0:
                        self.minactinfcnt[i] += 1
                    else:
                        self.maxactinfcnt[i] += 1

    def __update_minmaxact_upper(self, var, oldbkx, oldbux, newbkx, newbux, nzj, subj, valj):
        #nzj, subj, valj = self.getacol(var)
        #newbkx, _, newbux = self.getvarbound(var)
        
        # upper from infinite to finite
        if not self.__isufinite(oldbkx) and self.__isufinite(newbkx):
            du = newbux
            #print(var, 'new upper', du)
            # since upper bound now is finite, update counter
            for k in range(nzj):
                if valj[k] > 0.0:
                    self.maxactinfcnt[subj[k]] -= 1
                else:
                    self.minactinfcnt[subj[k]] -= 1
        else: # tighter upper bound
            du = newbux - oldbux # < 0
            #print(var, 'tighter upper', du)

        for k in range(nzj):
            if valj[k] > 0.0:
                self.maxact[subj[k]] += valj[k] * du
            else:
                self.minact[subj[k]] += valj[k] * du

    def __update_minmaxact_lower(self, var, oldbkx, oldblx, newbkx, newblx, nzj, subj, valj):
        #nzj, subj, valj = self.getacol(var)
        #newbkx, newblx, _ = self.getvarbound(var)
        
        if not self.__islfinite(oldbkx) and self.__islfinite(newbkx):
            dl = newblx
            #print(var, 'new lower', dl)
            # since lower bound now is finite, update counter
            for k in range(nzj):
                if valj[k] > 0.0:
                    self.minactinfcnt[subj[k]] -= 1
                else:
                    self.maxactinfcnt[subj[k]] -= 1
        else:
            dl = newblx - oldblx # > 0
            #print(var, 'tighter lower', dl)
        
        for k in range(nzj):
            if valj[k] > 0.0:
                self.minact[subj[k]] += valj[k] * dl
            else:
                self.maxact[subj[k]] += valj[k] * dl

    def __domain_propagation(self, con, var, aij, bkc, blc, buc, bkx, blx, bux):
        upperboundchanged, lowerboundchanged = False, False

        if bkx == mosek.boundkey.fx: # already fixed
            return upperboundchanged, lowerboundchanged

        # calculate inf and sup
        inf, sup = self.getinfsuprow(con, var, aij, bkx, blx, bux)

        # if boundupdate is well-defined
        if self.__isufinite(bkc) and inf != self.neginf:
            newbound = (buc - inf)/aij

            if aij > 0.0: # if coeff is positive update upper
                if self.getvartype(var) == mosek.variabletype.type_int:
                    newbound = floor(newbound)

                if self.__isufinite(bkx) and newbound >= bux:
                    upperboundchanged = False
                else:
                    print(var, 'upper tighter', bux - newbound)
                    self.chgvarbound(var, 0, 1, newbound)
                    upperboundchanged = True
                
            else: # if coeff is negative update lower
                if self.getvartype(var) == mosek.variabletype.type_int:
                    newbound = ceil(newbound)

                if self.__islfinite(bkx) and newbound <= blx:
                    lowerboundchanged = False
                else:
                    print(var, 'lower tigher', newbound - blx)
                    self.chgvarbound(var, 1, 1, newbound)
                    lowerboundchanged = True

        # if boundupdate is well-defined
        if self.__islfinite(bkc) and sup != self.posinf:
            newbound = (blc - sup)/aij
            if aij > 0.0: # if coeff is positive update lower
                if self.getvartype(var) == mosek.variabletype.type_int:
                    newbound = ceil(newbound)

                if self.__islfinite(bkx) and blx >= newbound:
                    lowerboundchanged = False
                else:
                    self.chgvarbound(var, 1, 1, newbound)
                    lowerboundchanged = True
                
            else: # if coeff is negative update upper
                if self.getvartype(var) == mosek.variabletype.type_int:
                    newbound = floor(newbound)

                if self.__isufinite(bkx) and bux <= newbound:
                    upperboundchanged = False
                else:
                    self.chgvarbound(var, 0, 1, newbound)
                    upperboundchanged = True

        return upperboundchanged, lowerboundchanged

    def __isufinite(self, bk):
        return bk in {mosek.boundkey.ra, mosek.boundkey.up, mosek.boundkey.fx}

    def __islfinite(self, bk):
        return bk in {mosek.boundkey.ra, mosek.boundkey.lo, mosek.boundkey.fx}

    def getinfsuprow(self, con, var, aij, bkx, blx, bux):
        if self.minactinfcnt[con] > 1 and self.maxactinfcnt[con] > 1:
            return self.neginf, self.posinf

        if self.minactinfcnt[con] > 1:
            inf = self.neginf
        elif self.minactinfcnt[con] == 0:
            if aij > 0.0:
                inf = self.minact[con] - aij * blx
            else: # aij < 0
                inf = self.minact[con] - aij * bux
        else: # == 1
            if aij > 0.0:
                if self.__islfinite(bkx): # meaning another variable is inf
                    inf = self.neginf
                else:
                    inf = self.minact[con]
            else: # aij < 0
                if self.__isufinite(bkx): # meaning another variable is inf
                    inf = self.neginf
                else:
                    inf = self.minact[con]

        if self.maxactinfcnt[con] > 1:
            sup = self.posinf
        elif self.maxactinfcnt[con] == 0:
            if aij > 0.0:
                sup = self.maxact[con] - aij * bux
            else: # aij < 0
                sup = self.maxact[con] - aij * blx
        else: # == 1
            if aij > 0.0:
                if self.__isufinite(bkx): # meaning another variable is inf
                    sup = self.posinf
                else:
                    sup = self.maxact[con]
            else: # aij < 0
                if self.__islfinite(bkx): # meaning another variable is inf
                    sup = self.posinf
                else:
                    sup = self.maxact[con]
        
        return inf, sup

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
        print('Started removing parallel rows')

        m = self.getnumcon()
        ptrb, ptre, sub, val = self.getarowslice(0, m)
        rowpattern = dict() # dictionary for rows with same non-zero elements
        for i in range(m):
            subi = sub[ptrb[i]:ptre[i]]
            vali = val[ptrb[i]:ptre[i]]
            normvali = [element / vali[0] for element in vali]
            key = myhash(subi, normvali, 1e-5) # near lin dep tolerance
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

                    if s > 0.0:
                        newlr = max(blq/s, blr)
                        newur = min(buq/s, bur)

                        lfinite = self.__islfinite(bkq) or self.__islfinite(bkr)
                        ufinite = self.__isufinite(bkq) or self.__isufinite(bkr)

                    else:
                        newlr = max(buq/s, blr)
                        newur = min(blq/s, bur)

                        lfinite = self.__isufinite(bkq) or self.__islfinite(bkr)
                        ufinite = self.__islfinite(bkq) or self.__isufinite(bkr)

                    # Check feasibility
                    if not newlr <= newur:
                        raise Exception('Infeasible')

                    self.chgconbound(r, 1, lfinite, newlr)
                    self.chgconbound(r, 0, ufinite, newur)

                    deletedrows.append(q)
                    rows.remove(q)
                    j += 1
                i += 1

        #print(deletedrows)
        self.removecons(deletedrows)
        print('Finished removing parallel rows')
        return len(deletedrows) == 0 # True if not removed

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