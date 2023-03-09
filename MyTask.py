import mosek
from helpfunctions import myhash
from math import isclose, floor, ceil
from collections import deque

class MyTask(mosek.Task):
    def __init__(self):
        super().__init__()
        self.I = [list() for i in range(7)] # index sets
        self.p = int() # parameter for cone
        self.neginf = float('-inf')
        self.posinf = float('inf')
        self.feas_tol = 1e-6
        self.size_tol = 1e8
        self.near_lin_dep_tol = 1e5

    def getarownumnzlist(self, start, stop):
        return [self.getarownumnz(i) for i in range(start, stop)]

    def updatep(self):
        self.p = sum(len(self.I[i]) for i in range(4))

    def removeemptyrows(self):
        m = self.getnumcon()
        nz = self.getarownumnzlist(0, m)

        #self.removecons(tuple(i for i in range(m) if nz[i] == 0))
        self.removecons((i for i in range(m) if nz[i] == 0))

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

    def presolve(self):
        self.removeemptyrows()
        self.remove_redundant()
        self.presolve_singleton()
        self.activity_bound()
        self.presolve_domain(1e8)
        self.remove_redundant()
        self.presolve_lindep()

    def presolve_domain(self, maxiters):
        """Performs domain propagation."""
        print('Domain propagation started')
        m = self.getnumcon()
        n = self.getnumvar()

        noboundchanged = True
        self.__calculate_minmaxact()

        # Only cons with 1 or 0 variables which are infinate need to be
        # considered at the start
        addpropcons = tuple(con for con in range(m) 
                if self.minactinfcnt[con] <= 1 or self.maxactinfcnt[con] <= 1)
        propcons = deque(addpropcons, m) #optimised for pop and append
        isinpropcons = [False,] * m #boolean for which are present in propcons
        for con in addpropcons:
            isinpropcons[con] = True

        # A matrix on row and col format
        ptrb_row, ptre_row, sub_row, val_row = self.getarowslice(0, m)
        ptrb_col, ptre_col, sub_col, val_col = self.getacolslice(0, n)

        # constraint bounds are unchanged
        bkc, blc, buc = self.getconboundslice(0, m)

        iters = 0
        while len(propcons) > 0 and iters < maxiters:
            iters += 1
            con = propcons.popleft()
            isinpropcons[con] = False

            subi = sub_row[ptrb_row[con]:ptre_row[con]]
            vali = val_row[ptrb_row[con]:ptre_row[con]]
            for var, aij in zip(subi, vali):
                bkx, blx, bux = self.getvarbound(var) #get bounds for updating
                ubndchg, lbndchg = self.__dom_prop(con, var, bkx, blx, bux, 
                                            bkc[con], blc[con], buc[con], aij)

                if ubndchg or lbndchg:
                    newbkx, newblx, newbux = self.getvarbound(var)
                    subj = sub_col[ptrb_col[var]:ptre_col[var]]
                    valj = val_col[ptrb_col[var]:ptre_col[var]]
                if ubndchg:
                    self.__update_minmaxact_upper(var, bkx, bux, 
                    newbkx, newbux, subj, valj)
                if lbndchg:
                    self.__update_minmaxact_lower(var, bkx, blx, 
                    newbkx, newblx, subj, valj)

                # if bound changed, add all cons with that variable
                if ubndchg or lbndchg:
                    noboundchanged = False

                    # Only cons (i) not present in propcons and (ii) have
                    # a possibility to be bounded in a domain propagation
                    # need to be added
                    addprop = tuple(con for con in subj 
                    if not isinpropcons[con] and (self.minactinfcnt[con] <= 1 
                    or self.maxactinfcnt[con] <= 1))
                    
                    propcons.extend(addprop)
                    for temp in addprop:
                        isinpropcons[temp] = True
        
        print('Domain propagation finished')
        print(iters)
        return noboundchanged

    def activity_bound(self):
        """Detects hidden fixed constraints and
        infeasibility based on activity."""
        m = self.getnumcon()
        bkc, blc, buc = self.getconboundslice(0, m)
        
        self.__calculate_minmaxact()

        for con in range(m):
            if self.maxactinfcnt[con] == 0 and self.__islfinite(bkc[con]):
                if self.maxact[con] <= blc[con] - self.feas_tol:
                    raise Exception('Infeasible')
                if (self.maxact[con] - blc[con] < self.feas_tol
                    and bkc[con] != mosek.boundkey.fx):
                    # if maxact == blc -> fixed constraint
                    # i.e. chg upper to blc
                    self.chgconbound(con, 0, 1, blc[con])
                    # update vectors in order to not call get() which
                    # flushes updates -> slow down
                    bkc[con], buc[con] = mosek.boundkey.fx, blc[con]

            if self.minactinfcnt[con] == 0 and self.__isufinite(bkc[con]):
                if self.minact[con] >= buc[con] + self.feas_tol:
                    raise Exception('Infeasible')
                elif (abs(self.minact[con] - buc[con]) < self.feas_tol 
                      and bkc[con] != mosek.boundkey.fx):
                    # if minact == buc -> fixed constraint
                    # i.e. chg lower to buc
                    self.chgconbound(con, 1, 1, buc[con])
                    # update vectors in order to not call get() which
                    # flushes updates -> slow down
                    bkc[con], blc[con] = mosek.boundkey.fx, buc[con]

    def remove_redundant(self):
        """Removes redundant bounds and constraints based on constraint
        activity."""
        print('Removing redundant constraints started')
        m = self.getnumcon()
        bkc, blc, buc = self.getconboundslice(0, m)

        self.__calculate_minmaxact()

        deletedrows = list()
        for con in range(m):
            # Check redundancy:
            upper_redundant, lower_redundant = False, False
            if self.maxactinfcnt[con] == 0 and self.__isufinite(bkc[con]):
                if self.maxact[con] <= buc[con] + self.feas_tol:
                    upper_redundant = True
                    #self.chgconbound(con, False, False, self.posinf)

            if self.minactinfcnt[con] == 0 and self.__islfinite(bkc[con]):
                if self.minact[con] >= blc[con]  - self.feas_tol:
                    lower_redundant = True
                    #self.chgconbound(con, True, False, self.neginf)

            if bkc[con] == mosek.boundkey.fx:
                if lower_redundant and upper_redundant:
                    deletedrows.append(con)
            else:
                if lower_redundant:
                    self.chgconbound(con, True, False, self.neginf)
                if upper_redundant:
                    self.chgconbound(con, False, False, self.posinf)

        bkc, _, _ = self.getconboundslice(0,m)
        self.removecons([con for con in range(m) if bkc[con] == mosek.boundkey.fr] + deletedrows)

        print(len(deletedrows))
        print('Removing redundant constraints finished')

        # return len(deletedrows) == 0 # True if no removed

    def presolve_singleton(self):
        """Detects singleton cons and turn them into variable cons."""
        print('Removing singletons started')
        m = self.getnumcon()
        n = self.getnumvar()
        singlecons = tuple(con for con in range(m) if self.getarownumnz(con) == 1)
        
        ptrb, ptre, sub, val = self.getarowslice(0, m)
        _, blc, buc = self.getconboundslice(0, m)
        _, blx, bux = self.getvarboundslice(0, n)
        vartype = self.getvartypelist(range(n))

        for con in singlecons:
            var = sub[ptrb[con]:ptre[con]][0]
            a =  val[ptrb[con]:ptre[con]][0]

            if a > 0:
                newu = min(buc[con]/a, bux[var])
                newl = max(blc[con]/a, blx[var])

            else: # aij < 0
                newu = min(blc[con]/a, bux[var])
                newl = max(buc[con]/a, blx[var])

            if vartype[var] == mosek.variabletype.type_int:
                newu = floor(newu + self.feas_tol)
                newl = ceil(newl - self.feas_tol)

            self.chgvarbound(var, False, newu < self.size_tol, newu)
            self.chgvarbound(var, True, newl > -self.size_tol, newl)
        
        self.removecons(singlecons)
        #print(len(singlecons))
        print('Removing singletons finished')

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

        subi, subj, val = self.getatrip()
        bkx, blx, bux = self.getvarboundslice(0,n)

        for i, j, aij in zip(subi, subj, val):
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

    def __update_minmaxact_upper(self, var, oldbkx, oldbux, newbkx, newbux, subj, valj):
        #_, subj, valj = self.getacol(var)
        #newbkx, _, newbux = self.getvarbound(var)
        
        # upper from infinite to finite
        if not self.__isufinite(oldbkx) and self.__isufinite(newbkx):
            du = newbux
            for aij, i in zip(valj, subj):
                if aij > 0.0:
                    self.maxactinfcnt[i] -= 1
                else:
                    self.minactinfcnt[i] -= 1
        else: # tighter upper bound
            du = newbux - oldbux # < 0

        for aij, i in zip(valj, subj):
            if aij > 0.0:
                self.maxact[i] += aij * du
            else:
                self.minact[i] += aij * du

    def __update_minmaxact_lower(self, var, oldbkx, oldblx, newbkx, newblx, subj, valj):
        #_, subj, valj = self.getacol(var)
        #newbkx, newblx, _ = self.getvarbound(var)
        
        if not self.__islfinite(oldbkx) and self.__islfinite(newbkx):
            dl = newblx
            for aij, i in zip(valj, subj):
                if aij > 0.0:
                    self.minactinfcnt[i] -= 1
                else:
                    self.maxactinfcnt[i] -= 1
        else:
            dl = newblx - oldblx # > 0
        
        for aij, i in zip(valj, subj):
            if aij > 0.0:
                self.minact[i] += aij * dl
            else:
                self.maxact[i] += aij * dl

    def __dom_prop(self, con, var, bkx, blx, bux, bkc, blc, buc, aij):
        ubndchg, lbndchg = False, False
        tol = 1e3 * self.feas_tol # minimum increase in bound to be enough

        #aij = self.getaij(con, var)
        #bkc, blc, buc = self.getconbound(con)
        #bkx, blx, bux = self.getvarbound(var)

        if bkx == mosek.boundkey.fx: # already fixed
            return ubndchg, lbndchg

        # calculate inf and sup
        inf, sup = self.getinfsuprow(con, var, aij, bkx, blx, bux)

        # if boundupdate is well-defined
        if self.__isufinite(bkc) and inf != self.neginf:
            newbound = (buc - inf)/aij

            if aij > 0.0: # if coeff is positive update upper
                if self.getvartype(var) == mosek.variabletype.type_int:
                    # add some slack to rounding
                    newbound = floor(newbound + self.feas_tol)

                if (self.__isufinite(bkx) and newbound + tol >= bux 
                    or newbound > self.size_tol): # max size
                    ubndchg = False
                else:
                    #print(var, 'upper tighter', bux - newbound)
                    # max size on bound
                    #self.chgvarbound(var, 0, newbound <= self.size_tol, newbound)
                    self.chgvarbound(var, False, True, newbound)
                    ubndchg = True
                
            else: # if coeff is negative update lower
                if self.getvartype(var) == mosek.variabletype.type_int:
                    newbound = ceil(newbound - self.feas_tol)

                if (self.__islfinite(bkx) and newbound - tol <= blx 
                    or newbound < -self.size_tol):
                    lbndchg = False
                else:
                    #print(var, 'lower tigher', newbound - blx)
                    #self.chgvarbound(var, 1, newbound >= -self.size_tol, newbound)
                    self.chgvarbound(var, True, True, newbound)
                    lbndchg = True

        # if boundupdate is well-defined
        if self.__islfinite(bkc) and sup != self.posinf:
            newbound = (blc - sup)/aij
            if aij > 0.0: # if coeff is positive update lower
                if self.getvartype(var) == mosek.variabletype.type_int:
                    newbound = ceil(newbound - self.feas_tol)

                if (self.__islfinite(bkx) and blx >= newbound - tol
                    or newbound < -self.size_tol):
                    lbndchg = False
                else:
                    #self.chgvarbound(var, 1, newbound >= -self.size_tol, newbound)
                    self.chgvarbound(var, True, True, newbound)
                    lbndchg = True
                
            else: # if coeff is negative update upper
                if self.getvartype(var) == mosek.variabletype.type_int:
                    newbound = floor(newbound + self.feas_tol)

                if (self.__isufinite(bkx) and bux <= newbound + tol 
                    or newbound > self.size_tol):
                    ubndchg = False
                else:
                    #self.chgvarbound(var, 0, newbound <= self.size_tol, newbound)
                    self.chgvarbound(var, False, True, newbound)
                    ubndchg = True

        return ubndchg, lbndchg

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
        """Detects parallell rows in the matrix A and turn them into a single
        con."""
        print('Started removing parallel rows')

        m = self.getnumcon()
        ptrb, ptre, sub, val = self.getarowslice(0, m)
        rowpattern = dict() # dictionary for rows with same non-zero elements
        for i in range(m):
            subi = sub[ptrb[i]:ptre[i]]
            vali = val[ptrb[i]:ptre[i]]
            normvali = [element / vali[0] for element in vali]
            key = myhash(subi, normvali, self.near_lin_dep_tol) # near lin dep tolerance
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
                    if abs(newlr-newur) < self.feas_tol: # set to eq
                        self.chgconbound(r, 1, lfinite, newlr)
                        self.chgconbound(r, 0, lfinite, newlr)
                    elif not newlr <= newur + self.feas_tol:
                        raise Exception('Infeasible')
                    else:
                        self.chgconbound(r, 1, lfinite, newlr)
                        self.chgconbound(r, 0, ufinite, newur)

                    deletedrows.append(q)
                    rows.remove(q)
                    j += 1
                i += 1

        #print(deletedrows)
        self.removecons(deletedrows)
        print(len(deletedrows))
        print('Finished removing parallel rows')
        return len(deletedrows) == 0 # True if not removed

    # def chgconbound(self, con, lower, finite, newbound):
    #     """Adapted method to check feasibility and round for fixed cons."""
        
    #     if not finite:
    #         super().chgconbound(con, lower, finite, newbound)
    #     else:
    #         _, blc, buc = self.getconbound(con)

    #         if lower:
    #             if not newbound <= buc + self.feas_tol:
    #                 raise Exception('Infeasible')

    #             if abs(newbound - buc) < self.feas_tol:
    #                 super().chgconbound(con, lower, finite, buc)
    #             else:
    #                 super().chgconbound(con, lower, finite, newbound)

    #         else: # change upper
    #             if not blc <= newbound + self.feas_tol:
    #                 raise Exception('Infeasible')

    #             if abs(newbound - blc) < self.feas_tol:
    #                 super().chgconbound(con, lower, finite, blc)
    #             else:
    #                 super().chgconbound(con, lower, finite, newbound)

    def chgvarbound(self, var, lower, finite, newbound):
        """Adapted method to round for fixed vars."""

        if not finite:
            super().chgvarbound(var, lower, finite, newbound)
        else:
            _, blx, bux = self.getvarbound(var)
            
            # if new lower bound is close enough to upper -> fixed
            if lower and abs(newbound - bux) < self.feas_tol:
                super().chgvarbound(var, lower, finite, bux)
            # if new upper bound is close enought to lower -> fixed
            elif not lower and abs(newbound - blx) < self.feas_tol:
                super().chgvarbound(var, lower, finite, blx)
            else: # otherwise changed bound to newbound
                super().chgvarbound(var, lower, finite, newbound)


        # if not finite:
        #     super().chgvarbound(var, lower, finite, newbound)
        # else:
        #     _, blx, bux = self.getvarbound(var)

        #     if lower:
        #         assert(newbound <= bux + self.feas_tol, 'Infeasible')

        #         if abs(newbound - bux) < self.feas_tol:
        #             super().chgvarbound(var, lower, finite, bux)
        #         else:
        #             super().chgvarbound(var, lower, finite, newbound)

        #     else:
        #         assert(blx - self.feas_tol <= newbound, 'Infeasible')

        #         if abs(newbound - blx) < self.feas_tol:
        #             super().chgvarbound(var, lower, finite, blx)
        #         else:
        #             super().chgvarbound(var, lower, finite, newbound)

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