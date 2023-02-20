import mosek

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
        return bk == mosek.boundkey.ra or bk == mosek.boundkey.up

    def __islfinite(self, bk):
        return bk == mosek.boundkey.ra or bk == mosek.boundkey.lo 

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
            key = hash(tuple(sorted(sub[ptrb[i]:ptre[i]])))
            if key in rowpattern.keys():
                rowpattern[key].append(i)
            else:
                rowpattern[key] = [i,]

        #for i in range(m):
        #    _, subi, vali = self.getarow(i)
        #    key = hash(tuple(sorted(subi)))
        #    if key in rowpattern.keys():
        #        rowpattern[key].append(i)
        #    else:
        #        rowpattern[key] = [i,]

        parallellrows = list() # pairs of parallell rows
        spairs = list() # normalisation constant for each parallell row
        for rows in rowpattern.values():
            nrows = len(rows)
            if nrows == 1:
                break #only one row

            hashlist = list()
            slist = list()
            for row in rows:
                subi = sub[ptrb[row]:ptre[row]]
                vali = val[ptrb[row]:ptre[row]]

                #_, subi, vali = self.getarow(row)
                vali = self.__sortsparselist(subi, vali)
                s = vali[0]
                slist.append(s)
                vali = [element/s for element in vali]
                hashlist.append(hash(tuple(vali)))

            temp = sorted(zip(hashlist, rows, slist))
            hashlist = [element[0] for element in temp]
            rows = [element[1] for element in temp]
            slist = [element[2] for element in slist]

            i = 0
            while i < nrows-1:
                if hashlist[i] == hashlist[i+1]:
                    parallellrows.append((rows[i], rows[i+1]))
                    spairs.append((slist[i], slist[i+1]))
                    i += 2
                else:
                    i += 1
        
        n = len(parallellrows)
        deletedrows = list()
        for i in range(n):
            pair = parallellrows[i]
            spair = spairs[i]

            bkc, blc, buc = self.getconboundlist(pair)

            # Refers to ZIB: "Presolve Reductions in Mixed Integer Programming"
            # If one of the constraints is free, then discard case
            if mosek.boundkey.fr in bkc:
                break

            # If both are equality constraint, must have same b
            # See case 1 p. 27-28 in article
            elif bkc == [mosek.boundkey.fx, mosek.boundkey.fx]:
                if blc[0] == (spair[0]/spair[1])*blc[1]:
                    deletedrows.append(pair[1])
                else:
                    return 'Infeasible' # possible to raise mosek error?

            # If one of the constraints is fixed
            # see case 2 on p. 28 in article
            elif mosek.boundkey.fx in bkc:
                # 0 or 1 is fixed constraint, other is inequality
                fxIndex = bkc.index(mosek.boundkey.fx)
                ineqIndex = 1-fxIndex

                if self.__isufinite(bkc[ineqIndex]):
                    bq = blc[fxIndex]
                    br = buc[ineqIndex]
                    s = spair[fxIndex] / spair[ineqIndex]
                    
                    if ((bq <= s*br and s > 0)
                        or (bq >= s*br and s < 0)):
                        deletedrows.append(pair[ineqIndex])
                    else:
                        return 'Infeasible'
                
                if self.__islfinite(bkc[ineqIndex]):
                    bq = blc[fxIndex]
                    br = -blc[ineqIndex]
                    s = - spair[fxIndex] / spair[ineqIndex]

                    if ((bq <= s*br and s > 0)
                        or (bq >= s*br and s < 0)):
                        deletedrows.append(pair[ineqIndex]) # this index will be duplicate!
                    else:
                        return 'Infeasible'
            
            # Case 3 on p. 28
            # How to structure this? Might not be able to remove row if only one constraint is removed
            else:
                bkq, bkr = bkc[0], bkc[1] # boundtype for each inequality
                if self.__isufinite(bkq) and self.__isufinite(bkr):
                    bq, br = buc[0], buc[1]
                    s = spair[0] / spair[1]

                if self.__isufinite(bkq) and self.__islfinite(bkr):
                    bq, br = buc[0], -blc[1]
                    s = - spair[0] / spair[1]

                if self.__islfinite(bkq) and self.__isufinite(bkr):
                    bq, br = -blc[0], buc[1]
                    s = - spair[0] / spair[1]

                if self.__islfinite(bkq) and self.__islfinite(bkr):
                    bq, br = -blc[0], -blc[1]
                    s = spair[0] / spair[1]

        self.removecons(deletedrows)

    def __checkineqcase(self,  bq, br, s, r, q):
        if bq <= s*br and s > 0:
            return [r,]
        
        elif bq >= s*br and s > 0:
            return [q,]
        
        elif bq > s*br and s < 0:
            return [r, q]

        elif bq == s*br and s < 0:
            self.appendcons(1)
            m = self.getnumcon()
            self.putconbound(m-1, mosek.boundkey.fx, bq, bq)
            return [r, q]

        elif bq < s*br and s < 0:
            return 'Infeasible'


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