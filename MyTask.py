import mosek

class MyTask(mosek.Task):
    def __init__(self):
        super().__init__()
        self.I = [list() for i in range(7)] # index sets
        self.p = int() # parameter for cone

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