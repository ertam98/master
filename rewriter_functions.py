import mosek

def getindexset(n, m, bkc, blc, buc, bkx, blx, bux):
    """Calculates index sets.
    n: int - number of variables
    m: int - number of constraints
    bkc: mosek.boundkey[] - type of bound for cons
    blc: float[] - lower bounds for cons
    buc: float[] - upper bounds for cons
    bkx: mosek.boundkey[] - type of bound for vars
    blx: float[] - lower bounds for vars
    bux: float[] - upper bounds for vars

    RETURN: I - list[list[]]
    """

    # Calculate index sets
    I = [list() for i in range(6)]
    for i in range(m):
        boundtype = bkc[i]
        if boundtype == mosek.boundkey.up:
            (I[0]).append(i)
        elif boundtype == mosek.boundkey.lo:
            I[1].append(i)
        elif boundtype == mosek.boundkey.fx:
            I[4].append(i)
    
    for j in range(n):
        boundtype = bkx[j]
        if boundtype == mosek.boundkey.up:
            I[2].append(j)
        elif boundtype == mosek.boundkey.lo:
            I[3].append(j)
        elif boundtype == mosek.boundkey.fx:
            I[5].append(j)

    return I