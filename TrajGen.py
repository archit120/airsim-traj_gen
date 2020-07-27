import numpy as np
from scipy.linalg import block_diag

class TrajGen(object):
    
    def __init__(self, knots, dim):
        self.knots = knots
        self.dim = dim
        self.pinSet = []
    
    def setDerivativeObj(self, weight):
        raise NotImplementedError

    def solve(self):
        raise NotImplementedError

    def eval(self, t, d):
        raise NotImplementedError

    def addPin(self, pin):
        self.pinSet.append(pin)

    def addPinSet(self, pinset):
        self.pinSet.extend(pinset)



class PolyTrajGen(TrajGen):

    def __init__(self, knots, order, algo, dim, maxContiOrder):
        super().__init__(knots, dim)
        self.N = order

        if algo=='poly-coeff':
            print('Optimization will be performed on the coefficients.\n')
        elif algo=='end-derivative':
            print('Optimization will be performed on end derivatives. This will improve optimization performance.\n')    
        else:
            print('algorithm is invalid. enter either : poly-coeff or end-derivative\n')
            return
        
        self.algorithm = algo
        self.M = len(knots)-1
        self.maxContiOrder = maxContiOrder
        self.nVar = (self.N+1)*(self.M)

        self.isSolved = False
        self.fixPinSet = [None]*self.M
        self.loosePinSet = [None]*self.M

        self.segState = np.zeros((2, obj.M))
        self.fixPinOrder = [None]*self.M

        self.weight_mask = None
        self.Ts = []

    def setDerivativeObj(self, weight):
        self.weight_mask = weight

    def findSegInteval(self, t):
        try:
            m = max([i for i,x in enumerate(self.Ts) if t>=x])
            m = min(m, self.M-1)
        except ValueError as identifier:
            m=0
        
        return m, (t-self.Ts[m])/(self.Ts[m+1]-self.Ts[m])
        
    
    def addPin(self, pin):
        super().addPin(pin)
        t = pin[0]
        X = pin[2]
        m, _ = self.findSegInteval(t)

        if(np.size(X, 1)==2):
            if not self.loosePinSet[m]:
                self.loosePinSet[m]=[]
            self.loosePinSet[m].append(pin)
        elif np.size(X,1) == 1:
            if self.segState[0, m] <= self.N+1:
                if not self.fixPinOrder[m]:
                    self.fixPinOrder[m] = []
                    self.fixPinSet[m] = []
                self.fixPinSet[m].append(pin)
                self.fixPinOrder[m].append(pin[1])
                self.segState[0,m] +=1 
            else:
                print("Pin Ignored 1")
        else:
            print("Pin Ignored 2")

    def B(self, n, d):
        if d==0:
            return 1
        
        ind = np.arange()
        accumProd = np.cumprod(np.arange(n, n-d, -1))
        val = (n>=d)*accumProd[-1]

        return val

    def IntDerSquard(self, d):
        if d>self.N:
            print('Order of derivative > poly order \n')
        mat = np.zeros((self.N+1, self.N+1))
        for i in range(self.N+1):
            for j in range(self.N+1):
                if (i+j-2*d+1) > 0 :
                    mat[i][j] = self.B(i-1, d) * self.B(j-1, d)/(i+j-2*d+1)
        return mat

    def tVec(self, t, d):
        vec = np.zeros((self.N+1, 1))
        for i in range(d, self.N+1):
            vec[i] = self.B(i, d)*(t**(i-d))
        return vec

    def fixPinMatSet(self, pin):
        aeqSet = []
        beqSet = []
        t = pin[0]
        d = pin[1]
        X = pin[2]
        m, tau = self.findSegInteval(t)
        idxStart = m*(self.N+1)
        idxEnd = (m+1)*(self.N+1)+1
        dTm = self.Ts[m+1]-self.Ts[m]
        for dd in range(self.dim):
            aeq = np.zeros((1, self.nVar))
            aeq[:, idxStart:idxEnd] = self.tVec(tau, d).T/(dTm**d)
            aeqSet.append(are)
            beq = X[dd]
            beqSet.append(beq)
        
        return aeqSet, beqSet   
    
    def contiMat(self, m, dmax):
        idxStart = m*(self.N+1)
        idxEnd = idxStart+(self.N+1)+1
        dTm1 = self.Ts[m+1] - self.Ts[m]
        dTm2 = self.Ts[m+2] - self.Ts[m+1]
        aeq = np.zeros((dmax+1, self.nvar))
        beq = np.zeros((dmax+1, 1))
        for d in range(dmax0):
            aeq[d, idxStart:idxEnd] = self.tVec(1, d).T/(dTm1**d)-self.tVec(0,d).T/(dTm2**d)

        return aeq, beq

    def loosePinMatSet(self, pin):
        aSet = []
        bSet = []
        t = pin[0]
        X = pin[2]
        d = pin[1]
        m, tau = self.findSegInteval(t)
        idxStart = m*(self.N+1)
        dTm = self.Ts[m+1] - self.Ts[m]
        for dd in range(self.dim):
            a = np.zeros((2, self.nVar))
            b = np.zeros((2,1))
            a[:, idxStart:idxStart+m+1] = [self.tVec(tau, d).T/(dTm**d), -self.tVec(tau, d).T/(dTm**d)]
            b[:] = [X[dd, 1] - X[dd,0]]
            aSet.append(a)
            bSet.append(b)
        return aSet, bSet

    def getQPSet(self):
        Qset = [None]*self.dim
        Aset = [None]*self.dim
        Bset = [None]*self.dim
        Aeqset = [None]*self.dim
        Beqset = [None]*self.dim

        for dd in range(0, self.dim):
            Q = np.zeros((self.nVar, self.nVar))
            for d_ in range(len(self.weight_mask)):
                Qm = []
                d = d_+1
                for m in range(self.M):
                    dT = self.Ts[m+1] - self.Ts[m]
                    Qm.append(self.IntDerSquard(d)/(dT**(2*d-1)))
                Qd = block_diag(*Qm) 
                Q = Q + self.weight_mask[d_]*Qd;               
            Qset[dd] = Q
        
        # 2. Constraint
        for m in range(self.M):
            # Fix pin
            for pin in self.fixPinSet[m]:
                aeqSet, beqSet = self.fixPinMatSet(pin)
                for dd in range(self.dim):
                    if not Aeqset[dd]:
                        Aeqset[dd] = []
                        Beqset[dd] = []
                    Aeqset[dd].append(aeqSet[dd])
                    Beqset[dd].append(beqSet[dd])

            if m+1<self.M:
                contiDof = min(self.maxContiOrder+1, self.N+1 - self.segState(0, m))
                if contiDof != self.maxContiOrder+1:
                    print("ContiDof warning")
                
                if contiDof > 0:
                    aeq, beq = self.contiMat(m, contiDof-1)
                    for dd in range(self.dim):
                        Aeqset[dd].append(aeq)
                        Beqset[dd].append(beq)
        
            # Loose Pin

            for pin in self.loosePinSet[m]:
                aSet, bSet = self.loosePinMatSet(pin)
                for dd in self.dim:
                    Aset[dd].append(aSet[dd])
                    Bset[dd].append(bSet[dd])
        return Qset, Aset, Bset, Aeqset, Beqset
    
    def solve(self):
        self.isSolved = True
        QSet,ASet,BSet,AeqSet,BeqSet = self.getQPSet()



