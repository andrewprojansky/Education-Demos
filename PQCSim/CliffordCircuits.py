"""
Code to understand simulation of quantum systems through a system of cliffords
and a system of exponentiated free fermion gates.

Project Start Date: 9/21/22
"""

import numpy as np
import random
import time

"""
Defines basis set of clifford gates, and all useful arrays/dictionarys for 
later simulation
"""

H = np.array([[1 / np.sqrt(2), 1 / np.sqrt(2)], [1 / np.sqrt(2), -1 / np.sqrt(2)]])
S = np.array([[1, 0], [0, 1j]])
CNotF = np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
CNotB = np.array([[1,0,0,0],[0,0,0,1],[0,0,1,0],[0,1,0,0]])
Id = np.array([[1, 0], [0, 1]])
zero = np.array([0,0])
rho_0 = np.tensordot(zero, zero,0)
gate_dict = {1: S, 2: H, 3: "CNOT"}
s_dict = {1: "S", 2: "H", 3: "CNOT"}
CNotMult= {("X", "X"): [1, "I"], ("X", "Y"): [1j, "Z"], ("X", "Z"): [-1j, "Y"],
           ("Z", "Z"): [1, "I"], ("Z", "Y"): [-1j, "X"], ("Z", "X"): [1j, "Y"]}

def cc(gate):
    """
    Parameters
    ----------
    gate : array
        input gate to be applied to state

    Returns
    -------
    array
        complex conjuagte transpose of gate

    """

    return np.transpose(np.matrix.conjugate(gate))

def initial_state(dim):
    """
    
    Parameters
    ----------
    dim : int
        2 to the number of qubits

    Returns
    -------
    state : array
        array/matrix of all zeros besides single eigenstate entry

    """

    state = np.zeros((dim,dim))
    state[0][0] = 1
    return state

class Experiment:
    """
    Can either run experiment creating clifford circuits by evolving density 
    matrix via gate application, or uses the stabalizer framework to evolve
    state by evolving the operators
    
    Parameters
    ----------
    dim: int
        number of qubits
    state: list, optional
        Initial state position. The default is the plus Z eigenstate
    num_steps : int, optional
        total number of steps. The default is 1.
    gd: dict, optional
        set of gates for non stabalizer formalism gate simulation
    stabalizer: bool, optional
        boolean as to whether run random gate circuit or random stabalizer circuit
    sgd: dict, optional
        list of gates for stabalizer circuit simulation
    
    """

    def __init__(
        self,
        dim,
        state,
        num_steps: int = 1,
        gd = gate_dict,
        stabalizer = False,
        sgd = s_dict

    ):
        self.num_steps = num_steps
        self.dim = dim
        self.state = state
        self.gd = gd
        self.stabc = stabalizer
        self.stab_dict=[]
        self.sgd = sgd
        self.time_r = 0
        
    ################# non-stabalizer gate functions ##################################
    """
    creates gates of correct size to multiply density operator
    """
    
    def OneQGate(self, gate, qubit):
        """
        Parameters
        ----------
        gate : int (keyword)
            keyword for base gate access from dictionary
        qubit : int
            qubit gate applied to

        Returns
        -------
        fgate : array
            returns gate of correct dimension for Hilbert space

        """

        if qubit == 0:
            fgate = gate
        else:
            fgate = Id
        for i in range(1, self.dim):
            if qubit == i:
                fgate = np.kron(gate, fgate)
            else:
                fgate = np.kron(Id, fgate)
        return fgate

    def CNOT(self, control, target):
        """
        Parameters
        ----------
        control : int
            arbitrary control qubit for CNot application
        target : int
            arbitrary target qubit for CNot application

        Returns
        -------
        fgate : array
            returns cnot gate of correct dimension for Hilbert space


        """
        
        if target > control:
            if target - 1 == 0: 
                fgate = CNotF
            else: 
                fgate = Id
            for i in range(1, self.dim - 1):
                if i == target-1:
                    fgate = np.kron(CNotF, fgate)
                else:
                    fgate = np.kron(Id, fgate)
            for i in range(target - 1, control, -1):
                if i-1 == 0:
                    sgate = CNotF
                else:
                    sgate = Id
                for j in range(1, self.dim-1):
                    if j == i-1:
                        sgate= np.kron(CNotF, sgate)
                    else:
                        sgate = np.kron(Id, sgate)                
                fgate = np.matmul(sgate, fgate)
                fgate = np.matmul(fgate, fgate)
        elif target < control: 
            if control == 1:
                fgate = CNotB
            else:
                fgate = Id
            for i in range(1, self.dim-1):
                if i + 1 == control:
                    fgate = np.kron(CNotB, fgate)
                else:
                    fgate = np.kron(Id, fgate)
            for i in range(target, control-1):
                if i == 0: 
                    sgate = CNotB
                else:
                    sgate = Id
                for j in range(1, self.dim-1):
                    if j == i: 
                        sgate = np.kron(CNotB, sgate)
                    else:
                        sgate = np.kron(Id, sgate)
                fgate = np.matmul(sgate, fgate)
                fgate = np.matmul(fgate, fgate)
        return fgate
    
    ################# stabalizer gate functions ##################################
    
    """
    Gate application in stabalizer formalism; updates stabalizers of state 
    depending on applied gate. Functions serve as look-up tables for update 
    rules
    """
    
    def stab_cnot(self, control, target):
        """
        Parameters
        ----------
        control : int
            arbitrary control qubit for CNot application
        target : int
            arbitrary target qubit for CNot application

        """
        #print((control, target))
        for stab in self.stab_dict:
            if target > control:
                if target in stab.keys():
                    if control not in stab.keys() and stab[target][1] != "X":
                        stab[control] = [1, "Z"]
                    elif control in stab.keys():
                        stab[control][0] = stab[control][0] * CNotMult[("Z", stab[control][1])][0]
                        stab[control][1] = CNotMult[("Z", stab[control][1])][1]
                        if stab[control][1] == "I":
                            del stab[control]
                elif control in stab.keys() and stab[control][1] != "Z":
                    stab[target] = [1, "X"]
            else:
                if control in stab.keys():
                    if target not in stab.keys() and stab[control][1] != "Z":
                        stab[target] = [1, "X"]
                    elif target in stab.keys():
                        stab[target][0] = stab[target][0] * CNotMult[("X", stab[target][1])][0]
                        stab[target][1] = CNotMult[("X", stab[target][1])][1]
                        if stab[target][1] == "I":
                            del stab[target]
                elif target in stab.keys() and stab[target][1] != "X":
                    stab[control] = [1, "Z"]
        
    def stab_oneq(self, gate, qubit):
        """
        Parameters
        ----------
        gate : str (keyword)
            keyword for base gate access from dictionary
        qubit : int
            qubit gate applied to

        """
        
        #print(gate)
        for stab in self.stab_dict:
            if qubit in stab.keys():
                if gate == "H":
                    if stab[qubit][1] == "Z":
                        stab[qubit][1] = "X"
                    elif stab[qubit][1] == "X":
                        stab[qubit][1] = "Z"
                    elif stab[qubit][1] == "Y":
                        stab[qubit][0] = -1 * stab[qubit][0]
                if gate == "S":
                    if stab[qubit][1] == "X":
                        stab[qubit][1] = "Y"
                        stab[qubit][0] = -1 * stab[qubit][0]
                    elif stab[qubit] == "Y":
                        stab[qubit] = "X"
                        stab[qubit][0] = -1 * stab[qubit][0]
                        
    ################# state functions ##################################   

    def final_state(self):
        """
        Based on stabalizer state, list unique state that each (set of )qubit
        is in... do I have to even do this? 
        """
                
    ################# run functions ##################################
    def experiment(self):
        """
        Depending on if a stabalizer circuit or not, applies number of
        randomly chosen gates to state of system

        """
        
        start = time.time()
        if self.stabc== False:
            for i in range(self.num_steps):
                gate = random.randint(1, len(self.gd))
                if self.gd[gate] == "CNOT":
                    control = random.randint(0, self.dim-1)
                    target = random.randint(0, self.dim-1)
                    while control == target: 
                        target = random.randint(0, self.dim-1)
                    Cgate = self.CNOT(control, target)

                else:
                    qubit = random.randint(0, self.dim-1)
                    Cgate = self.OneQGate(self.gd[gate], qubit)
                    
                    
                self.state = np.matmul(Cgate, np.matmul(self.state, cc(Cgate)))
        else:
            for i in range(self.dim):
                self.stab_dict.append({i: [1, "Z"]})
            for j in range(self.num_steps):
                gate = random.randint(1, len(self.gd))
                if self.sgd[gate] == "CNOT":
                    control = random.randint(0, self.dim-1)
                    target = random.randint(0, self.dim-1)
                    while control == target: 
                        target = random.randint(0, self.dim-1)
                    self.stab_cnot(control, target)
                else:
                    qubit = random.randint(0, self.dim - 1)
                    self.stab_oneq(self.sgd[gate], qubit)
                #print(self.stab_dict)
                #print('')
        end = time.time()
        self.time_r = end-start
                    
#%%
dim = 6
#exp = Experiment(dim, state = [initial_state(2**dim)], num_steps=10, stabalizer=False)
exp = Experiment(dim, state = [], num_steps=10, stabalizer=True)
exp.experiment()
#del exp
#%%

