#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 11:29:16 2022

@author: andrewprojansky
"""


import numpy as np

"""
Defines important gates for computation; clifford gates, T, and custom gates
"""

H = np.array([[1 / np.sqrt(2), 1 / np.sqrt(2)], [1 / np.sqrt(2), -1 / np.sqrt(2)]])
S = np.array([[1, 0], [0, 1j]])
T = np.array([[1, 0], [0, np.exp(1j*np.pi/4)]])
CNotF = np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
CNotB = np.array([[1,0,0,0],[0,0,0,1],[0,0,1,0],[0,1,0,0]])
Id = np.array([[1, 0], [0, 1]])

def cnot(dim, control, target):
    
    if target > control:
        if target - 1 == 0: 
            fgate = CNotF
        else: 
            fgate = Id
        for i in range(1, dim - 1):
            if i == target-1:
                fgate = np.kron(CNotF, fgate)
            else:
                fgate = np.kron(Id, fgate)
        for i in range(target - 1, control, -1):
            if i-1 == 0:
                sgate = CNotF
            else:
                sgate = Id
            for j in range(1, dim-1):
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
        for i in range(1, dim-1):
            if i + 1 == control:
                fgate = np.kron(CNotB, fgate)
            else:
                fgate = np.kron(Id, fgate)
        for i in range(target, control-1):
            if i == 0: 
                sgate = CNotB
            else:
                sgate = Id
            for j in range(1, dim-1):
                if j == i: 
                    sgate = np.kron(CNotB, sgate)
                else:
                    sgate = np.kron(Id, sgate)
            
            fgate = np.matmul(sgate, fgate)
            fgate = np.matmul(fgate, fgate)
                
                              
    return fgate
            
fgate = cnot(3, 2, 0)
print(fgate)