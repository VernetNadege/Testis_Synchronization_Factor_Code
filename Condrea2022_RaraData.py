import numpy as np
from SynchronizationFactor import *

#RaraGerm+/+
WT = np.array([
[63, 	88,	47,	61,	71,	50,	77,	68,	112,	73],
[9,  	29,	15,	14,	12,	8,  	8,	12,	20, 	15],
[40, 	43,	28,	40,	41,	37,	37,	35,	40,	34],
[77, 	75,	72,	65,	67,	71,	66,	67,	80,	64],
[21, 	13,	15,	16,	17,	20,	14,	20,	23,	17],
[6, 	13,	15,	11,	11,	5,	11,	11,	10,	10],
[19, 	30,	26,	29,	28,	20,	24,	34,	37,	26],
[16, 	16,	17,	25,	15,	14,	21,	9,	17,	18]]);

Ws = np.sum(WT,1)


print("Synchronization factor for each RaraGerm+/+")
for i in range(WT.shape[1]):
    a = WT[:,i].copy()   # the wild type under analysis is a 
    b = (Ws-WT[:,i]).copy() #we remove the contribution of the WT under analysis to compute the summation of the WT.
    SynchronizationFactor(a,b)

#RaraGerm-/-
Mutant1 = np.array([
[70,	62,	49],
[21,	18,	12],
[37,	33,	36],
[79,	53,	63],
[20,	19,	25],
[11,	7,	10],
[20,	17,	12],
[19,	17,	26]]);

#RaraD/L2
Mutant2 = np.array([
[60,	66,	82],
[24,	23,	30],
[33,	32,	54],
[81,	66,	53],
[16,	15,	29],
[10,	10,	15],
[17,	9,	21],
[39,	18,	26]]);


#RaraD/Germ-
Mutant3 = np.array([
[63,	69,	61],
[14,	16,	9],
[52,	46,	56],
[49,	59,	67],
[24,	22,	13],
[18,	14,	15],
[13,	11,	21],
[17,	22,	29]]);

print(" ")
print("Synchronization factor for RaraGerm-/-")
for i in range(Mutant1.shape[1]):
    a = Ws.copy()
    b = (Mutant1[:,i]).copy()
    SynchronizationFactor(b,a)

print(" ")
print("Synchronization factor for RaraD/L2")
for i in range(Mutant2.shape[1]):
    a = Ws.copy()
    b = (Mutant2[:,i]).copy()
    SynchronizationFactor(b,a)

print(" ")
print("Synchronization factor for RaraD/Germ-")
for i in range(Mutant3.shape[1]):
    a = Ws.copy()
    b = (Mutant3[:,i]).copy()
    SynchronizationFactor(b,a)
