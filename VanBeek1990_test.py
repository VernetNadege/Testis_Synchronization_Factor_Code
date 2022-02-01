import numpy as np
from SynchronizationFactor import *
#the file Synchronization factor is supposed to be in the same directory than this file.


#A small example
b = 68.26/3
a = (100-68.26)/7
#No need to normalize (sum to 1 or whatever) the data
WT1       = [ 0.04 ,0.04, 0.04, 0.02, 0.02, 0.02, 0.04, 0.04, 0.04, 0.04 ]
WT2       = [ 2.86 ,2.86, 2.86, 1.43, 1.43, 1.43, 2.86, 2.86, 2.86, 2.86 ]
# note that WT2 is proportional to WT1.

mutant   = [ a ,a, a, b, b, b, a, a, a, a  ]
#In this example, the stages 4, 5 and 6 contain 68.26 % of the stage tubule cross sections for the mutant, and there
#are the only enriched stages so that the synchronization windows
#starts from the beginning of stage 4 to the end of stage 6.
#This windows contains only (1+1+1)/(1+1+1+2+2+2+2+2+2+2) of the stage tubule cross
#section (for the WT) so that the expected synchronisation factor is:
ExpectedRes = 0.6826/(3/17)

print("processing  mutant and WT1 as WT")
SynchronizationFactor(mutant,WT1)
print("\nprocessing  mutant and WT2 as WT")
SynchronizationFactor(mutant,WT2)
print(f"\nExpected result : {ExpectedRes} \n")

WT3       = [ 1.456 ,2.402, 3.321, 2.123, 2.456, 1.232, 2.666, 2.236, 1.156, 2.452 ]
print("Be carefull, we only consider the two first digits after the coma (for WT) so that a warning message is going to say that WT will be approximated")
print("processing  mutant and WT3 as WT")
res =SynchronizationFactor(mutant,WT3)


#Example of the article of Van Beek and Meistrich.1990
WT      = [11.9,7.6,5.0,5.6,4.9,6.6,20,10.4,4.4,2.1,2.2,7.8,6.4,5.0]
mutant1 = [0,0,0,3.7,1.8,37.2,47.6,8.8,0,0.3,0,0,0,0]
mutant2 = [0.4,2.2,1.1,3.7,2.6,41.3,29.9,10.3,3,3,1.1,1.1,0,0]
mutant3 = [0,0.4,1.4,1.4,2.5,20.7,51.4,18.6,2.5,0,0,0,0,0]
mutant4 = [12.3,7.5,2,2.3,1.4,1.8,9.4,9.9,4.1,2.4,3,15.4,15.8,12.8]
mutant5 = [2.8,2.3,1.5,1.2,2.2,5.5,18.6,25,6.3,4.3,4.6,13.6,9.5,3.2]
print("\nprocessing  PVA36-1")
res1 = SynchronizationFactor(mutant1,WT)
print("\nprocessing  PVA36-2")
res2 = SynchronizationFactor(mutant2,WT)
print("\nprocessing  PVA36-3")
res3 = SynchronizationFactor(mutant3,WT)
print("\nprocessing  PVA128-1")
res4 = SynchronizationFactor(mutant4,WT)
print("\nprocessing  PVA128-2")
res5 = SynchronizationFactor(mutant5,WT)

#we obtain synchronization factor values similar to the one in the article.
#we expect our value to be slightly greater since we test several windows.
#However, for 128-2, the value in the article is 1.77 and ours is 1.75. It may be explained
#by the fact that we use the rounded value from the article.  


