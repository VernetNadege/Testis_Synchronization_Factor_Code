# Repository of reproducible code for calculating Synchronization Factor
The synchronization factor is used as a method for quantitative analysis of the synchrony of spermatogenesis.

Published in Cells:
- Condrea, D.; Souali-Crespo, S.; Féret, B.; Klopfenstein, M.; Faisan, S.; Mark, M.; Ghyselinck, N.B.; Vernet, N. Retinoic Acid Receptor Alpha Is Essential in Postnatal Sertoli Cells but Not in Germ Cells. [Cells 2022, 11, 891](https://doi.org/10.3390/cells11050891)

Funded by grants from:
- CNRS, INSERM, UNISTRA, and Agence Nationale pour la Recherche ([ANR-17-CE14-0010-01](https://anr.fr/Projet-ANR-17-CE14-0010) and [ANR-20-CE14-0022](https://anr.fr/Projet-ANR-20-CE14-0022)).

# Mat and meth
The comparison of the stage frequencies between the wild types (WT) and the mutants was performed using the synchronization factor described by van Beek and Meistrich ([van Beek and Meistrich, 1990](https://github.com/VernetNadege/Testis_Synchronization_Factor_Code/blob/main/VanBeek1990_VitA%20_Synchronisation.pdf)). This factor was computed from a synchronization window. For the mutant, this window was positioned inside the cycle of the seminiferous epithelium so as to contain 68,26% of the stage tubule cross sections (first criterium) and so as to contain the most enriched stages with respect to the wild types (second criterium). The window width for the mutant was then defined as the percentage of the cycle of the seminiferous epithelium (determined in the WT) that elapses during the part of the cycle defined by the positioned window. Finally, the synchronization factor was defined as 0,6826 divided by the window width of the mutant and is expected greater than 1 (since the synchronization window contains the most enriched stages in the mutant). Positioning the window inside the cycle is not a well-defined problem. Whereas the first criterium can easily be checked, the second one is a subjective criterion (especially when all stages are present) so that results are subject to experimenter interpretation bias. That is why, we developed an algorithm that automatically estimates the position of the synchronization window inside the cycle using objective criteria: we computed for the mutant the width of all windows that contain 68,26% of the stage tubule cross sections and retained the one that lead to the highest synchronization factor. 

# How to use Synchronization Factor calculator
Download ["How to use Files" .docx](https://github.com/VernetNadege/Testis_Synchronization_Factor_Code/blob/main/How%20To%20Use%20Files.docx)

# Information about .py files
Following files are supposed to be in the same directory.
-	Condrea2022_RaraData.py
-	VanBeek1990_test.py
-	SynchronizationFactor.py
