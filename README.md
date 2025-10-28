
# Mumford-Type Shimura Curves Contained in the Torelli Locus
This repository accompanies the paper *Mumford-Type Shimura Curves Contained in the Torelli Locus* by Thomas Bouchet, Jeroen Hanselman, Andreas Pieper and Sam Schiavone.  It contains files giving the equations of two families of abelian 4-folds of Mumford type, as well as computations that are part of the proofs. The paper can be found on arXiv: https://arxiv.org/abs/2510.00093.

The Shimura curves were found by essentially interpolating CM points on the families, for which we were able to reconstruct the curves. This repository does **not** contain code to do this reconstruction. The code used for this can be found in:

 - https://github.com/Andreas-Pieper/Genus-4-RM-CM (Applying our reconstruction methods to construct CM curves.)
 - https://github.com/JHanselman/reconstructing-g4 (Reconstructing curves of genus 4 from their Theta constants.)
 - https://github.com/Thittho/Genus-4 (Computing invariants of genus 4 curves and reconstructing genus 4 curves from their invariants.)


### Contents of the repository

Here is a short overview of the files in our repository:

 - **Data**
	 - *latex_table_script.sage* (Script to generate LaTeX tables for the paper. TODO: Change paths in script)
	 - *mumford_jac.py*
- **Future** 
	- hypergeometric.m
	- quatalg.m
- **Paranjape**
	- correspondences237.m
	- correspondences239.m
	- elliptic.m
	- rational.m
- **Semistable**
   	- **237.m** (Proposition 4.5 of [BHPS]) A script to determine the semistable reduction of the 2,3,7 family.
	- **239.m** (Proposition 4.5 of [BHPS]) A script to determine the semistable reduction of the 2,3,9 family.
	- **discriminant.m**: A script to compute the discriminant of the families
- **Shimura**
	- **Shimura_237_param.m** (Table 1 of [BHPS]) Putative CM points of the 237 family.
	- **Shimura_239_param.m** (Table 2 of [BHPS]) Putative CM points of the 239 family.
	- **Shimura_models.m** (Theorem 1.2 of [BHPS]) Contains equations for the families, as well as functions that allow you to create members of the families and check if a specified genus 4 curve lies in the family. (Some of the methods depend on https://github.com/Thittho/Genus-4)
- **old-group-search**
	- group-scratchwork.m
	- reflex-field-script.sage

## Installation instructions

First make sure you have the latest version of Magma installed. (The code was written for Magma V2.28-7).

Then use the command
```
git clone https://github.com/SamSchiavone/Mumford-type.git
```

in the folder you want to download the files to.

For some of the functionality it is necessary to install Bouchet's genus 4 invariants package: https://github.com/Thittho/Genus-4
