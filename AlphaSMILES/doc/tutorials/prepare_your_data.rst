Prepare your data
=================

Here we will present you how to configurate the project. Think about installing all the library needed see `requirement <../about.html>`_ .


In order to train the RNN you need data. For this project, we use SMILES to represent the molecules.

Your SMILES file have to be in the `data_in/` repertory. One SMILES per line.

The currents data files are :

+ *all_smiles* : SMILES extracted from the Riken data base:

	+ 3,346,716 SMILES in *all_smiles* file (80 Mo)

	+ 3,042,980 SMILES in *all_smiles_unique* file (without doubles) (73.2 Mo)

+ *emolecules* : SMILES from `eMolecules <https://www.emolecules.com/>`_ :

	+ 22,396,732 SMILES in *emolecules* file (870 Mo)

	+ 22,165,087 SMILES in  *emolecules_clean* (without SMILES with a dot) (860 Mo)

+ *pubchem* : SMILES for PubChemQC

	+ 90,131,772 SMILES in *pubchem* file (4.6 Go)

	+ 76,182,195 SMILES in *pubchem_unique* file (3.8 Go)


You can use the *clean_data* function to delete the "." from the SMILES.

.. autofunction:: data_in.clean_data.clean_data
