TODO
====

Parser
------

Code a parser to separate tokens in SMILES

Links:

http://www.dalkescientific.com/writings/diary/archive/2007/06/25/smiles_states.html

https://docs.eyesopen.com/toolkits/python/oechemtk/SMILES.html

https://open-babel.readthedocs.io/en/latest/FileFormats/SMILES_format.html

https://metamolecular.com/cheminformatics/smiles/

http://opensmiles.org/

https://www.rdocumentation.org/packages/rcdk/versions/3.4.7.1/topics/parse.smiles

https://github.com/ottohahn/SMILES

https://github.com/pyparsing/pyparsing

https://github.com/tsudalab/ChemGE/blob/master/zinc_grammar.py

https://github.com/rdkit/rdkit

http://opensmiles.org/spec/open-smiles-2-grammar.html


asynchronous task and Boinc
---------------------------

Adding the possibility to realise asynchronous task to use Boinc. See in the file *'mcts/properties/dft.py'* to generate the *.inp* file to send to boinc.
After generation, you have to send the *.inp* file to boinc then  adding the SMILES to a data base of "SMILES on boinc" to check if they are ready.

A better way would be hidding the use of Boinc to AlphaSMILES by using an API to calculate the properties of a SMILES.

see the `boinc doc <https://quchempedia.univ-angers.fr/boinc/docs/index.html>`_.

API
---

Building an API to ask the calculations of properties for a SMILES (or RDKit Mol) would be very interesting and reusable.

When using AlphaSMILES or any other machine learning/generator/code that give you SMILES, you can send the SMILES (or RDKit Mol) to the API
with a list of properties you need, if there is DFT calculations, ask for which degree of theory you need. Eventually a kind of score or a function to use for scoring.

The API check in the data base if the properties have already been calculated. If yes, the results are returned to the program (score or properties)
otherwise, the API can decide if the calculation can be launch on Boinc or on a server directly.


