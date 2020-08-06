Create your properties
======================

AlphaSMILES can accept new properties. You can find it and code new ones in *'mcts/properties/properties.py'*.

Whenever you add properties or not, the validity of the SMILES will always be checked with this function from the class *Property*:

.. code-block:: python
    :linenos:

    def calculate(self):
        try:
            m = MolFromSmiles("".join(p.config['long_prefix']) + "".join(self.smiles.element))
            self.smiles.properties[p.s_valid] = False
            if m is not None:
                self.smiles.properties["InChI"] = MolToInchi(m)
                m = AddHs(m)
                AllChem.EmbedMolecule(m)
                AllChem.UFFOptimizeMolecule(m)
        except Exception as e:
            print("Error rdkit : " + repr(e))
            m = None
        if m is not None:
            self.smiles.properties[p.s_valid] = True
            with p.lock_update_data:
                p.tree_info[p.info_good] += 1
                self.smiles.properties[p.s_id] = p.tree_info[p.info_good]
        else:
            with p.lock_update_data:
                p.tree_info[p.info_bad] += 1
                self.smiles.properties[p.s_id] = p.tree_info[p.info_bad]
        return m

If the RDKit function *'MolFromSmiles'* success to convert the SMILES into mol then the molecule is valid.
The ID of the molecule depends on the current number of good SMILES or bad SMILES.

The function return the RDKit mol that all decorator of the class Property must return at the end of their *'calculate'* function.

There are already 4 properties coded :

+ *'SAScoreProperty2DDecorator'* : Calculate the Synthetic Accessibility Score.

+ *'CycleProperty2DDecorator'* : Calculate the size of the longest cycle.

+ *'LogPProperty2DDecorator'* : Calculate the LogP of the molecule.

+ *'DFTPropertyDecorator'* : Calculate the wavelength and oscillator strength of a molecule.


Build your new property
-----------------------

To create a new property you have to implement a new class extended from the *PropertyDecorator*.

This is a template of your new property:

.. code-block:: python
    :linenos:

    class YourNewProperty2DDecorator(PropertyDecorator):
    def __init__(self, property_decorated):
        super().__init__(property_decorated)

    def calculate(self):
        m = self.property.calculate()
        if m is not None:
            try:
                # your code to calculate your new property
                # dont forget to manage the exception
                self.smiles.properties["new property"] = calcul_new_property(self.smiles)
            except Exception as e:
                print("Error [your new property] : " + repr(e))
                self.smiles.properties[p.s_valid] = False
                m = None
        return m

Use it in your project
----------------------

After coding the new class, you can add it to your project by adding it to the configuration.

.. code-block:: python

    config['properties'] = [("mcts.properties.properties", "SAScoreProperty2DDecorator"),
                            ("mcts.properties.properties", "CycleProperty2DDecorator"),
                            ("mcts.properties.properties", "LogPProperty2DDecorator"),
                            ("mcts.properties.properties", "DFTPropertyDecorator"),
                            ("mcts.properties.properties", "YourNewProperty2DDecorator"),
                            ]

The first element of the tuple is the packages and module where the class is located, the second is the name of the class. Like :

.. code-block:: python

    from mcts.properties.properties import YourNewProperty2DDecorator
