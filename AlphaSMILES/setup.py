import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="AlphaSMILES-cg",
    version="0.2.0",
    author="Cyril-Grl",
    description="SMILES generator using MCTS and RNN",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Cyril-Grl/AlphaSMILES",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
)
