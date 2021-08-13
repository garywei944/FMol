## Updates
### 08/12/2021
*I noticed that people are staring this repo recently. It might because deepmind released [alphafold](https://github.com/deepmind/alphafold) as an independent project last month. But this is a project I worked on when I was a sophomore and some parts of the project doesn't really work well. So I decide to resume working on it this summer and next semester. Updates coming soon.*


# FMol
**A simplified drug discovery pipeline** -- generating SMILE molecular with AlphaSMILES, predicting protein structure with AlphaFold, and checking the druggability with fPocket/Amber.

![FMol](doc/images/FMol.png)

## Requirements
- 64-bit Linux, we will use `mamba` for package management so distribution is not a problem.

### Install python environment
***Updates coming soon***


### Fix the issue with Theanos
If the default framework used by keras is Theanos, use the following line to switch to TensorFlow print `Using TensorFlow backend.` / `Using Theanos backend.` when you launch the program:
```
export KERAS_BACKEND='tensorflow'
```

### Third party library
* AlphaSMILES uses 3D calculation(DFT) library [Gaussian 09](https://gaussian.com/) by default. If you want this functionality works well, here are some [guides](https://sites.google.com/site/liuxiaogang0206xg/blog/installinggaussian09inubuntu1404) how to set up Gaussian 09 on Ubuntu.
* I use [RECONSTRUCT](http://www.bioinformatics.org/owl/reconstruct/) to reconstruct protein tertiary structure in `.pdb` format from contact map. This software does not works as expected so far, it's still a beta version and the organization is working on it. It's expected to provide an easy way to reconstruct protein tertiary structure. For chemistry professionals, see [Recovery of protein structure from contact maps](https://reader.elsevier.com/reader/sd/pii/S1359027897000412?token=696157F6372B29B840D7878D5304082B2B16BA96A57049537A997E74FBA51C9E8C4C30BB7CA056C204D930072F126D21). They use [Tinker](https://dasher.wustl.edu/tinker/) to reconstruct the protein tertiary structure.


## Usage
### Quick Start
1. Download AlphaFold weight data from [here](https://storage.googleapis.com/alphafold_casp13_data/alphafold_casp13_weights.zip).
2. Install Gaussian 09 and make sure `g09` works well in your terminal
3. Extract the sample input data in `AlphaSMILES/data_in` provided in `.tar.xz` and `.tar.gz` format.
4. Make a new subfolder `alphafold_pytorch/model` and extract the weight folders into `model`.
5. Modify the variable in `fmol.py` according to your PC.
6. Run `./fmol.py`

### If you only want to use AlphaSMILES or AlphaFold
#### AlphaSMILES
Please check [doc](http://forge.info.univ-angers.fr/~cgrelier/AlphaSMILES/index.html) for usage tutorial. [Cyril-Grl](https://github.com/Cyril-Grl) has made a brilliant documentation for it. I provide some additional input data, sample configurations for `rnn` and `mcts`, and a sample output using the sample configurations. There is also a local version of the documentation if Cyril's website shuts down, it's in `AlphaSMILES/doc/_build/html/index.html`

##### Quick start
If you have Gaussian 09 set up and `g09` works well in your terminal and just want a quick start:
1. Extract the sample input data in `AlphaSMILES/data_in` provided in `.tar.xz` and `.tar.gz` format.
2. Change the options in `AlphaSMILES/main.py`
3. Simply run `AlphaSMILES/main.py`

#### alphafold_pytorch
1. To run the project, you need to firstly download [pre-trained weights](http://bit.ly/alphafold-casp13-weights) from Deepmind repos.
2. Create a folder named `model` under `alpha_fold_pytorch`
3. Extract the weights downloaded in step 1 and move `873731`, `916425`, and `941521` 3 folders into the `model` folder.
4. The samples inputs is provided, so simply run `./alphafold_pytorch/alphafold.sh` to run the project.

##### Remarks
1. Technically we can use original deepmind [AlphaFold](https://github.com/deepmind/deepmind-research/tree/master/alphafold_casp13) rather than alphafold_pytorch. But I got too many error warnings when I run their code and they didn't provide a good way to visualize the output. So I choose [alphafold_pytorch](https://github.com/Urinx/alphafold_pytorch) at last.
2. For more details, check [alphafold_pytorch readme](alphafold_pytorch/README.md)
3. If you encounter issue that says out of GPU memory, uncomment line 16 of `alphafold_pytorch/alphafold.sh`. That allows you to run 3 trainings at a time, not all 8 trainings by default.


### molxalpha
#### utils.py
##### utils.rr_to_cm(*input*)
I provide a method to convert CAPS13-RR file to contact map file that [RECONSTRUCT](http://www.bioinformatics.org/owl/reconstruct/) accepts. It create a contact map file in `.cm` file format within the same folder as the input `.rr` file.
###### Params
* **input**(*string*) - path to the input file
###### Returns
* None

#### fpocket
* Use the `install_fpocket.sh` shell script under `scripts` folder to install fpocket on your machine.
* For more information check their [repo](https://github.com/Discngine/fpocket)

#### Amber
***Updates coming soon***


## What's stucked
* The output file of alphafold comes in `.rr` [casp13-rr format](https://predictioncenter.org/casp13/index.cgi?page=format). It stores the probability of two atoms on the protein chain could contact within 8 angstroms. But fpocket only accept input file in `.pdb` format, which basically stores the 3-D coordinate information of each atom. Reconstructing reliable PDB file from the CASP13-RR file is still an unsolved problem in academic circles. [RECONSTRUCT](http://www.bioinformatics.org/owl/reconstruct/) is a third-party software using [TINKER](https://dasher.wustl.edu/tinker/) package aiming to reconstruct PDB file from `.cm` contact map file format, but does not work well. I wrote a tool to convert CASP13-RR format into contact map format(see `utils.rr_to_cm`).
* Deepmind didn't open-source the procedure of protein tertiary structure prediction, especially the part of training model from CASP PDB dataset. However, it's essential to the accuracy of prediction of arbitrary protein structure. 


## Todo
* Reconstruct the project.
* Using [Tinker](https://dasher.wustl.edu/tinker/) to reconstruct protein tertiary structure is a classical approach.


## Copyright Claim
To make the project easier to deploy on the cloud, I copied and merged some repos into this project according to their licence.
* [AlphaSMILE](https://github.com/Cyril-Grl/AlphaSMILES)
* [alphafold_pytorch](https://github.com/Urinx/alphafold_pytorch)
