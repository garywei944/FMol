# FMol
A simplified drug discovery pipeline -- generating SMILE molecular with AlphaSMILES, predicting protein structure with AlphaFold, and checking the druggability with fPocket/Amber.
![FMol](doc/images/FMol.png)


## Todo
As we can tell from the descrption, there would be 3 parts of the project(molecule, protein, and mol-x-protein). Now I stucked at reconstructing pdb protein tertiary structure from contact map.
* Research on [Tinker](https://dasher.wustl.edu/tinker/) to reconstruct protein tertiary structure.
* Add functions to reconstruct protein CASP-RR files in `fmol.py`.
* Maybe, create visualized configurator to config `rnn` and `mcts` used in AlphaSMILES


## Requirements
In short, AlphaSMILES relies on TensorFlow but alphafold relies on PyTorch, so the better way to run this project is to set up separate virtual environments according to their original documentation and run them separately. However, to save disk space and make things more automatic, here are 2 ways to set up an overall environment.

### Choice 1: Manually install with conda
I do recommend you to use [Anaconda](https://www.anaconda.com/) to manage your packages and environments. A reason is AlphaSMILES uses [rdkit](http://www.rdkit.org/), but they do not provide a way to install via `pip`.   
  
Use the following lines to create a new environment named `fmol` with some packages already installed and then activate it.
```
conda create -n fmol python=3 anaconda
conda activate fmol
```

Then we need to install TensorFlow and PyTorch. Anaconda may stuck at some point solving environments. Don't worry, `Ctrl + C` and try more times the issue would be solved usually.

* If you only want to use CPU for this project
```
conda install tensorflow
conda install pytorch torchvision cpuonly -c pytorch
```

* If you also want to use GPU to shorten the runtime
```
conda install tensorflow-gpu
conda install pytorch torchvision cudatoolkit=10.2 -c pytorch
```
Change the version of `cudatoolkit` correspondingly. More details check [PyTorch](https://pytorch.org/get-started/locally/) website.  
  
Then install other libraries with the following lines.
```
conda install -c rdkit rdkit
conda install -c conda-forge keras
conda install -c omnia cclib
pip install pptree
conda install -c anaconda joblib
conda install -c conda-forge tensorboard
```

Don't forget to install third party libraries if you want the project works as expected. After that, we are done with environments.

### Choice 2: Use spec-file.txt
It's also possible to clone an existing environment from a specification file I provided:
```
conda create -n fmol --file spec-file.txt
```
Then activate it with
```
conda activate fmol
```

### Troubles with Theanos
If the default framework used by keras is Theanos, use the following line to switch to TensorFlow print `Using TensorFlow backend.` / `Using Theanos backend.` when you launch the program:
```
export KERAS_BACKEND='tensorflow'
```
It's configured in PyCharm configure file `.idea/workspace.xml`, but need to set up manually if u don't run this project via PyCharm.

### Third party library
* AlphaSMILES uses 3D calculation(DFT) library [Gaussian 09](https://gaussian.com/) by default. If you want this functionality works well, here are some [guides](https://sites.google.com/site/liuxiaogang0206xg/blog/installinggaussian09inubuntu1404) how to set up Gaussian 09 on Ubuntu.
* I use [RECONSTRUCT](http://www.bioinformatics.org/owl/reconstruct/) to reconstruct protein tertiary structure in `.pdb` format from contact map. This software does not works as expected so far, it's still a beta version and the organization is working on it. It's expected to provide an easy way to reconstruct protein tertiary structure. For chemistry professionals, see [Recovery of protein structure from contact maps](https://reader.elsevier.com/reader/sd/pii/S1359027897000412?token=696157F6372B29B840D7878D5304082B2B16BA96A57049537A997E74FBA51C9E8C4C30BB7CA056C204D930072F126D21). They use [Tinker](https://dasher.wustl.edu/tinker/) to reconstruct the protein tertiary structure.

### Remark
Personally I develop and run this project on an Ubuntu 20.04 instance with CUDA 10.2 + cudnn 7. I didn't test it on mac OS or Windows since my macbook does not have a graphic card and running the bash script on windows via WSL is obviously inefficient. Feel free to open an issue page if you test this project on other platforms but encounter compatibility issues.


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
Please check [doc](http://forge.info.univ-angers.fr/~cgrelier/AlphaSMILES/index.html) for usage tutorial. [Cyril-Grl](https://github.com/Cyril-Grl) has made an brilliant documentation for it. I provide some additional input data, sample configurations for `rnn` and `mcts`, and a sample output using the sample configurations. There is also a local version of the documentation if Cyril's website shuts down, it's in `AlphaSMILES/doc/_build/html/index.html`

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
* I have not include any part of amber in this project. But it's a powerful and useful library in chemistry.


## What's stucked
* The output file of alphafold comes in `.rr` [casp13-rr format](https://predictioncenter.org/casp13/index.cgi?page=format). It stores the probability of two atoms on the protein chain could contact within 8 angstroms. But fpocket only accept input file in `.pdb` format, which basically stores the 3-D coordinate information of each atom. Reconstructing reliable PDB file from the CASP13-RR file is still an unsolved problem in academic circles. [RECONSTRUCT](http://www.bioinformatics.org/owl/reconstruct/) is a third-party software using [TINKER](https://dasher.wustl.edu/tinker/) package aiming to reconstruct PDB file from `.cm` contact map file format, but does not work well. I wrote a tool to convert CASP13-RR format into contact map format(see `utils.rr_to_cm`).
* Deepmind didn't open-source the procedure of protein tertiary structure prediction, especially the part of training model from CASP PDB dataset. However, it's essential to the accuracy of prediction of arbitrary protein structure. 


## Copyright Clain
To make the project easier to deploy on the cloud, I copied and merged some repos into this project according to their licence.
* [AlphaSMILE](https://github.com/Cyril-Grl/AlphaSMILES)
* [alphafold_pytorch](https://github.com/Urinx/alphafold_pytorch)
