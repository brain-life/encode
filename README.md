# Multidimensional encoding of brain connectomes

![alt tag](https://cloud.githubusercontent.com/assets/11638664/18020005/00d2fb9e-6bad-11e6-861a-d746cca35f0d.png)

# About
This software implements a method to encode both, structural brain connectomes and diffusion-weighted magnetic resonance data, into multidimensional arrays (tensors). The encoding method allows implementing Connectoem Evaluation Models (Pestilli et al., 2015) with dramatic reduction in storage requirements with up to 40x compression factor. 

## Application.
* Encoding of brain conenctome and associated phenotypes into multidimensional arrays.
* Evaluate the evidence supporting white-matter connectomes generated using [magnetic resonance diffusion imaging](http://en.wikipedia.org/wiki/Diffusion_MRI) and [computational tractography ](http://en.wikipedia.org/wiki/Tractography).
* Perform statistical inference on white-matter connectomes: Compare white-matter connectomes, show the evidence for white-matter tracts and connections between brain areas.

## License.
#### Copyright (2016), [Franco Pestilli](http://francopestilli.com/), frakkopesto@gmail.com, [Cesar Caiafa](http://web.fi.uba.ar/~ccaiafa), ccaiafa@gmail.com
 
## [Documentation](TBA).

## [Stable code release](TBA).

## How to cite the software.
[Caiafa, C. and Pestilli, F.](Multidimensional encoding of brain connectomes)

## How to cite the Linear Fascicle Evaluation (LiFE) method.
[Pestilli, Franco, Jason D. Yeatman, Ariel Rokem, Kendrick N. Kay, and Brian A. Wandell. Evaluation and statistical inference for human connectomes. Nature methods 11, no. 10 (2014): 1058-1063.](http://www.nature.com/nmeth/journal/v11/n10/abs/nmeth.3098.html)

## Funding.
This work was supported by grants by the Indiana Clinical and Translational Institute (CTSI, NIH ULTTR001108).

## Installation.
1. Download (https://github.com/brain-life/life).
2. [Start MatLab](http://www.mathworks.com/help/matlab/startup-and-shutdown.html).
3. Add repository to the [matlab search path](http://www.mathworks.com/help/matlab/ref/addpath.html).

## Dependencies.
* [MatLab](http://www.mathworks.com/products/matlab/).
* [vistasoft](https://github.com/vistalab/vistasoft).
* [Matlab Brain Anatomy (MBA)](https://github.com/francopestilli/mba).

## Getting started.
Learn about LiFE by using [file_name.m](github.io/file) in [MatLab](http://www.mathworks.com/help/matlab/startup-and-shutdown.html).

### 1. [Download the repository](https://github.com/brain-life/life).
* Download the LiFE repository from the TAR/ZIP files linked [here](https://github.com/brain-life/life/archive/master.zip).
* UNZIP/UNTAR the file.
* Add the life folder to your matlab search path. To do so in the MatLab prompt type: 
```
   >> addpath(genpath('/my/path/to/the/life/folder/'))
```

### 2. [Download the vistasoft repository](https://github.com/vistalab/vistasoft).
* Download the VISTASOFT repository from the TAR/ZIP files linked [here](https://github.com/vistalab/vistasoft/archive/master.zip).
* UNZIP/UNTAR the file.
* Add the VISTASOFT folder to your matlab search path. To do so in the MatLab prompt type: 
```
   >> addpath(genpath('/my/path/to/the/VISTASOFT/folder/'))
```

### 3. [Download the Demo Data](http://PURL-to-DEMO-data).
* Download the LiFE demo data set from the repository [here](https://demo-data/demo_data.tar.gz).
* UNZIP/UNTAR the file.
* Add the unzipped/untarred Data folder to your matlab search path. To do so in the MatLab prompt type:
```
   >> addpath(genpath('/my/path/to/the/data_demo/folder/'))
```

### 4. [Read the demo_connectome_evalaution documentation](http://URL).
Read the description of the calculations in the documentation inside the file, demo_connectome_evaluation.m by typing the following in the matlab prompt: 
```
  >>  edit demo_connectome_evaluation.m
```

### 5. [Run the demo_connectome_evalaution code](URL).
This final step will run the life_demo code. The code will perform the operations described [here](http://URL). 
```
  >>  demo_connectome_evalaution.m
```
demo_connectome_evalaution.m runs in about 3 hours on a modern Intel processor with 8GB of RAM. This code has been tested with MatLab 2015b on Ubuntu 15+ and Mac OSX 10.11.

