# Encode: Multidimensional encoding of brain connectomes

![alt tag](https://cloud.githubusercontent.com/assets/11638664/18526917/8a73562c-7a90-11e6-93d6-7bd5055b1f32.png)

# About
This software implements a framework to encode structural brain connectomes into multidimensional arrays. These arrays are commonly referred to as [tensors](https://arxiv.org/abs/1403.4462). Encoding Connectomes provides an agile framework for computing over connectome edges and nodes efficiently. We provide several examples of operations that can be performed using the framework.

One major application of the tensor encoding is the implementaion of the [Linear Fascicle Evaluation method](http://francopestilli.github.io/life/), in short [LiFE](http://www.nature.com/nmeth/journal/v11/n10/abs/nmeth.3098.html). The tensor encoding method allows implementing LiFE with dramatic reduction in storage requirements, up to 40x compression factors. Furtheremore, connectome encoding allows performing multiple computational neuroanatomy operations such as tract-dissections, virtual lesions, and connectivity estimates very efficiently using the machine-friendly array operators. 

We provide demos to expain how to:
 (1) Load and encode diffusion-weighted data and tractography models of white matter fascicles, as well as perform multidimensional arrays operations. 
 (2) Build and optimize a Linear Fascicle Evaluation model. 
 (4) Perform neuronatomical segmentations, computational neuroanatomy operations and virtual lesions using the connectome encoding framework.
 (4) Reproduce some fo the figures of article describing the method implemented in thsi toolbox: Caiafa and Pestilli, forthcoming.

## Application.
* Encoding of brain conenctome and associated phenotypes into multidimensional arrays.
* Evaluate the evidence supporting white-matter connectomes generated using [magnetic resonance diffusion-weighted imaging](http://en.wikipedia.org/wiki/Diffusion_MRI) and [computational tractography ](http://en.wikipedia.org/wiki/Tractography).
* Perform statistical inference on white-matter connectomes: Compare white-matter connectomes, show the evidence for white-matter tracts and connections between brain areas.

## License.
#### Copyright (2016), [Franco Pestilli](http://francopestilli.com/), frakkopesto@gmail.com, [Cesar Caiafa](http://web.fi.uba.ar/~ccaiafa), ccaiafa@gmail.com
 
## [Documentation](TBA).

## [Stable code release](TBA).

## How to cite the software.
[Caiafa, C. and Pestilli, F.](Multidimensional encoding of brain connectomes) Multidimensional encoding of brain connectomes (forthcoming.)

## Funding.
This work was supported by grants by the Indiana Clinical and Translational Institute (CTSI, NIH ULTTR001108).

## Installation.
1. Download (https://github.com/brain-life/encode).
2. [Start MatLab](http://www.mathworks.com/help/matlab/startup-and-shutdown.html).
3. Add repository to the [matlab search path](http://www.mathworks.com/help/matlab/ref/addpath.html).

## Dependencies.
* [MatLab](http://www.mathworks.com/products/matlab/).
* [vistasoft](https://github.com/vistalab/vistasoft).
* [Matlab Brain Anatomy (MBA)](https://github.com/francopestilli/mba).

## Getting started.

### 1. [Download the repository](https://github.com/brain-life/encode).
* Download the Encode repository from the TAR/ZIP files linked [here](https://github.com/brain-life/encode/archive/master.zip).
* UNZIP/UNTAR the file.
* Add the encode folder to your matlab search path. To do so in the MatLab prompt type: 
```
   >> addpath(genpath('/my/path/to/the/encode/folder/'))
```

### 2. [Download the vistasoft repository](https://github.com/vistalab/vistasoft).
* Download the VISTASOFT repository from the TAR/ZIP files linked [here](https://github.com/vistalab/vistasoft/archive/master.zip).
* UNZIP/UNTAR the file.
* Add the VISTASOFT folder to your matlab search path. To do so in the MatLab prompt type: 
```
   >> addpath(genpath('/my/path/to/the/VISTASOFT/folder/'))
```
### 3. [Download the MBA repository](https://github.com/francopestilli/mba).
* Download the MBA repository from the TAR/ZIP files linked [here](https://github.com/francopestilli/mba/archive/master.zip).
* UNZIP/UNTAR the file.
* Add the MBA folder to your matlab search path. To do so in the MatLab prompt type: 
```
   >> addpath(genpath('/my/path/to/the/MBA/folder/'))
```

### 4. [Download the Demo Datasets](http://purl.dlib.indiana.edu/iusw/data/2022/20995/Demo_Data_for_Multidimensional_Encoding_of_Brain_Connectomes.tar.gz).
* Download the demo datasets from the repository [doi:10.5967/K8X63JTX](http://purl.dlib.indiana.edu/iusw/data/2022/20995/Demo_Data_for_Multidimensional_Encoding_of_Brain_Connectomes.tar.gz).
* UNTAR the main file Demo_Data_for_Multidimensional_Encoding_of_Brain_Connectomes.tar.gz
* Go inside the folder Demo_Data_for_Multidimensional_Encoding_of_Brain_Connectomes/ and UNZIP the following files: Figs_data.zip, HCP3T.zip, HCP7T, and STN. You can deleted the original .zip files once they are unziped.
* The structures of files and foldes under the main folder should looks like as follows
* feDemoDataPath.m
* Figs_data/
* HCP3T/
* HCP7T/
* README.txt
* STN/
* 
* Add the main data folder (Demo_Data_for_Multidimensional_Encoding_of_Brain_Connectomes/) to your matlab search path. To do so in the MatLab prompt type:
```
   >> addpath(genpath('/my/path/to/the/Demo_Data_for_Multidimensional_Encoding_of_Brain_Connectomes/'))
```
### 5. [Run the demo_connectome_encoding code](/scripts/demos/demo_connectome_encoding.m).
Here you will learn about creating the tensor representation of a connectoms and perform basic operations such as identifying fascicles having a particular spatial orientation in a small voxel area. 
```
  >>  demo_connectome_encoding.m
```
### 6. [Run the demo_connectome_data_comparison code](/scripts/demos/demo_connectome_data_comparison.m).
This code reproduce Fig. 3 of the paper "Multidimensional encoding of brain connectomes", by C. Caiafa and F. Pestilli. 
```
  >>  demo_connectome_data_comparison.m
```
### 7. [Run the demo_virtual_lesion code](/scripts/demos/demo_virtual_lesion.m).
This code allows you to compute virtual lesions on a particular brain dataset and visualize particular major tracts together with their path-neighborhood, i.e. fascicles sharing same voxels. 
```
  >>  demo_virtual_lesion.m
```
### 8. [Run the demo_LiFE code](/scripts/demos/demo_LiFE.m).
This code allows you to compute compute the fascicles weights for two different tractography methods, probabilistic and deterministic tractographies, on a same brain. This is similar to the original LiFE demo in  https://github.com/francopestilli/life but here a full brain dataset is used. The optimization (fitting fascicles weights) runs in about 3 hours on a modern Intel processor with 8GB of RAM. This code has been tested with MatLab 2015b on Ubuntu 15+ and Mac OSX 10.11.
```
  >>  demo_LiFE.m
```

