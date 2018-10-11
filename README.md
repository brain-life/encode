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
#### Copyright (2017), [Cesar Caiafa](http://web.fi.uba.ar/~ccaiafa), ccaiafa@gmail.com and [Franco Pestilli](https://psych.indiana.edu/directory/faculty/pestilli-franco.html), frakkopesto@gmail.com

## [Stable code release](https://github.com/brain-life/encode/releases/tag/v0.45).

## How to cite the software.
[Caiafa, C. and Pestilli, F. Multidimensional encoding of brain connectomes. Scientific Reports. volume 7, Article number: 11491 (2017)](https://www.nature.com/articles/s41598-017-09250-w).

[Caiafa, C., Sporns, O., Saykin, A., and Pestill, F., Unified representation of tractography and diffusion-weighted MRI data using sparse multidimensional arrays. Advances in Neural Information Processing Systems 30 (NIPS 2017)](http://papers.nips.cc/paper/7021-unified-representation-of-tractography-and-diffusion-weighted-mri-data-using-sparse-multidimensional-arrays)

## Funding.
[![NIH-5UL1TR-001108_05](https://img.shields.io/badge/NIH_5UL1TR-001108_05-green.svg)](https://projectreporter.nih.gov/project_info_details.cfm?aid=9283642&icde=41600065&ddparam=&ddvalue=&ddsub=&cr=4&csb=FY&cs=DESC&pball=)
[![NSF-BCS-1636893](https://img.shields.io/badge/NSF_BCS-1636893-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1636893)

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

### 9. Description of the output data structure (fe.mat)

#### fe fields
```
• fe.name: [text], customized name of the structure.
• fe.type: [text], type description.
• fe.life: [1 × 1 struct], results of the LiFE method, see fe.life fields description below.
• fe.fg: [1 × 1 struct], input fiber group (connectome) information, see fe.fg fields description below.
• fe.roi: [1 × 1 struct], input fiber group (connectome) information, see fe.roi fields description below.
• fe.path: [1 × 1 struct], paths to input data (connectome, dMRI data, etc.).
• fe.rep: [1×1 struct], dMRI data for a repeated measurement if available. 
```

#### fe.life fields
```
• fe.life.M.Phi: [Na × Nv × Nf sptensor], sparse array Φ encoding the connectome.
• fe.life.M.Nphi: Discretization number in azimuth, default = 360.
• fe.life.M.Ntheta: Discretization number in elevation, default = 360.
• fe.life.M.orient: [3 × Na double], matrix containing in its columns the orientations for each Dictionary element.
• fe.life.M.DictSig: [Nθ × Na double], Dictionary matrix containing in its columns the canonical diffusion kernel (demeaned) at different orienta- tions.
• fe.life.xform.img2acpc: [4 × 4 double], Image to ACPC affine trans- form.
• fe.life.xform.acpc2img: [4 × 4 double], ACPC to image affine trans- form.
• fe.life.fibers: not used.
• fe.life.diffusion_signal_img: [Nv × Nθ double], dMRI data.
• fe.life.diffusion_S0_img: [Nv × 10 double], diffusion signal measured at b = 0 (ten times). Usually, we compute S0(v) by averaging over the ten available values.
• fe.life.bvecs: [Nθ × 3 double], each row indicates the gradient 3D di- rection.
• fe.life.bvals: [Nθ × 1 double], each row indicates the used b-value for each gradient.
• fe.life.bvecsindices: [Nθ × 1 int], indices for measurements with non- zero b-value.
• fe.life.imagedim: [Nx , Ny , Nz , Nm ], size of the input 4D dMRI dataset where (Nx,Ny,Nz) are the the sizes in each coordinate, x, y and z, respec- tivelly; and Nm correspond to the number of measurements with b = 0 and b ̸= 0.
• fe.life.voxel2FNpair: not used.
• fe.life.modelTensor: [λ1 , λ2 , λ3 ], parameters of the Diffusion Tensor (DT) model used to generate the Dictionary of diffusion kernels.
• fe.life.fit: this field contains the results of applying the LiFE method, see description below.
```

#### fe.life.fit fields
```
• fe.life.fit.randState: saved by the optimization algorithm.
• fe.life.fit.results: saved by the optimization algorithm.
• fe.life.fit.weights: [Nf × 1 double], fascicles weights obtained as a result of the optimization.
• fe.life.fit.params.fitMethod: [text], name of used optimization method, default is “bbnls”.
```

#### fe.fg fields
```
• fe.fg.name: [text], name of Fiber Group (connectome).
• fe.fg.colorRgb: [R, G, B], color specification for visualization.
• fe.fg.thickness: [double], thickness specificattion for visualization.
• fe.fg.visible: [binary], parameter for visualization.
• fe.fg.seeds: Seeds used by MRTRIX software to generate the connec- tome.
• fe.fg.seedRadius: MRTRIX parameter.
• fe.fg.seedVoxelOffsets: MRTRIX parameter.
• fe.fg.params: additional MRTRIX parameters.
• fe.fg.fibers: {Nf ×1 cell}, set of fascicles (streamlines), generated with a tractography algorithm. Each cell contains the x,y,z coordinates a list of 3D points describing the trajectory of a single fascicle.
```

#### fe.roi fields
```
• fe.roi.name: [text], name of ROI.
• fe.roi.color: [R, G, B], color specification for visualization.
• fe.roi.coords: [Nv ×3 double], voxel coordinates, each row specifies the spatial location of a single voxel.
• fe.roi.visible: [binary], parameter for visualization.
• fe.roi.mesh: not used.
• fe.roi.dirty: parameter for visualization
• fe.roi.query_id: not used.
```

### Notes
This software is not compatible with Matlab 2014b, but only with later versions of MatLab.
