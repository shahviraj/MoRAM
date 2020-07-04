# MoRAM - Modulo Recovery using Alternating Minimization
### This repository provides matlab implementation for the MoRAM algorithm as discussed in [Signal Reconstruction from Modulo Observations](https://arxiv.org/pdf/1812.00557.pdf).

## Requirements
* The code uses SPGL ([SPGL1: A solver for large-scale sparse reconstruction](https://www.cs.ubc.ca/~mpf/spgl1/index.html)). All the required matlab files are included in the repo (under folder `spgl1-1.9`), however, in case if you experience difficulties using SPGL, kindly install SPGL package on your system using the instructions provided on their website.
* Code is developed and tested with Matlab version 2017a.
* A sample test image of lovett hall (as used in the paper) is provided in the `test_images` folder. Test iamge is loaded, cropped to appropriate size, and transformed using Haar wavelets within the code.

## Instructions
* Add all folders and sub-folders to matlab path.
* The mean-reconstruction-error plots as provided in the paper can be reproduced by running the file `mod_reconst.m`. 
* The image reconstruction results as provided in the paper can be reproduced by running the file `img_mod_reconst.m`.
* Results would be stored in the `./results` directory.
* Following parameters can be tuned in both of the above matlab files:
  - `pr.n` => Length of the input signal
  - `pr.b` => Number of sparse blocks. Needs to be kept at 1. Other values are not supported.
  - `pr.max_iter` => Total number of iterations for Alternating Minimization 
  - `pr.R` => Period of the modulo function 
  - `pr.rho` => Spread of the true measurements, y =A*z
  - `pr.mspan` => List of different values of number of measurements (m) to be considered.
  - `pr.s_span` => List of different sparsity values to be considered.
  - `pr.method` =>  Compressed Recovery method to be used. Use `justice-pursuit` to reproduce the results from the paper.
  - `pr.init_method` => Initial estimation method to be used. Use `simple_rcm` to reproduce the results from the paper. Other methods are experimental, and not to be considered.

## Support
Kindly contact [Viraj Shah](http://virajshah.me) (vjshah3@illinois.edu) in case any help is needed.

## Acknowledgements
Parts of this code is adopted from codes written by Gauri Jagatap. We also acknowledge the developers of SPGL library. 

## Cite
Following bibtex can be used to cite our manuscript:
```
@article{shah2018signal,
  title={Signal reconstruction from modulo observations},
  author={Shah, Viraj and Hegde, Chinmay},
  journal={arXiv preprint arXiv:1812.00557},
  year={2018}
}
```
