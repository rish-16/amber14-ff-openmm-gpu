# OpenMM AMBER14 on GPU 

This repo allows you to run the AMBER14 forcefield on a pre-run MD simulation via OpenMM without having to run _more_ MD. It extracts vector information from the trajectory, stores them in an OpenMM-friendly data format, and the OpenMM simulator parses these coordinates to give frame-wise energies and forces. 

## Installation and Setup

Run the following to download the CUDA-accelerated version of OpenMM. Take note of your CUDA version (you can find out using `nvidia-smi`) and make appropriate changes when downloading `cudatoolkit`.

```bash
conda create -n amber14 python=3.6 -y
conda activate amber14
conda install -c conda-forge cudatoolkit=12.5 -y
conda install -c conda-forge openmm mdanalysis -y
```

Check whether OpenMM can access CUDA:
```bash
python -m openmm.testInstallation
```

You should get benchmarking results that look like this:

```bash
OpenMM Version: 7.6
Git Revision: ad113a0cb37991a2de67a08026cf3b91616bafbe

There are 4 Platforms available:

1 Reference - Successfully computed forces
2 CPU - Successfully computed forces
3 CUDA - Successfully computed forces
1 warning generated.
1 warning generated.
1 warning generated.
1 warning generated.
4 OpenCL - Successfully computed forces

Median difference in forces between platforms:

Reference vs. CPU: 6.29276e-06
Reference vs. CUDA: 6.73166e-06
CPU vs. CUDA: 7.39089e-07
Reference vs. OpenCL: 6.74399e-06
CPU vs. OpenCL: 7.80542e-07
CUDA vs. OpenCL: 2.19129e-07

All differences are within tolerance.
```

## Running OpenMM AMBER14

The complete source code is in `openmm_amber14.py`. They take in a crystal structure (`.pdb`) and a GROMACS trajectory (`.xtc`) for the protein obtained from MD, and compute frame-wise energies and forces. You can save them as `.npy` tensors downstream. If you have a CUDA GPU, OpenMM will automatically detect it and run the forcefield on it.

Relevant demo files for protein `16pk` are included to run the code. It has 400+ residues and the code is scalable over sequence lengths.

> Take note that this code can only extract energies and forces from conformations that are within a trajectory. Pointwise evaluations of new conformations are not possible and may require some form of relaxation (WIP). 

## License
[MIT](https://github.com/rish-16/amber14-ff-openmm-gpu/blob/main/LICENSE)