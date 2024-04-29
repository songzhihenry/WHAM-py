# Weighted Histogram Analysis Method (WHAM) for Replica Molecular Dynamics Simulations

## **Introduction**


Welcome to the repository dedicated to the implementation of the Weighted Histogram Analysis Method (WHAM) tailored explicitly for the analysis of outputs from replica discrete molecular dynamics simulations (rxDMD)<sup>[1](#myfootnote1)</sup>. In rxDMD simulations, a series of neighboring replicas spanning a temperature spectrum is utilized to effectively sample the conformational free energy landscape. 

To overcome local free energy barriers and facilitate diverse conformational state sampling, temperature exchanges between adjacent replicas occur repeatedly, following the Metropolis criterion. The Andersen thermostat regulates the temperatures for each replica, ensuring efficient exploration of phase space<sup>[2](#myfootnote2)</sup>.

Subsequently, WHAM emerges as a powerful statistical method capable of reconstructing the underlying free energy landscape from such biased sampling data obtained in molecular dynamics simulations<sup>[3](#myfootnote3),[4](#myfootnote4)</sup>. Moreover, WHAM provides the ensemble average value of specific parameters as a function of temperature, offering valuable insights into the thermodynamics of the system under study.

This repository serves as a comprehensive toolkit for researchers and practitioners involved in replica molecular dynamics simulations, empowering them with the means to analyze and interpret their simulation data effectively.


## **Dependencies**
Python 3.6 or later versions are recommended. The version of dependencies may vary depending on the Python version.
- [Numpy](https://pypi.org/project/numpy/)
- [Pandas](https://pypi.org/project/pandas/)
- [Argparse](https://pypi.org/project/argparse/)

Make sure you have these dependencies installed before running the scripts.


## **Usage**

### **Download**

1. Clone this repository to your local machine:

```bash
git clone https://github.com/songzhihenry/WHAM-py
```

### **File configuration**

2. Edit the file named *rex.cond* based on the temperature distribution of your replica simulation setup. Note that the temperature values in the *rex.cond* file are in Kelvin.

3. Edit the *flist* file, which stands for "file list". This file should include the file extensions of the trajectory data you intend to analyze, such as *.rg, *.rmsd, etc.*

4. Edit the *fmatrix* file, which stands for "file matrix". This file is useful if you have a large matrix of trajectory data, where the dimensions of the matrix correspond to the number of frames and parameters. This method of data analysis is beneficial for residue-wise or atom-wise interaction contacts that involve extensive parameters.

5. Edit the *wham-ave-task* file to specify and customize WHAM parameters. Ensure to review the file and adjust it according to your simulation case. The file contains comments to guide users in the correct input format.

6. Ensure that the trajectory data is located in the same directory as the script.

### **Running the Script**

```bash
python wham-1.2.py wham-ave-task
```

## **Important Notes**

1. Ensure that the order of your trajectory data starts from 0.
2. The total number of trajectories or replicas should match the number of replicas assigned in the *wham-ave-task* file.


## **Reference**
<a name="myfootnote1">1</a>: Sugita, Y.; Okamoto, Y. *Replica-Exchange Molecular Dynamics Method for Protein Folding.* Chem. Phys. Lett. 1999, 314 (1–2), 141–151.  
<a name="myfootnote2">2</a>: Andersen, H. C. *Molecular Dynamics Simulations at Constant Pressure and/or Temperature.* J. Chem. Phys. 1980, 72 (4), 2384–2393.  
<a name="myfootnote3">3</a>: Kumar, S., Rosenberg, J. M., Bouzida, D., Swendsen, R. H. & Kollman, P. A. *THE weighted histogram analysis method for free‐energy calculations on biomolecules*. I. The method. J Comput Chem 13, 1011–1021 (1992).  
<a name="myfootnote4">4</a>: John D. Chodera, William C. Swope, Jed W. Pitera, Chaok Seok, and Ken A. Dill, *Use of the weighted histogram analysis method for the analysis of simulated and parallel tempering simulations*, Chem. Theory Comput. 2007, 3, 1, 26–41.  
