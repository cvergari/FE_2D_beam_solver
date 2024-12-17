# FE beam solver

A simple finite element solver in Matlab for beam problems in 2D AND 3D. It allows for forces as boundary conditions, but also imposed displacements and rotations.


### Getting started
obj = FEM_2D_Solver;

obj = obj.demo();


obj = FEM_3D_Solver;

obj = obj.demo();


### Description
The software is implemented as a class: beam properties can be accessed through the object's properties, while methods are provided to launch the solver and plot results.
 
Two demo methods are provided: demo() and demo2(). They can be used as a starting point to use the software. The results of the demos were validated with ANSYS; APDL scripts to reproduce the results are provided in the folders "Ansys demo".

### TODO
The output of reaction forces is not implemented yet, just because I have not needed it yet.

### Authors
* **Claudio Vergari** - *Initial work* - [cvergari on Github](https://github.com/cvergari/)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

## Version
Version 0.3


## Matlab
[![View FE_beam_solver on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://fr.mathworks.com/matlabcentral/fileexchange/74764-fe_2d_beam_solver)
