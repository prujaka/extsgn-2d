# extsgn-2d
A second-order 2D solver for the Extended Lagrangian Serre-Green-Naghdi equations, as described in [1].

Available numerical methods:

* First-order splitting with the exact ODE part 
* Second-order impicit-explicit ARS(2, 2, 2) method 

The only available numerical flux function is the HLLC solver.

## Requirements
I have not rigorously checked the exact minimal versions of the dependencies but everything should work already with rather old versions since there are no complex libraries used in the project.

Simulations:

* gfortran

Postprocessing:

* python 3
* gnuplot
* matplotlib



## Setup

Install all the requirements on your machine and open the project with any IDE or text editor.



## Source file short description

Each file in `scr/`  except `main.f90` is a module. 

> `src/`
>
> > `aux.f90`  contains procedures for data outputs 
> >
> > `main.f90`  main program, only calls procedures implemented in `methods.f90` and `aux.f90`
> >
> > `methods.f90` a major module containing all the initial & boundary conditions, and numerical methods
> >
> > `model.f90` contains model-specific functions: pressure, sound speed and time step
> >
> > `parameters.f90` the input parameters module, the main interface for a user.
>
> `img/`
>
> > `contourplot.py` numerical schlieren representation of the 2D output
> >
> > `pic-data.py ` a fun script to import 2D images as an initial condition.
> >
> > `plot.py` extracts and plots the horizontal 1D section from the 2D output data



## Input parameters interface `parameters.f90`

* Set the parameters of the model `LAMBDA` and `GG`

* Set final physical time `TFIN` and `CFL` number

* Set mesh parameters

    * Number of cells in x and y directions `NX` and `NY` 

    * Domain size `XL`, `XR`, `YL`, `YR`

* Set initial condition flag `SELECTOR_IC`. Available options:
    * 1D Riemann problems in x or y directions `IC_RP_X`, `IC_RP_Y`
    * Cylindrical and square 2D Riemann problems `IC_RP_CYL`, `IC_RP_SQR`
    * Generalized Riemann problems: piecewise initial data perturbed with sine functions `IC_GENRP_SINX` and  `IC_GENRP_SINX_SINY`
* Set initial discontinuity position `XMID` `YMID` and cylunder `RADIUS` for the Riemann problem. `RADIUS` is Also used as the HALF square side for the square 2D RP

* Set boundary condition coefficients of the velocities in the ghost cells:
    * 1.0 corresponds to the transparent BC, and -1.0 corresponds to the Wall BC
* Select numerical method flag `SELECTOR_METHOD`. Available options:
    * First-order splitting `METHOD_FIRST_ORDER_SPLITTING`
    * Second-order IMEX ARS(2, 2, 2) `METHOD_IMEX_ARS_222`



## Run the simulations

A simple makefile is used. All the modules are compiled together but the project is small so the building time is very short.

* `make compile` to compile the project with the optimization flags for faster computations. For the moment, after the compilation, the `.mod` files will appear directly in the root directory. I will soon fix this behaviour to keep the project tree neat.
* `make run` to run execute the code
* `make debug` for debugging with debug flags
* `make makerun` is a sequence of `make compile` and `make run`
* `make clean` to remove all the `.mod` and bin files. 



## 1D and 2D plots

I have not yet implemented the comfortable authomatic plotting, so some things should be done manually:

* Go to the `postprocessing/` directory:

    ```shell
    cd postprocessing
    ```

* To draw a full 2D contourplot of the water depth `h` with the numerical schlieren filter `log(1 + log(1 + 25*|grad h|))` applied, execute the `contourplot.py` script with python:

    ```shell
    python contourplot.py
    ```

    * You will see the pop-up matplotlib window which will also be saved as `schlieren-2d.png` in the `postprocessing/` directory.

* To draw a 1D cross-section of the 2D data along the x axis:

    * Specify the number of cells you used in the `parameters.f90` module inside the `plot.py` script. I will soon add the authomatic reading funtionality, but for the moment let's do it manually.

        Example: 

        ```python
        n_x = 500
        n_y = 500
        ```

    * Execute the `plot.py` script:

        ```py
        python plot.py
        ```

    * You will see the pop-up matplotlib window which will also be saved as `huvetaw-1d.png` in the `postprocessing/` directory. You can disable the pop-up matplotlib figure by commenting 



Reference:

[1] Sergey Tkachenko, Sergey Gavrilyuk, Jacques Massoni, Extended Lagrangian approach for the numerical study of multidimensional dispersive waves: Applications to the Serre-Green-Naghdi equations, Journal of Computational Physics, Volume 477, 2023, 111901. https://doi.org/10.1016/j.jcp.2022.111901
