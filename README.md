# extsgn-2d
A finite volume-based second-order 2D solver for the Extended Lagrangian Serre-Green-Naghdi equations on rectangular Cartesian mesh, as described in https://doi.org/10.1016/j.jcp.2022.111901. Available numerical methods are the first-order splitting with the exact ODE part, and the second-order impicit-explicit ARS(2, 2, 2) method. The only available numerical flux function is the HLLC solver.

## Getting started
These instructions will get you a copy of the project up and running on your local machine.

### Prerequisites

I have not rigorously checked the exact minimal versions of the dependencies but everything should work already with rather old versions since there are no complex libraries used in the project. Install these dependencies with your favourite package manager on your system.

For simulations:

* gfortran

For postprocessing:

* python 3
* matplotlib
* plotly
* cv2 (for making output videos)
* Pillow (optional, for reading 2D images as initial data)

### Installing

Fork/clone the repository, install all the requirements on your local machine, and open the project with any IDE or text editor.

### Source file short description

Each file in `scr/`  except `main.f90` is a module.

> `src/`
>
> > `m_aux.f90`  contains procedures for data outputs
> >
> > `main.f90`  main program, only calls procedures implemented in `m_methods.f90` and `m_aux.f90`
> >
> > `m_methods.f90` a major module containing all the initial & boundary conditions, and numerical methods
> >
> > `m_model.f90` contains model-specific functions: pressure, sound speed and time step
> >
> > `m_parameters.f90` the input parameters module, the main interface for a user.
>
> `img/`
>
> > `make-video.py` a script to generate a video from the sequence of images
> >
> > `pic-data.py` a fun script to import 2D images as an initial condition.
> >
> > `plot-sequence.py` for generating a sequence of images from simulation data
> >
> > `plot.py` for plotting 1D and 2D results



## Using the code for your simulations

### Input parameters interface `m_parameters.f90`

* Set the parameters of the model `LAMBDA` and `GG`

* Set final physical time `TFIN` and `CFL` number

* Set mesh parameters

    * Number of cells in x and y directions `NX` and `NY`

    * Domain size `XL`, `XR`, `YL`, `YR`

* Set initial condition flag `SELECTOR_IC`. Available options:
    * 1D Riemann problems in x or y directions: `IC_RP_X`, `IC_RP_Y`
    * Cylindrical and square 2D Riemann problems: `IC_RP_CYL`, `IC_RP_SQR`
    * Generalized Riemann problems: piecewise initial data perturbed with sine functions: `IC_GENRP_SINX`, `IC_GENRP_SINX_SINY`
    * Load initial data from a file: `IC_MATRIX`
* Set initial discontinuity position `XMID` `YMID` and cylunder `RADIUS` for the Riemann problem. `RADIUS` is Also used as the HALF square side for the square 2D RP

* Set boundary condition coefficients of the velocities in the ghost cells:
    * 1.0 corresponds to the transparent BC, and -1.0 corresponds to the Wall BC
* Select numerical method flag `SELECTOR_METHOD`. Available options:
    * First-order splitting `METHOD_FIRST_ORDER_SPLITTING`
    * Second-order IMEX ARS(2, 2, 2) `METHOD_IMEX_ARS_222`
* Choose whether you want to produce intermediate output files by setting `SELECTOR_PERC_OUTPUT` to `PERC_OUTPUT_ON` or to `PERC_OUTPUT_OFF` otherwise. The `N_FILES` parameter controls the number of intermediate output files.



### Runing the simulations

A simple makefile is used:

* `make` to compile the project with the optimization flags for faster computations
* `make run` to run execute the code
* `make clean` to remove all the `.mod` and bin files.
* `make plot` to plot a full 2D contour plot, its horizontal 1D cross-section, and an artsy version of the 2D contour plot



## Postprocessing
### Solution plots

There are three available kinds of plots:
1. Cross-section plot of the 2D data corresponding to the points sampled from the horizontal axis (`y = 0`)
2. Two-dimensional contour plot of the water depth `h` with the numerical schlieren filter `log(1 + log(1 + 25*|grad h|))` applied
3. Artsy canvas version of the water depth `h` 2D contour plot
4. 3D surface plot of the water depth `h`

In order to plot the solution:

* First, make sure that in your terminal, you are in the root directory of the git repository so that you could use the Makefile.

* Then, to draw all the plots mentioned above, use the following command:

    ```shell
    make plot
    ```

  They will appear as `plot_1d_sections.png`, `plot_2d_schlieren.png`, `plot_2d_artsy.png` and `plot_3d_surface.png` in the `img/` directory.

### Intermediate output plots and videos
If you chose to output intermediate simulation data, you can generate images from it:
```zsh
make generate_images
```
which will appear in the `img/` directory, and also, to compile an `.mp4` video out of them:

```zsh
make generate_video
```
which will be saved in the `vid/` folder.

### Bonus: custom image of your choice as initial data for simulations

Imagine possessing the ability to freeze time and utilize a 3D printer to generate an initial configuration of a fluid surface. Upon uploading any picture, the printer converts it to grayscale, with darker pixels indicating higher elevations on the surface. This grayscale image serves as the starting point for a numerical simulation. And then you press the play button and let it gracefully disintegrate over time.

Good news is that you actually have this 3D printer in the `pic-data.py` script. You can choose whatever picture you want, put it to the `img/` directory under the name `input.jpg`, and generate the initial data by executing the command:

```zsh
make generate_matrix
```

Then, set the `SELECTOR_IC` to `IC_MATRIX` in the `m_parameters.f90` module and run the simulation as described in the above sections.


## Guidelines for repository contributors

For the moment, I do not accept the contributions. To use the project for your goals, fork it on your local machine. If you face any problems, feel free to open an issue, and if you want to discuss anything, open a discussion in the Discussions pane.
