# extsgn-2d
2D code for the Extended Lagrangian Serre-Green-Naghdi equations



## Requirements
Simulations:

* gfortran

Postprocessing:

* python
* gnuplot
* matplotlib



## Setup

Install all the requirements on your machine and open the project with any IDE or text editor.



## Source file short description

Each file in `scr/` is a module. 

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
> > `contourplot.py`
> >
> > `pic-data.py`
> >
> > `plot1d.gnu`
> >
> > `plot2d.gnu`



