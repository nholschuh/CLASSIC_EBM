# CLASSIC_EBM
Python library for the classical Energy Balance Model (EBM) (Budyko, Sellers, North, etc.)

## Overview

This repository contains some Python routines for exploring the solutions to the classic EBM of Budyko et. al. The model is a one-dimensional representation of a steady-state single hemisphere climate,

<a href="https://www.codecogs.com/eqnedit.php?latex=-\frac{d}{dx}D(1-x^2)\frac{dT}{dx}&space;&plus;&space;A&space;&plus;&space;BT&space;=&space;QS(x)a(x,x_\mathrm{i})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?-\frac{d}{dx}D(1-x^2)\frac{dT}{dx}&space;&plus;&space;A&space;&plus;&space;BT&space;=&space;QS(x)a(x,x_\mathrm{i})," title="-\frac{d}{dx}D(1-x^2)\frac{dT}{dx} + A + BT = QS(x)a(x,x_\mathrm{i})," /></a>

where the coordinate _x_ is the sine of latitude, _T_ is the zonally-averaged surface temperature, the term _A_ + _BT_ represents outgoing longwave radiation to space and the term _QS_(_x_)_a_(_x_,_x<sub>i</sub>_) represents the solar forcing. The coalbedo _a_ depends on the ice-cap edge _x<sub>i</sub>_ (representing the ice albedo feedback). Heat transport is parameterised as diffusion along the mean meridional temperature gradient with parameter _D_ tuned such that the heat transport profile matches observations.

## Dependencies
  * Python 2.7.14
  * NumPy 1.14.3
  * MatPlotLib 2.2.2
