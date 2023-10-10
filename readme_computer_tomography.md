# Readme for Ipython Notebook `ComputerTomography.iypnb`

## Overview 

Some experiments with  stuff related more or less to Computer Tomography

## Notebooks

Here only the more *important* Jupyter notebooks (along with supporting material) are summarized. 

### Notebook: `projections_d_phi_eq.ipynb`

The notebook reviews methods for the computation of projections of a 2D images.

Projections are used to compute the `Radon Transform` from which the original images can be recovered (at least approximately) . The inversion procedure is not part of this notebook. The focus is on the following topics:

1) Definition of *coordinate system* often used with the Radon Transform

2) A method for the computation of a projection along a single line which intersects an image is presented

3) A simple test image is provided and used to demonstrate the computational procedure

4) A full sinogram (a set of projections for various tuples (d, theta)) is computed and displayed

5) The *sinogram* of a single point image is computed analytically. The resulting sinogram has a sinusoidal shape; hence the name `sinogram`

6) The sinogram of a synthetic image with 3 points is computed. The result are 3 sinusoids. It it shown that the computed sinogram of a single point closely matches the result from analytical computation (just to prove the validity of the underlying `Python` functions)

7) The computation of sinograms is demonstrated for two *real world* images.

**Note**

A PDF document `projections_d_phi_e1.pdf` gives an overview.

## Literature

The literature I could find is mostly focused on the inversion of the Radon Transform to reconctruct the original image. There is little guidance in how to compute the Radon transform from an image.

The Radon Transform
Carsten Høilund
Aalborg University, VGIS, 07gr721
November 12, 2007

Tomographic Image Reconstruction.
An Introduction.
Milan Zvolský
Lecture on Medical Physics
by Erika Garutti & Florian Grüner
28.11.2014

2012
The Radon Transform and the Mathematics of Medical Imaging
Jen Beatty
Colby College

The Radon transform and its inverse
August 2011
Alain C. France, CEA/Saclay

---

## Sub-directories

### Sub-directory: `images`

test images used in the notebook

`images/testImgRect1.npy` : simple test image

`images/tree_dublin.png` : image of a tree

`images/lamp_dublin.png` : image of lamp

### Sub-directory: `modules`

`modules/intersections.py` : `Python` functions to:

1) compute the intersections of a projection line with the boundaries of an image

2) projection along a single projection line

3) projections along many projection lines with common angle theta





