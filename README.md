# Readme for Ipython Notebook `ComputerTomography.iypnb`

## Overview 

Some experiments with  stuff related more or less to Computer Tomography

## Notebooks

Here only the more *important* Jupyter notebooks (along with supporting material) are summarized. 

### Notebook: `radonTransform_theory_1.ipynb`

Mathematical definition of an image (2D function, real valued)

Definition of a straight line through an image. Based on the definition of a line a projection can be computed.

Mathematical definition of the continous Radon transform.

Mathematical definition of the continous 2D Fourier transform. Definition of periodic repetitions of the 2D image and representation as 2D Fourier series.

Evaluation of the continous image $f(x,\ y)$ at discrete points and introduction of *aliased* Fourier coefficients.

Derivation of the `Central / Fourier Slice Theorem` and its relation to the 1D Fourier transform of a projection for a fixed projection angle.

How to reconstruct the image from its projections / Radon transform. Some possible approaches are mentioned:

1) reconstruction using the `Fourier Slice Theorem` directly.

2) reconstruction using filtered backprojections

3) model based approaches which solve a system of equations iteratively

[radonTransform_theory_1.pdf](radonTransform_theory_1.pdf)


### Notebook: `annotations_article_Mersereau_Oppenheim.ipynb`

While attempting to understand `Digital Reconstruction of Multidimensional Signals from their Projections``

Authors: R.M. Mersereau , A.V. Oppenheim,
Proceeding of the IEEE Vol. 62, No. 10 October 1974

I made some notes ...

**Note**

---

### Notebook: `makeTestImages.ipynb`

The notebook generates some simple test images which could be used to compute projections. Other notebooks will use these test images.

[makeTestImages.pdf](makeTestImages.pdf)

---

### Notebook: `radon_transform_computation_2.ipynb`

Some definitions to explain how to compute projections from an image.

The lines along which projections are computed are defined in the $t, \theta$ domain. The projection line is *discretized* in order to pass through points of the image grid. The projection is computed by summing over all image points which belong to the discretized projection line.

Using *discretization* for the projection lines allows to use `Numpy` for the computation. Thus the method of computation is reasonably fast.

Some projections are provided as a demonstration. For some images the full set of projections for angles in the range $0 \le \theta \lt \pi$ are computed, displayed and stored for later usage.

[radon_transform_computation_2.pdf](radon_transform_computation_2.pdf)

---

### Notebook: `filtered_projections_ex1.ipynb`

A set of projections (aka sinogram) is transformed via DFT, filtered by a `Ram-Lak filter` and transformed back into the signal/projection domain. Using backprojection the original image is reconstructed.

Some test images which have been processed in this way show that the algorithm is capable to restore the images with good accuracy.

(there is still an issue with *low contrast* of the filtered backprojection. More work needed to understand the possible cause and how to remedy it)

[filtered_projections_ex1.pdf](filtered_projections_ex1.pdf)

---

### Notebook: `interpolate_polar_to_cartesian.ipynb`

This notebook is still **experimental** .

From the 'Fourier Slice Theorem' a relationship of the 1D Fourier transform of a projection for some projection angle $\theta$ with the 2D Fourier transform along a polar line is established. To reconstruct the image from the 2D Fourier transform it is necessary to interpolate the 2D transform data from polar coordinates into cartesian coordinates.

The notebooks provides some hints how to do that. 

[interpolate_polar_to_cartesian.pdf](interpolate_polar_to_cartesian.pdf)

---

### Notebook: `learning_2D_interpolation.ipynb`

The `Scipy`numerical library provides methods to interpolate multi-dimensional data.

The notebook works through some examples to get familiar with some methods provides by package `Scipy.interpolate`.

A possible application is the interpolation of data in polar coordates $r, \theta$ into cartesian coordinates. The polar data can be represented on rectangular grid (row: -> angles, cols: radius). The cartesian data points to be interpolated can then be transformed into polar coordinates and use the rectangular grid of polar data as a look-up table.

[learning_2D_interpolation.pdf](learning_2D_interpolation.pdf)

---

### Notebook: `projections_d_phi_e1.ipynb`

**Note**

This notebook is an earlier version of notebook `radon_transform_computation_2.ipynb`. It will no longer be maintained.

The notebook reviews methods for the computation of projections of a 2D images.

Projections are used to compute the `Radon Transform` from which the original images can be recovered (at least approximately) . The inversion procedure is not part of this notebook. The focus is on the following topics:

1) Definition of *coordinate system* often used with the Radon Transform

2) A method for the computation of a projection along a single line which intersects an image is presented

3) A simple test image is provided and used to demonstrate the computational procedure

4) A full sinogram (a set of projections for various tuples (d, theta)) is computed and displayed

5) The *sinogram* of a single point image is computed analytically. The resulting sinogram has a sinusoidal shape; hence the name `sinogram`

6) The sinogram of a synthetic image with 3 points is computed. The result are 3 sinusoids. It it shown that the computed sinogram of a single point closely matches the result from analytical computation (just to prove the validity of the underlying `Python` functions)

7) The computation of sinograms is demonstrated for two *real world* images.

[projections_d_phi_e1.pdf](projections_d_phi_e1.pdf)

---

## Literature

The literature I could find is mostly focused on the inversion of the Radon Transform to reconctruct the original image. There is little guidance in how to compute the Radon transform from an image.

1) The Radon Transform
Carsten Høilund
Aalborg University, VGIS, 07gr721
November 12, 2007

2) Tomographic Image Reconstruction.
An Introduction.
Milan Zvolský
Lecture on Medical Physics
by Erika Garutti & Florian Grüner
28.11.2014

3) 2012
The Radon Transform and the Mathematics of Medical Imaging
Jen Beatty
Colby College

4) The Radon transform and its inverse
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






