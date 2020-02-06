# HCl/DCl IR Ro-Vibrational Spectroscopy Lab
 Code for Pitt's CHEM 1430 HCl/DCl IR Spectroscopy Lab. This lab had a lot of data analysis, but a lot of it was essentially programming the same thing over again for a different transition. So I took the code I had written and condensed it into a few functions that you can use for any transition. I will clean and add the remaining code I have written, which fit the peak wavenumbers vs their transition number to a function, and then returns the fitting parameters, along with computing other interesting/necessary values needed for this lab. Examples of how to use the functions can be found in either the `Examples.py` or the `Examples.ipynb`, which will render in the browser. 
 
 ### Features
 #### Entire Spectrum:
 * Read exported csv spectroscopy file into python
 * Plot entire spectrum
 * Better documentation
 #### For a Single Transition:
 * Automatically detect peaks
 * Label peaks with transition number, wavenumber
 * Label Branches
 * Fit functions to wavenumber vs transition number
 * Plot fitted function, wavenumber vs transition number with title and proper axis labels
 
 ### New in 1.1.0
 * Automatic Peak Detection
 * Functions to fit peak wavenumbers/associated quantum numbers to functions
 * Updated plotting function to only label peaks with integer values to prevent overlap
 * Centered labels overtop their associated peaks.
 
 
 ## Dependencies
 This library is written in `Python 3` and depends on `Numpy`, `Scipy`, and `Matplotlib`. 

## Awknowledgements
Special thanks to my lab partner Forrest.
