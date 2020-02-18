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
 
 ### New in Upgrade 1.1.0
 * Automatic Peak Detection
 * Functions to fit peak wavenumbers/associated quantum numbers to functions
 * Updated plotting function to only label peaks with integer values to prevent overlap
 * Centered labels overtop their associated peaks.
 ### New in Bugfix 1.1.1
 * Tested with another data set.
 * Fixed Peak Detection algorithm to account for noise potentially larger than prominance. 
   * This could have ruined filtering algorithm, as <sup>37</sup>Cl peaks would be higher than noise, meaning they would not be filtered. 
 * Identify true initial P and R-Branch peaks. (J = 1, J = 0)
   * These peaks were previously at risk at being filtered out since they can be lower than peaks on either side.
 
 ## Dependencies
 This library is written in `Python 3` and depends on `Numpy`, `Scipy`, and `Matplotlib`. 

## Acknowledgements
Special thanks to my lab partner Forrest, and colleague Erik who provided data for testing.
