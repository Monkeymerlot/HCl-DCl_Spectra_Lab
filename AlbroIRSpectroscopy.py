import numpy as np
import matplotlib.pyplot as plt
from scipy import signal, fftpack, optimize


def readSpectrum(filename):
    """Reads an HCl/DCl ro-vibrational spectrum from a text file. 
    
    Parameters
    ----------
    filename: string
        A string/path of the file containing the spectra.
    
    Returns
    -------
    absorption: array
        array of absorption magnitude data
    wavenumber: array
        array of corresponding data to absorption
    """
    data = np.genfromtxt(filename, skip_header = 2, delimiter = ',').T
    absorption = data[1]
    wavenumber = data[0]
    return absorption, wavenumber


def plotfull(absorption, wavenumber):
    """Plot the full ro-vibrational spectrum of HCl/DCl.
    
    This function is useful to see where the different individual spectra lie, so you can then crop
    and plot them using "detectpeaks".
    
    Parameters
    -------
    absorption: array
        array of absorption magnitude data
    wavenumber: array
        array of corresponding data to absorption
    """
    plt.plot(wavenumber, absorption)
    plt.xlim(np.max(wavenumber), np.min(wavenumber))
    plt.xlabel(" $\\tilde{\\nu}$ $(cm^{-1})$")
    plt.ylabel("Relative Absorption")


def detectpeaks(absorption, wavenumber, height, distance):
    """Detect and label peaks based on their quantum numbers and series for HCl/DCl IR ro-vibrational
    spectroscopy data.
    
    This function takes 2 one-dimensional arrays of the absorption data and the wavenumber data for a 
    single transition spectrum, the minimum height of the peaks, and the minimum distance between two 
    peaks, and then will label the peaks with wavenumber, quantum number, and what branch they belong 
    to (P-Branch or R-Branch). 
    
    Parameters
    ----------
    absorption: one-dimensional array
        Absorption signal for one transition (eg HCl: v = 0 -> 1). Should be the minimum length and 
        correctly trimmed. This function will not work if there is more than 1 spectrum.
    wavenumber: one-dimensional array
        Corresponding wavenumbers to absorption data. Must be the same length as absorption.
    height: number
        The minimal height that an absorption peak must have for detection by this function. Height 
        and distance are the two ways that allow for only detecting signals created by Cl-35 instead
        of signals created by Cl-37.
    distance: number
        The minimal distance between neighboring absorption peaks. See height for more information.
        
    Returns
    -------
    rbranch: list
        List of the wavenumbers where absorption peaks occur for the R-Branch.
    jr: list
        The corresponding quantum numbers for the transitions in rbranch.
    pbranch: list
        List of the wavenumbers where absorption peaks occur for the P-Branch.
    jp: list
        The corresponding quantum numbers for the transitions in pbranch.
        
    Warnings
    --------
    This function must have absorption and wavenumber arrays provided. Otherwise assertion error 
    will be raised.
    Absorption and wavenumber arrays must be the same length.
    
    """
    
    assert len(absorption) == len(wavenumber), ("absorption and wavenumber arrays must be the same length")
    assert len(absorption) != 0, ("No absorption data provided. Must provide an array of absorption data.")
    assert len(wavenumber) != 0, ("No wavenumber data provided. Must provide an array of wavenumber data.")
    
    #calculate the range of absorption values
    absrange = np.ptp(absorption)
    
    #automatic ranges for x, y axis.
    plt.ylim(np.min(absorption) - 0.2*absrange,  np.max(absorption) + 0.2*absrange)
    plt.xlim(np.max(wavenumber), np.min(wavenumber))
    
    plt.xlabel("$\\tilde{\\nu}$ $cm^{-1}$")
    plt.ylabel("Relative Absorption")
    
    #find peaks, using height and distance.
    peaks, _ = signal.find_peaks(absorption, height = height, distance = distance)
    
    #want to separate out the P and R branches, largest spacing is in the middle
    spacing = np.abs(peaks[:-1]-peaks[1:])
    split = np.where(spacing == np.max(spacing))[0][0]
    split = int(0.5*(peaks[split]+peaks[split+1]))
    split = wavenumber[split]
    
    #separate out the P and R branch, add lines and annotate.
    rbranch = []
    pbranch = []
    for peak in peaks:
        if wavenumber[peak] > split:
            rbranch.append(wavenumber[peak])
            plt.annotate(str(wavenumber[peak]), (wavenumber[peak], absorption[peak]),
                horizontalalignment='right')
        else:
            pbranch.append(wavenumber[peak])
            plt.annotate(str(wavenumber[peak]), (wavenumber[peak], absorption[peak]),
                horizontalalignment='left')
        plt.axvline(wavenumber[peak], color = "gainsboro")
        
    #create arrays for quantum numbers
    jr = np.arange(0, len(rbranch), 1)
    jp = np.arange(1, len(pbranch)+1, 1)
    jr = np.flip(jr)
    
    #label each peak with quantum number
    labelheight = np.min(absorption) - 0.05*absrange
    for i in range(len(rbranch)):
        plt.annotate(str(jr[i]), (rbranch[i], labelheight))
    for i in range(len(pbranch)):
        plt.annotate(str(jp[i]), (pbranch[i], labelheight))
    
    #label each peak series
    plt.annotate("R-Branch", (0.5*(rbranch[0] - rbranch[-1])+rbranch[-1], labelheight - 0.075*absrange), 
                horizontalalignment = 'center')
    plt.annotate("P-Branch", (0.5*(pbranch[0] - pbranch[-1])+pbranch[-1], labelheight - 0.075*absrange), 
                horizontalalignment = 'center')
    
    plt.plot(wavenumber, absorption)
    
    return rbranch, jr, pbranch, jp


def maketitle(compound, initial, final):
    """Quick helper function to generate a title for HCl/DCl IR ro-vibrational spectroscopy data.
    
    Takes the compound name, the initial energy level, and the final energy level and generates a 
    properly formatted LaTeX title for a graph.
    
    Parameters
    ----------
    compound: string
        The compound name for the spectrum of the graph. (eg HCl or DCl)
    initial: number
        Initial energy level of the transition.
    final: number
        Final energy level of the transition.
    """
    title = compound + " ($\\nu = " + str(initial) + " \\rightarrow{} " + str(final) + "$)"
    plt.title(title)


def waveindex(value):
    """converts a wavenumber to index. Assuming there are 6000 cm^-1 and 8 points/cm^-1
    """
    index = int((6000-value)*8)
    return index
