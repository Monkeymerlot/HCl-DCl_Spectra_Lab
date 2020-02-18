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
    
    This function is useful to see where the different individual 
    spectra lie, so you can then crop and plot them using "detectpeaks".
    
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


def detectpeaks(absorption, wavenumber, filter_peaks = True, height = None, distance = 40):
    """Detect and label peaks based on their quantum numbers and series
    for HCl/DCl IR ro-vibrational spectroscopy data.
    
    This function takes 2 one-dimensional arrays of the absorption data
    and the wavenumber data for a single transition spectrum, the 
    minimum height of the peaks, and the minimum distance between two 
    peaks, and then will label the peaks with wavenumber, quantum 
    number, and what branch they belong to (P-Branch or R-Branch). 
    
    Parameters
    ----------
    absorption: one-dimensional array
        Absorption signal for one transition (eg HCl: v = 0 -> 1). 
        Should be the minimum length and correctly trimmed. This 
        function will not work if there is more than 1 spectrum.
    wavenumber: one-dimensional array
        Corresponding wavenumbers to absorption data. Must be the same
        length as absorption.
    filter_peaks: Boolean
        Filter peaks will detect all peaks, both HCl-35 and HCl-37, and
        then filter out peaks which are smaller than the peaks on either
        side of it. Peaks on the edges (largest j values in either the
        P-Branch or R-Branch respectively) are not filtered nor checked.
    height: number
        The minimal height that an absorption peak must have for 
        detection by this function. Height and distance are the two ways
        that allow for only detecting signals created by Cl-35 instead
        of signals created by Cl-37. Must be provided if 
        filter_peaks = False, otherwise will be ignored and automatic 
        filtering used.
    distance: number
        The minimal distance between neighboring absorption peaks. See 
        height for more information.
        
    Returns
    -------
    rbranch: list
        List of the wavenumbers where absorption peaks occur for the 
        R-Branch.
    jr: list
        The corresponding quantum numbers for the transitions in 
        rbranch.
    pbranch: list
        List of the wavenumbers where absorption peaks occur for the 
        P-Branch.
    jp: list
        The corresponding quantum numbers for the transitions in 
        pbranch.
        
    Warnings
    --------
    This function must have absorption and wavenumber arrays provided.
        Otherwise assertion error will be raised.
    Absorption and wavenumber arrays must be the same length.
    
    """
    
    assert len(absorption) == len(wavenumber), ("absorption and wavenumber arrays must be the same length")
    assert len(absorption) != 0, ("No absorption data provided. Must provide an array of absorption data.")
    assert len(wavenumber) != 0, ("No wavenumber data provided. Must provide an array of wavenumber data.")
    
    #calculate the range of absorption values
    absrange = np.ptp(absorption)
    
    #automatic ranges for x, y axis.
    plt.ylim(np.min(absorption) - 0.2*absrange,  
             np.max(absorption) + 0.2*absrange)
    plt.xlim(np.max(wavenumber), np.min(wavenumber))
    
    plt.xlabel("$\\tilde{\\nu}$ $cm^{-1}$")
    plt.ylabel("Relative Absorption")

    
    if filter_peaks == False and height != None:
        #find peaks, using height and distance.
        peaks, _ = signal.find_peaks(absorption, height = height, distance = distance)
    else:
        if filter_peaks == False and height == None:
            print("Missing Filter Parameter: height not specified.")
            print("Using Automatic Filtering method instead.")
        #find all peaks (Cl35 and Cl37) and then filter out Cl37 peaks
        peaks, props = signal.find_peaks(absorption, height = 1.005*np.mean(absorption), prominence = 0.015)
        peak_prom = props['prominences']
        
        #set initial variables for filtering
        filt_peaks = []
        
        #calculate P and R Branch split
        spacing = np.abs(peaks[:-1]-peaks[1:])
        split = np.where(spacing == np.max(spacing))[0][0]
        #assign which peak is the first P R-Branch peak
        p_start = peaks[split+1]
        r_start = peaks[split]
        
        #The first Peaks for P and R should be similar. If not, then we have 
        #used the Cl-37 Peak for r_start, and want to use 1 to the left instead
        if abs(peak_prom[split+1] - peak_prom[split])/(peak_prom[split+1]) > 0.20:
            r_start = peaks[split-1]
        
        for i in range(1, len(peaks)-1):
            
            peak = peaks[i]
            
            prev_peak = peaks[i-1]
            
            next_peak = peaks[i+1]
            
            #REDO using prominances instead of peak heights.
            #calculate percentage of peak absorption for next peak and previous peak
            prev_perc = abs(absorption[peak] - absorption[prev_peak])/(absorption[prev_peak])
            next_perc = abs(absorption[peak] - absorption[next_peak])/(absorption[next_peak])
            
            if next_perc > 0.9:
                next_peak = peaks[i+2]
                next_perc = abs(absorption[peak] - absorption[next_peak])/(absorption[next_peak])

            #Filtering peaks.
            if absorption[peak] < absorption[prev_peak] and absorption[peak] < absorption[next_peak] and peak != p_start and peak != r_start:
                #edge case correction for lowest j peaks in R and P
                filt_peaks.append(i)

        peaks = np.delete(peaks, filt_peaks)

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
            plt.annotate(str(int(wavenumber[peak])), (wavenumber[peak], absorption[peak]),
                        horizontalalignment='center')
        else:
            pbranch.append(wavenumber[peak])
            plt.annotate(str(int(wavenumber[peak])), (wavenumber[peak], absorption[peak]),
                        horizontalalignment='center')
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
    """Quick helper function to generate a title for HCl/DCl IR 
    ro-vibrational spectroscopy data.
    
    Takes the compound name, the initial energy level, and the final
    energy level and generates a properly formatted LaTeX title for a
    graph.
    
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
    """converts a wavenumber to index. Assuming there are 6000 cm^-1 and
    8 points/cm^-1.
    """
    index = int((6000-value)*8)
    return index


def fit_p_branch(pbranch, jp, fundamental = False):
    """Fits a specified function to P-Branch transition wavenumbers to
    their corresponding quantum numbers. Also will plot function and
    raw data.
    
    The specified function is 
    nu = c -(2*b - 3*a) * j - 2*a * j**2 + 4*d * j**3 for the overtone 
    transition, where c = 2 * ve - 6 * ve * xe, or 
    nu = c - (2*b - 2*a) * j - a * j**2 + 4*d * j**3 for the fundamental
    transition, where c = ve - 2 * ve * xe. 
    
    Parameters
    ----------
    pbranch: list
        List of the wavenumbers where absorption peaks occur for the 
        P-Branch.
    jp: list
        The corresponding quantum numbers for the transitions in 
        pbranch.
    fundamental: Boolean
        Used to determine fitting function to use. True will use the
        function for the fundamental transitions, while False will use
        the function for the first overtone transitions.
        
    Returns
    -------
    a: Float
        Rotational-vibrational coupling constant. Usually denoted a_e.
    b: Float
        Rotational constant. Usually denoted B_e
    c: Float
        Intercept of the function. For the fundemental transition this 
        value is equal to c = ve - 2 * vexe, where ve is the harmonic
        frequency and vexe is the anharmonicity constant.
    d: Float
        Centrifugal Distortion Constant, not to be confused with D_o, 
        which is the Bond Energy.
    pcov: list
        List of corresponding errors at one standard deviation. Found by
        taking the square root of the diagonal of the covariant matrix.
    """
    #choose function to fit to depending on transition
    if fundamental == False:
        def f(j, c, b, a, d):
            nu = c -(2*b - 3*a) * j - 2*a * j**2 + 4*d * j**3
            return nu
    if fundamental == True:
        def f(j, c, b, a, d):
            nu = c - (2*b - 2*a) * j - a * j**2 + 4*d * j**3
            return nu
        
    #fit the data to the function.
    fit, pcov = optimize.curve_fit(f, jp, pbranch)
    
    #find the standard error from the covariance matrix
    pcov = np.sqrt(np.diag(pcov))
    c_err, b_err, a_err, d_err = pcov
    #reorder errors to match returned values
    pcov = [a_err, b_err, c_err, d_err]
    
    #unpack the data returned by the fitting function.
    c, b, a, d = fit
    
    #create an array of values for Jp fitted function.
    jp_int = np.linspace(jp[0], jp[-1], 100)
    
    #plot fitted function, and raw values
    plt.plot(jp, pbranch, "b.")
    plt.plot(jp_int, f(jp_int, c, b, a, d))
    
    plt.title("P-Branch Transitions")
    plt.ylabel("$\\tilde{\\nu}$ $cm^{-1}$")
    plt.xlabel("$J_P$")
    
    return a, b, c, d, pcov


def fit_r_branch(rbranch, jr, fundamental = False):
    """Fits a specified function to R-Branch transition wavenumbers to
    their corresponding quantum numbers. Also will plot function and 
    raw data.
    
    The specified function is 
    nu = c+(2*b-5*a-4*d) + (2*b-7*a-12*d)*j - (2*a+12*d)*j**2 - 4*d*j**3
    for the overtone transition, where c = 2*ve - 6*vexe, or 
    nu = c+(2*b-3*a-4*d) + (2*b-4*a-12*d)*j - (a+12*d)*j**2 - 4*d*j**3 
    for the fundamental transition, where c = ve - 2*vexe. 
    
    Parameters
    ----------
    rbranch: list
        List of the wavenumbers where absorption peaks occur for the 
        R-Branch.
    jr: list
        The corresponding quantum numbers for the transitions in 
        rbranch.
    fundamental: Boolean
        Used to determine fitting function to use. True will use the
        function for the fundamental transitions, while False will use
        the function for the first overtone transitions.
        
    Returns
    -------
    a: Float
        Rotational-vibrational coupling constant. Usually denoted a_e.
    b: Float
        Rotational constant. Usually denoted B_e
    c: Float
        Intercept of the function. For the fundemental transition this 
        value is equal to c = ve - 2 * vexe, where ve is the harmonic
        frequency and vexe is the anharmonicity constant.
    d: Float
        Centrifugal Distortion Constant, not to be confused with D_o, 
        which is the Bond Energy.
    pcov: list
        List of corresponding errors at one standard deviation. Found by
        taking the square root of the diagonal of the covariant matrix.
    """
    #choose function to fit to depending on transition
    if fundamental == False:
        def f(j, c, b, a, d):
            nu = c + (2*b - 5*a - 4*d) + (2*b - 7*a - 12*d) * j - (2*a + 12*d) * j**2 - 4*d * j**3
            return nu
    if fundamental == True:
        def f(j, c, b, a, d):
            nu = c + (2*b - 3*a - 4*d) + (2*b - 4*a - 12*d) * j - (a + 12*d) * j**2 - 4*d * j**3
            return nu
    
    #fit the data to the function.
    fit, pcov = optimize.curve_fit(f, jr, rbranch)
    #find the standard error from the covariance matrix
    pcov = np.sqrt(np.diag(pcov))
    c_err, b_err, a_err, d_err = pcov
    #reorder errors to match returned values
    pcov = [a_err, b_err, c_err, d_err]
    
    #unpack the data returned by the fitting function.
    c, b, a, d = fit

    #create an array of values for Jp fitted function.
    jr_int = np.linspace(jr[0], jr[-1], 100)
    
    #plot fitted function, and raw values
    plt.plot(jr, rbranch, "b.")
    plt.plot(jr_int, f(jr_int, c, b, a, d))

    plt.title("R-Branch Transitions")
    plt.ylabel("$\\tilde{\\nu}$ $cm^{-1}$")
    plt.xlabel("$J_R$")
    
    return a, b, c, d, pcov