from AlbroIRSpectroscopy import *

# Read file into the program
absorb, waveln = readSpectrum("AlbroFordham.csv")

# Plot full spectrum
plt.figure(figsize = (20, 4))
plotfull(absorb, waveln)
plt.show()

# Trim full spectrum down to a single transition:
peak1a = absorb[waveindex(5850):waveindex(5450)]
peak1w = waveln[waveindex(5850):waveindex(5450)]

# Plot and label spectrum of a single transition:
plt.figure(figsize = (20, 4))
rb1, jr1, pb1, jp1 = detectpeaks(peak1a, peak1w)
maketitle("HCl", 0, 2)
plt.show()

# Fit Funciton to Wavenumber vs Quantum Number
plt.figure(figsize = (20, 4))
plt.subplot(121)
ar1, br1, cr1, dr1, rerr1 = fit_r_branch(rb1, jr1, False)
plt.subplot(122)
ap1, bp1, cp1, dp1, perr1 = fit_p_branch(pb1, jp1, False)
plt.show()

# Print Out Fitting Parameters
print("a_e(P-Branch):      ", ap1, "\na_e(R-Branch):      ", ar1)
print("B_e(P-Branch):      ", bp1, "\nB_e(R-Branch):      ", br1)
print("2ve-6vexe(P-Branch):", cp1, "\n2ve-6vexe(R-Branch):", cr1)
print("D(P-Branch):        ", dp1, "\nD(R-Branch):        ", dr1)

## Plotted, labelled and fit spectra for other transitions:
peak3a = absorb[waveindex(3100):waveindex(2590)]
peak3w = waveln[waveindex(3100):waveindex(2590)]
plt.figure(figsize = (20, 4))
rb3, jr3, pb3, jp3 = detectpeaks(peak3a, peak3w)
maketitle("HCl", 0, 1)
plt.show()

plt.figure(figsize = (20, 4))
plt.subplot(121)
ar3, br3, cr3, dr3, rerr3 = fit_r_branch(rb3, jr3, True)
plt.subplot(122)
ap3, bp3, cp3, dp3, perr3 = fit_p_branch(pb3, jp3, True)
plt.show()

peak2a = absorb[waveindex(4240):waveindex(3990)]
peak2w = waveln[waveindex(4240):waveindex(3990)]
plt.figure(figsize = (20, 4))
rb2, jr2, pb2, jp2 = detectpeaks(peak2a, peak2w)
maketitle("DCl", 0, 2)
plt.show()

plt.figure(figsize = (20, 4))
plt.subplot(121)
ar2, br2, cr2, dr2, rerr2 = fit_r_branch(rb2, jr2, False)
plt.subplot(122)
ap2, bp2, cp2, dp2, perr2 = fit_p_branch(pb2, jp2, False)
plt.show()


peak4a = absorb[waveindex(2250):waveindex(1875)]
peak4w = waveln[waveindex(2250):waveindex(1875)]
plt.figure(figsize = (20, 4))
rb4, jr4, pb4, jp4 = detectpeaks(peak4a, peak4w)
maketitle("DCl", 0, 1)
plt.show()

plt.figure(figsize = (20, 4))
plt.subplot(121)
ar4, br4, cr4, dr4, rerr4 = fit_r_branch(rb4, jr4, True)
plt.subplot(122)
ap4, bp4, cp4, dp4, perr4 = fit_p_branch(pb4, jp4, True)
plt.show()

# Print out all fitting parameters:
# Uncomment to print all values...
#print(ar1, ar3, ar2, ar4)
#print(ap1, ap3, ap2, ap4)
#print(br1, br3, br2, br4)
#print(bp1, bp3, bp2, bp4)
#print(cr1, cr3, cr2, cr4)
#print(cp1, cp3, cp2, cp4)
#print(dr1, dr3, dr2, dr4)
#print(dp1, dp3, dp2, dp4)