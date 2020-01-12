import AlbroIRSpectroscopy

absorb, waveln = readSpectrum("AlbroFordham.csv")   

plt.figure(figsize = (20, 4))
plotfull(absorb, waveln)
plt.show()

peak1a = absorb[waveindex(5850):waveindex(5450)]
peak1w = waveln[waveindex(5850):waveindex(5450)]
plt.figure(figsize = (20, 4))
rb1, jr1, pb1, jp1 = detectpeaks(peak1a, peak1w, height = 0.432, distance = 40)
maketitle("HCl", 0, 2)
plt.show()

peak3a = absorb[waveindex(3100):waveindex(2590)]
peak3w = waveln[waveindex(3100):waveindex(2590)]
plt.figure(figsize = (20, 4))
rb3, jr3, pb3, jp3 = detectpeaks(peak3a, peak3w, height = 0.55, distance = 40)
maketitle("HCl", 0, 1)
plt.show()

peak2a = absorb[waveindex(4240):waveindex(3990)]
peak2w = waveln[waveindex(4240):waveindex(3990)]
plt.figure(figsize = (20, 4))
rb2, jr2, pb2, jp2 = detectpeaks(peak2a, peak2w, height = 0.46, distance = 40)
maketitle("DCl", 0, 2)
plt.show()

peak4a = absorb[waveindex(2250):waveindex(1875)]
peak4w = waveln[waveindex(2250):waveindex(1875)]
plt.figure(figsize = (20, 4))
rb4, jr4, pb4, jp4 = detectpeaks(peak4a, peak4w, height = 0.7, distance = 40)
maketitle("DCl", 0, 1)
plt.show()