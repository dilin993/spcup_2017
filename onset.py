 #/usr/bin/env python -W ignore::VisibleDeprecationWarning

from pylab import*
from scipy.io import wavfile
import numpy as np
import math
import warnings
from scipy import ndimage
from scipy import signal
#warnings.filterwarnings("ignore", category=DeprecationWarning) 

def getFFT(s1,sampFreq):
    #n = len(s1) 
    n = 2048;
    p = fft(s1) # take the fourier transform
    nUniquePts = int(ceil((n+1)/2.0))
    p = p[0:nUniquePts]
    p = abs(p)

    p = p / float(n) # scale by the number of points so that

    p = p**2  # square it to get the power 

    if n % 2 > 0: # we've got odd number of points fft
        p[1:len(p)] = p[1:len(p)] * 2
    else:
        p[1:len(p) -1] = p[1:len(p) - 1] * 2 # we've got even number of points fft

    #freqArray = arange(0, nUniquePts, 1.0) * (sampFreq / n);
    #plot(freqArray/1000, 10*log10(p), color='k')
    #xlabel('Frequency (kHz)')
    #ylabel('Power (dB)')

    return p

def mel(nfft,Fs,N_MEL,LOWER,UPPER):
    #LOWER = 300;
    #UPPER = 8000;
    #N_MEL = 10;

    lower_mel = 1125*math.log(1+LOWER/700.0);
    upper_mel = 1125*math.log(1+UPPER/700.0);

    #print lower_mel,upper_mel

    step = (upper_mel-lower_mel)/(N_MEL+1);
    #print step

    m = [lower_mel+i*step for i in range(0,N_MEL+2)]
    #print m

    f = [math.floor((nfft+1)*700*(math.exp(e/1125.0)-1)/Fs) for e in m]
    #print f

    return f


def traingle(lower,middle,upper,nfft):
    filter = np.zeros(nfft/2+1)

    for i in arange(lower,middle+1):
        filter[i] = (i - lower)/float(middle-lower);

    for i in arange(middle,upper+1):
        filter[i] = (upper - i)/float(upper-middle);

    #timeArray = arange(0, nfft, 1)
    #plot(timeArray, filter, color='k')
    #show()
    return filter;



def computeFilterBank(mel_f,signal,N_MEL,nfft):
    filter_out= [];

    for i in arange(0,N_MEL):
        mel_filter = traingle(mel_f[i],mel_f[i+1],mel_f[i+2],nfft)
        fil_out = np.multiply(signal,mel_filter);
        filter_out.append(np.sum(fil_out));

    return filter_out;




def get_spectrogram(signal,sampFreq):
    HOP_SIZE = 512.0;
    WINDOW_SIZE = 2048.0;
    nfft = 2048;
    
    N_MEL = 26;
    LOWER = 300;
    UPPER = 22000;
    hanning = np.hanning(WINDOW_SIZE)

    mel_f = mel(nfft,sampFreq,N_MEL,LOWER,UPPER)


    #i = 0;
    n = len(signal)
    hops = ceil(float(n)/HOP_SIZE)

    spectogram = np.zeros((N_MEL,hops));

    timeArray = [];
    for i in arange(0,hops):
        strip = [];
        if (i*HOP_SIZE+WINDOW_SIZE)<n:
            strip = signal[i*HOP_SIZE:i*HOP_SIZE+WINDOW_SIZE];
            strip = np.multiply(hanning,strip)
            timeArray.append((i*HOP_SIZE+WINDOW_SIZE/2.0)/float(sampFreq));
        else:
            strip = np.zeros(nfft);
            strip_c= signal[i*HOP_SIZE:];
            strip[0:len(strip_c)] = strip_c;
            hanning_cus = np.hanning(len(strip));
            strip = np.multiply(hanning_cus,strip)
            timeArray.append((i*HOP_SIZE+len(strip_c)/2.0)/float(sampFreq));

        fft_strip = getFFT(strip,sampFreq)
        #print len(fft_strip)

        mel_strip = computeFilterBank(mel_f,fft_strip,N_MEL,nfft)
        mel_strip = [x*(2.**15) for x in mel_strip]
        spectogram[:,i] = mel_strip;
    timeArray = np.array(timeArray)
    return spectogram, timeArray

def get_onset(y,fs):
    spectrogram, timeArray = get_spectrogram(y, fs)
    print 'spec: ', spectrogram.shape
    dy = ndimage.sobel(spectrogram, 1)  
    print 'dy: ', dy.shape
    dy[dy<0] = 0 # rectify
    onset = np.sum(dy,axis=0)
    onset = (onset - onset.min()) / onset.max() # normalize
    return onset, timeArray

    

#program Starts here

sampFreq, snd = wavfile.read('open_001.wav')
#snd = snd[:,1]
sampFreq = float(sampFreq)
print sampFreq

snd = snd / (2.**15)


onset, timeArray = get_onset(snd,sampFreq)

peakind = signal.find_peaks_cwt(onset, timeArray)

plot(timeArray,onset)

show()

peaks = timeArray[peakind]
np.savetxt('peaks.txt', peaks)

