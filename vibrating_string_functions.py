def ClampedString(c1, c2, c3, IC, t, nt, nx, nround=10):
    '''
    # IN #
    
    c1 & c2 & c3: int or float; cf. ClampedString_method()
    
    IC: list of two functions of x; [f(x, t=0), f_t(x, t=0)]
    
    t: int or float; simulation time in seconds
    
    nt: int; number of total time steps from T=0 to T=t-(t/nt)
        [ ht = round(t/nt, nround) ]
    
    nx: int; number of total points on the string with a length equal to one, including the boundary points
        [ hx = round(1/(nx-1), nround) ]
    
    nround: int; number of decimal places for rounding the data
    
    # OUT #
    
    Tdata: numpy.ndarray; a linear space from 0 to (and including) t-ht with nt points
    
    Xdata: numpy.ndarray; a linear space from 0 to (and including) 1 with nx points
    
    Wdata: numpy.ndarray; an array consisting of nt arrays of size nx, corresponding to the data of each time step
    
    # METHOD #
    
    cf. ClampedString_method()
    '''
    import numpy as np
    
    ht = round(t/nt, nround)
    hx = round(1/(nx-1), nround)
    
    Tdata = np.around(np.linspace(0, t-ht, nt), nround)
    Xdata = np.around(np.linspace(0, 1, nx), nround)
    Wdata = np.zeros((nt, nx))
    
    
    if c2==0:
        
        for j in range(nt):  # f(x=0, t) & f(x=nx-1, t)
            Wdata[j][0] = 0
            Wdata[j][nx-1] = 0
        
        for i in range(1, nx-1):  # f(x, t=0)
            Wdata[0][i] = round(IC(Xdata[i])[0], nround)
        
        for i in range(1, nx-1):  # f(x, t=1)
            Wdata[1][i] = round(((2*Wdata[0][i] + 2*ht*IC(Xdata[i])[1] + c1*((ht/hx)**2)*(
                                Wdata[0][i+1] - 2*Wdata[0][i] + Wdata[0][i-1]) + (c3*ht**2)*IC(Xdata[i])[1])/
                                (1 - c3*ht/2))/(1 + (1 + c3*ht/2)/(1 - c3*ht/2)), nround)
        
        for j in range(1, nt-1):  # f(x, t>1)
            for i in range(1, nx-1):
                Wdata[j+1][i] = round((2*Wdata[j][i] - Wdata[j-1][i] + c1*((ht/hx)**2)*(
                                Wdata[j][i+1] - 2*Wdata[j][i] + Wdata[j][i-1]) - (c3*ht/2)*Wdata[j-1][i])/
                                (1 - c3*ht/2), nround)
        
#         print('\nSuccessfully solved!\n')
        return Tdata, Xdata, Wdata
    
    
    else:
    
        for j in range(nt):  # f(x=0, t) & f(x=nx-1, t)
            Wdata[j][0] = 0
            Wdata[j][nx-1] = 0
        
        for i in range(1, nx-1):  # f(x, t=0)
            Wdata[0][i] = round(IC(Xdata[i])[0], nround)
        
        i = 1  # f(x=1, t=1)
        
        Wdata[1][i] = round(((2*Wdata[0][i] + 2*ht*IC(Xdata[i])[1] + ((ht/hx)**2)*(
                        c1*(Wdata[0][i+1] - 2*Wdata[0][i] + Wdata[0][i-1]) + 
                (c2/hx**2)*(Wdata[0][i+2] - 4*Wdata[0][i+1] + 6*Wdata[0][i] - 4*Wdata[0][i-1] + 
                (Wdata[0][i] - 2*hx*0))) + (c3*ht**2)*IC(Xdata[i])[1])/(1 - c3*ht/2))/
                        (1 + (1 + c3*ht/2)/(1 - c3*ht/2)), nround)
        
        for i in range(2, nx-2):  # f(x!=[1, nx-2], t=1)
        
            Wdata[1][i] = round(((2*Wdata[0][i] + 2*ht*IC(Xdata[i])[1] + ((ht/hx)**2)*(
                        c1*(Wdata[0][i+1] - 2*Wdata[0][i] + Wdata[0][i-1]) + 
                (c2/hx**2)*(Wdata[0][i+2] - 4*Wdata[0][i+1] + 6*Wdata[0][i] - 4*Wdata[0][i-1] + Wdata[0][i-2])) + 
                (c3*ht**2)*IC(Xdata[i])[1])/(1 - c3*ht/2))/(1 + (1 + c3*ht/2)/(1 - c3*ht/2)), nround)
        
        i = nx-2  # f(x=nx-2, t=1)
        
        Wdata[1][i] = round(((2*Wdata[0][i] + 2*ht*IC(Xdata[i])[1] + ((ht/hx)**2)*(
                        c1*(Wdata[0][i+1] - 2*Wdata[0][i] + Wdata[0][i-1]) + 
                (c2/hx**2)*((Wdata[0][i] + 2*hx*0) - 4*Wdata[0][i+1] + 6*Wdata[0][i] - 
                4*Wdata[0][i-1] + Wdata[0][i-2])) + (c3*ht**2)*IC(Xdata[i])[1])/(1 - c3*ht/2))/
                        (1 + (1 + c3*ht/2)/(1 - c3*ht/2)), nround)
        
        for j in range(1, nt-1):  # f(x, t>1)
        
            i = 1  # f(x=1, t>1)
        
            Wdata[j+1][i] = round((2*Wdata[j][i] - Wdata[j-1][i] + ((ht/hx)**2)*(
                        c1*(Wdata[j][i+1] - 2*Wdata[j][i] + Wdata[j][i-1]) + 
                (c2/hx**2)*(Wdata[j][i+2] - 4*Wdata[j][i+1] + 6*Wdata[j][i] - 4*Wdata[j][i-1] + 
                (Wdata[j][i] - 2*hx*0))) - (c3*ht/2)*Wdata[j-1][i])/(1 - c3*ht/2), nround)
        
            for i in range(2, nx-2):  # f(x!=[1, nx-2], t>1)
        
                Wdata[j+1][i] = round((2*Wdata[j][i] - Wdata[j-1][i] + ((ht/hx)**2)*(
                        c1*(Wdata[j][i+1] - 2*Wdata[j][i] + Wdata[j][i-1]) + 
                (c2/hx**2)*(Wdata[j][i+2] - 4*Wdata[j][i+1] + 6*Wdata[j][i] - 4*Wdata[j][i-1] + 
                Wdata[j][i-2])) - (c3*ht/2)*Wdata[j-1][i])/(1 - c3*ht/2), nround)
        
            i = nx-2  # f(x=nx-2, t>1)
        
            Wdata[j+1][i] = round((2*Wdata[j][i] - Wdata[j-1][i] + ((ht/hx)**2)*(
                        c1*(Wdata[j][i+1] - 2*Wdata[j][i] + Wdata[j][i-1]) + 
                (c2/hx**2)*((Wdata[j][i] + 2*hx*0) - 4*Wdata[j][i+1] + 6*Wdata[j][i] - 
                4*Wdata[j][i-1] + Wdata[j][i-2])) - (c3*ht/2)*Wdata[j-1][i])/(1 - c3*ht/2), nround)
        
#         print('\nSuccessfully solved!\n')
        return Tdata, Xdata, Wdata



def ClampedString_method():
    '''
    # METHOD #
    
    f_tt = c1 f_xx + c2 f_xxxx + c3 f_t  (c2 & c3 < 0)
    
    x_i := x[i]  (x_i+1 := x[i+1])
    t_j := t[j]
    
    f(x_i, t_j) := f(i, j)
    
    nx = len(x)
    nt = len(t)
    
    Clamped Boundary Condition:  f(0, j) = 0    &    f(nx-1, j) = 0
                                 f_x(0, j) = 0  &    f_x(nx-1, j) = 0
    
    
    if c2 == 0:
    
    
    o ... o     o     o     o     o ... o
                                            .
                                            .
    o ... o     o     o     o     o ... o   .


    o ... o     o     o     o     o ... o   t_j+1


    o ... o     o     o     o     o ... o   t_j


    o ... o     o     o     o     o ... o   t_j-1
    

    o ... o     o     o     o     o ... o   .
                                            .
                                            .
    o ... o     o     o     o     o ... o   t_0
    
    
    o ... o     o     o     o     o ... o   t_-1

    x_0       x_i-1  x_i   x_i+1        x_nx-1
    
    
    x_i+1 - x_i = hx (for all i)  &  t_j+1 - t_j = ht  (for all j)
    
    f_x(i, j) = (f(i+1, j) - f(i-1, j))/(2 hx)
    f_xx(i, j) = (f(i+1, j) - 2 f(i, j) + f(i-1, j)) / hx^2
    
    f_t(i, j) = (f(i, j+1) - f(i, j-1))/(2 ht)
    f_tt(i, j) = (f(i, j+1) - 2 f(i, j) + f(i, j-1))/ht^2
    
    
    (f(i, j+1) - 2 f(i, j) + f(i, j-1))/ht^2 = c1 (f(i+1, j) - 2 f(i, j) + f(i-1, j)) / hx^2 + 
                                                c3 (f(i, j+1) - f(i, j-1))/(2 ht)
    
    ==>  f(i, j+1) = (2 f(i, j) - f(i, j-1) + c1 (ht/hx)^2 (f(i+1, j) - 2 f(i, j) + f(i-1, j)) - 
                                                c3 ht/2 f(i, j-1)) / (1 - c3 ht/2)  (I)
    
    (f(i, 1) - f(i, -1))/(2 ht) = f_t(i, 0)  ==>  f(i, -1) = f(i, 1) - 2 ht * f_t(i, 0)  (II)
    
    
    else:
    

    o     o ... o     o     o     o     o ... o     o
                                                        .
                                                        .
    o     o ... o     o     o     o     o ... o     o   .


    o     o ... o     o     o     o     o ... o     o   t_j+1


    o     o ... o     o     o     o     o ... o     o   t_j


    o     o ... o     o     o     o     o ... o     o   t_j-1
    

    o     o ... o     o     o     o     o ... o     o   .
                                                        .
                                                        .
    o     o ... o     o     o     o     o ... o     o   t_0
    
    
    o     o ... o     o     o     o     o ... o     o   t_-1

  x_-1   x_0       x_i-1  x_i  x_i+1        x_nx-1  x_nx
    
    
    x_i+1 - x_i = hx (for all i)  &  t_j+1 - t_j = ht  (for all j)
    
    f_x(i, j) = (f(i+1, j) - f(i-1, j))/(2 hx)
    f_xx(i, j) = (f(i+1, j) - 2 f(i, j) + f(i-1, j)) / hx^2
    f_xxxx(i, j) = (f(i+2, j) - 4 f(i+1, j) + 6 f(i, j) - 4 f(i-1, j) + f(i-2, j)) / hx^4
    
    f_t(i, j) = (f(i, j+1) - f(i, j-1))/(2 ht)
    f_tt(i, j) = (f(i, j+1) - 2 f(i, j) + f(i, j-1))/ht^2
    
    
    (f(i, j+1) - 2 f(i, j) + f(i, j-1))/ht^2 = c1 (f(i+1, j) - 2 f(i, j) + f(i-1, j)) / hx^2 + 
                                c2 (f(i+2, j) - 4 f(i+1, j) + 6 f(i, j) - 4 f(i-1, j) + f(i-2, j)) / hx^4 + 
                                c3 (f(i, j+1) - f(i, j-1))/(2 ht)
    
    ==>  f(i, j+1) = 2 f(i, j) - f(i, j-1) + (ht/hx)^2 (c1 (f(i+1, j) - 2 f(i, j) + f(i-1, j)) + 
                                c2/hx^2 (f(i+2, j) - 4 f(i+1, j) + 6 f(i, j) - 4 f(i-1, j) + f(i-2, j))) + 
                                c3 ht/2 (f(i, j+1) - f(i, j-1))
    
    ==>  f(i, j+1) = (2 f(i, j) - f(i, j-1) + (ht/hx)^2 (c1 (f(i+1, j) - 2 f(i, j) + f(i-1, j)) + 
                                c2/hx^2 (f(i+2, j) - 4 f(i+1, j) + 6 f(i, j) - 4 f(i-1, j) + f(i-2, j))) - 
                                c3 ht/2 f(i, j-1)) / (1 - c3 ht/2)  (I)
    
    (f(i, 1) - f(i, -1))/(2 ht) = f_t(i, 0)  ==>  f(i, -1) = f(i, 1) - 2 ht * f_t(i, 0)  (II)
    
    f(-1, j) = f(1, j) - 2 hx * f_x(0, j)        (III)
    f(nx, j) = f(nx-2, j) + 2 hx * f_x(nx-1, j)  (IV)
    '''
    return



def StringSound(Tdata, Xdata, Wdata, sources_indices, mic_coor, v_air=1, 
                fmin=0, fmax=22050, hlim=0.05, nround=10):
    '''
    # IN #
    
    Tdata, Xdata, Wdata: cf. ClampedString()
    sources_indices: list or numpy.ndarray; indices of Xdata points acting as sound sources in the microphone environment
    mic_coor: list of two ints or floats; coordinates of the microphone recording the sound
    v_air: int or float; the propagation speed of sound from the string to the microphone
    fmin: int or float; lower frequency limit of the FFT (fast Fourier transform)
    fmax: int or float; upper frequency limit of the FFT
    hlim: 0<=float<=1; setting the lower amplitude limit of each FFT to its maximum amplitude times hlim
    nround: int; number of decimal places for rounding the data
    
    # OUT #
    
    sources: list of numpy.ndarrays; containing the vibrations of each source as an array
    sources_peaks: list of numpy.ndarrays; containing the FFT peaks of each source in the shape of [[frequencies], [amplitudes]]
    mic: numpy.ndarray; the superposition of the propagated vibrations of the sources (sound waves) in the microphone coordinates
    '''
    import numpy as np
    
    nt = len(Tdata)
    nsources = len(sources_indices)
    nft = nt//2 + 1
    
    freqs = np.around(np.fft.rfftfreq(nt, 1/(nt//(Tdata[-1] + Tdata[1]))), nround)
    
    sources = []
    for i, j in enumerate(sources_indices):
        sources.append(Wdata[:, j])
    
    sources_fts = []
    for i in range(nsources):
        sources_fts.append(np.around(np.abs(np.fft.rfft(sources[i])), nround))
    
    sources_peaks = []
    for i in range(nsources):
        source_hlim = np.max(sources_fts[i])*hlim
        indices = np.array([j for j in range(nft) 
                            if (sources_fts[i][j]>=source_hlim and fmin<=freqs[j] and freqs[j]<=fmax)])
        sources_peaks.append(np.array([freqs[indices], sources_fts[i][indices]]))
    
    radii = np.around(((Xdata[sources_indices]-mic_coor[0])**2 + mic_coor[1]**2)**0.5, nround)
    
    mic = np.around(sum([sum([(sources_peaks[j][1][i]/radii[j]**2)*
                                np.cos(2*np.pi*sources_peaks[j][0][i]*(radii[j]/v_air - Tdata)) 
                                for i in range(len(sources_peaks[j][0]))]) 
                                for j in range(nsources)]), nround)
    
    return sources, sources_peaks, mic



def ReconstructSound(mic, sound_f0, sr, c1, fmin=0, fmax=22050, hlim=0.05, 
                      sound_l=3, sound_sr=44100, nround=10):
    '''
    Takes the FFT of the microphone signal (the input 'mic'), extracts its peaks, and 
    accordingly reconstructs a sound signal with the desired fundamental frequency (the input 'sound_f0'), 
    length (the input 'sound_l') and sampling rate (the input 'sound_sr').
    '''
    import numpy as np
    
    nt = len(mic)
    nft = nt//2 + 1
    
    freqs = np.around(np.fft.rfftfreq(nt, 1/sr), nround)
    
    mic_ft = np.around(np.abs(np.fft.rfft(mic)), nround)
    hlim = np.max(mic_ft)*hlim
    indices = np.array([j for j in range(nft) if (mic_ft[j]>=hlim and fmin<=freqs[j] and freqs[j]<=fmax)])
    mic_peaks = np.array([freqs[indices], mic_ft[indices]])

    mic_f = [np.around(mic_peaks[0]*(sound_f0/(c1**0.5/2)), nround), 
             np.around(mic_peaks[1]/1e6, nround)]  # the string's fundamental frequency is (c1**0.5)/2
                                                   # large values of mic_peaks[1] result in distortions
    mic_sound = sines(mic_f, l=sound_l, sr=sound_sr, normalize=False)  
    
    return mic_peaks, mic_f, mic_sound



def sines(f, phi=0, l=3, sr=44100, options=0, normalize=True, nround=10):
    '''
    Generate sinusoidal signals
    
    f: integer or float: signal's frequency
       [[frequencies_1], [amplitudes_1], ...] (could also consist of one element)
    phi: 0: no phase difference
         [phases in radians]
    l: signal's length in seconds
    sr: signal's sampling rate
    options: 0: hearing the components simultaneously
             1: hearing the components individually and successively
    normalize: whether or not to normalize the signal to 1
    nround: number of decimal places for rounding the samples
    '''
    import numpy as np
    
    t = np.around(np.arange(int(l*sr))/sr, nround)
    
    if (type(f)==int or type(f)==float or type(f)==np.float64):
        return np.around(np.sin(2*np.pi*f*t + phi), nround)
    
    n = len(f[0])
    
    if not phi:
        phi = [0]*n
        
    if options==0:
        s = np.zeros(int(sr*l))
        for i in range(n):
            s += f[1][i]*sines(f=f[0][i], phi=phi[i], l=l, sr=sr, options=0, normalize=normalize, nround=nround)
        if normalize:
            return np.around(s/np.max(np.abs(s)), nround)
        else:
            return np.around(s, nround)
    
    elif options==1:
        s = np.zeros(int(sr*l*n))
        for i in range(n):
            seg = f[1][i]*sines(f=f[0][i], phi=phi[i], l=l, sr=sr, options=0, normalize=normalize, nround=nround)
            for j in range(int(sr*l)):
                s[int(i*sr*l)+j] = seg[j]
        if normalize:
            return np.around(s/np.max(np.abs(s)), nround)
        else:
            return np.around(s, nround)



def save(s, name, sr=44100, notice=True):
    '''
    Write wav file
    
    s: mono or stereo signal (to make a stereo signal use 'np.column_stack((s_left, s_right))')
    notice: whether or not to print message
    '''
    from scipy.io import wavfile
    import numpy as np
    
    s = np.array(s)
    wavfile.write(name, sr, s)
    
    if notice:
        print(f'{name} is saved.\n')
    
    return



def adjust(s, als=[0.02]*2, sr=44100, r=1, nround=10):
    '''
    Fade signal in OR out
    
    s: mono signal
    als: [al_beg, al_end]; lengths from the signal's beginning and end (in seconds) to be adjusted according to x^r; >= 0
    r: >= 1; cf. als
    '''
    import numpy as np
    
    s = np.array(s)
    
    ratios = np.array(list(np.linspace(0, 1, int(als[0]*sr))**r) + 
                      [1]*(len(s)-int(sum(als)*sr)) + 
                      list(np.linspace(0, 1, int(als[1]*sr))**r)[::-1])
    
    s *= ratios
    
    return np.around(s, nround)
