### Simulation 1: different x0s; c2 & c3 = 0 ###

'''
To combine the video and audio outputs:
ffmpeg -i strings_simulation_1.mp4 -i strings_simulation_1.wav -c copy -map 0:v:0 -map 1:a:0 strings_simulation_1_combined.mp4
'''

from vibrating_string_functions import *
import numpy as np
from matplotlib import pyplot as plt, animation

# PARAMETERS #

c1 = 16  # v_string**2
c2 = 0
c3 = 0

x0s = [0.05, 0.15, 0.25, 0.35, 0.5]
y0 = 0.015

def IC_tri(x):
    if x<x0:
        return [(y0/x0)*x, 0]
    else:
        return [(y0/(1-x0))-(y0/(1-x0))*x, 0]

sr = 500
nx = 100 + 1
nround = 10

nsteps = len(x0s)  # number of animation steps

nf_cycle = int(sr*(2/c1**0.5))  # it takes the pulse 2/v seconds to travel back to its first place [nf --> ns]
dt_anim = 1  # int; the time of each animation step in seconds (every second is exactly half a cycle)
nt_anim = dt_anim*nsteps*(nf_cycle//2)  # should obviously <= nt
fr_anim = 25  # frame rate of the animation
skip_frames = (nf_cycle//2)//fr_anim  # saves a frame every (skip_frames)th frame

sources_indices = np.arange(5, nx-1, 5)
mic_coor = [0.5, 2]
v_air = 1

sources_fmin = 0
sources_fmax = sr/2  # max = sr/2
sources_hlim = 0.01

mic_fmin = 0
mic_fmax = 15*(c1**0.5)/2
mic_hlim = 0.02

sound_f0 = 200

sound_sr_recon = 44100  # for the reconstructed sound in ReconstructSound()

# in order to have the original sound (with 'sound_f0' fundamental and 'dt_anim' time intervals)
sound_sr_orig = int(nf_cycle*sound_f0)
t = (dt_anim*sound_sr_orig)/sr
nt = int(sr*t)

colors = ['#d7d7d7', '#d44d0d', '#283bb5', '#355c1f']

# WAVE DATA #

Wdatas = []
for i in range(nsteps):
    x0 = x0s[i]
    Tdata, Xdata, Wdata = ClampedString(c1, c2, c3, IC_tri, t, nt, nx, nround)
    Wdatas.append(Wdata)

# SOUND DATA #

mics = []

for i in range(nsteps):
    sources, sources_peaks, mic = StringSound(Tdata, Xdata, Wdatas[i], sources_indices, mic_coor, v_air, 
                                              sources_fmin, sources_fmax, sources_hlim, nround)
    mics.append(mic)

# MIC SPECTRUMS #

mics_peaks = []
mics_sound = []

for i in range(nsteps):
    mic_peaks, mic_f, mic_sound = ReconstructSound(mics[i], sound_f0, sr, c1, mic_fmin, mic_fmax, mic_hlim, 
                                                   dt_anim, sound_sr_recon, nround)
    mics_peaks.append(mic_peaks)
    mics_sound.append(mic_sound)

# ANIMATION #

fig, axes = plt.subplots(nsteps, 2, figsize=(40, 5*nsteps), facecolor=colors[0], 
                                                gridspec_kw={'width_ratios':[1.5, 1]})

fig.suptitle(
    f'$c_1={c1}\;,\;c_2={c2}\;,\;c_3={c3}\;;\;string\_x\in[0,1]\;;\;mic_{{(x,y)}}={tuple(mic_coor)}\;;\;\
        mic_{{hlim}}={mic_hlim}$', fontsize=40)

for i in range(nsteps):
    axes[i][0].axis('off')
    axes[i][0].text(0.5, y0*3, f'$x_0={x0s[i]}$', fontsize=25, horizontalalignment='center')
    axes[i][0].set_xlim(-0.05, 1.05)
    axes[i][0].set_ylim(-.05, .05)

mics_peaks_maxfreq = 0
for peak in mics_peaks:
    if peak[0][-1]>mics_peaks_maxfreq:
        mics_peaks_maxfreq = peak[0][-1]

mics_peaks_maxft = 0
for peak in mics_peaks:
    if np.max(peak[1])>mics_peaks_maxft:
        mics_peaks_maxft = np.max(peak[1])

f0 = c1**0.5/2
if f0==int(f0):
    f0 = int(f0)

axes[0][1].set_title('$FT: amplitude(frequency)$', y=1.1, fontsize=25)

for i in range(nsteps):
    
    axes[i][1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    axes[i][1].tick_params(axis='both', which='major', labelsize=20)
    axes[i][1].yaxis.offsetText.set_fontsize(20)
    axes[i][1].set_ylim(0, (mics_peaks_maxft/mics_peaks_maxft)*1.1)
    axes[i][1].set_xlim(0, mics_peaks_maxfreq+f0)  # extending the axis
    
    j = 1
    while (j-1)*f0<=mics_peaks_maxfreq:  # j-1 in order to add one more harmonic to the plot
        axes[i][1].axvline(j*f0, lw=0.5, c=colors[3])
        j += 1
    
    xticks = [k*f0 for k in range(1, j)]
    axes[i][1].set_xticks(xticks)
    axes[i][1].set_xticklabels(xticks)
    
    axes[i][1].scatter(mics_peaks[i][0], mics_peaks[i][1]/mics_peaks_maxft, marker='o', s=150, c=colors[1])
    
waves = []
for i in range(nsteps):
    waves.append(axes[i][0].plot([], [], c=colors[3], lw=4)[0])

def update(i):
    if i%((nf_cycle//2)*dt_anim)==0:
        j = i//((nf_cycle//2)*dt_anim)
        waves[j].set_color(colors[1])
        if j!=0:
            waves[j-1].set_color(colors[3])
    for j in range(nsteps):
        waves[j].set_data(Xdata, Wdatas[j][i])
    return waves

anim = animation.FuncAnimation(fig=fig, func=update, frames=range(0, nt_anim, skip_frames), interval=1e3/fr_anim)
anim.save('strings_simulation_1.mp4')
print('\nAnimation is saved.\n')
plt.close(fig)

# SOUND #

sound = []
for mic in mics:
    sound += list(adjust(mic, als=[0.02]*2, sr=sound_sr_orig, r=2))
sound = np.array(sound)
sound /= np.max(np.abs(sound))
save(sound, 'strings_simulation_1.wav', sr=sound_sr_orig, notice=False)
print('Audio is saved.\n')
