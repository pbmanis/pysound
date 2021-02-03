"""
Test sounds and plot waveforms.

This script tests the sound waveform generator for a variety of sounds

"""
import argparse
from pathlib import Path
import sys
import wave  # pythom module
from collections import OrderedDict

import numpy as np
import pyqtgraph as pg
import scipy.signal
from matplotlib import colors as colors
from matplotlib import pyplot as mpl

from pysound import pystim as pystim
from pysound import sound as sound

# define available waveforms:

stims = OrderedDict(
    [
        ("pip", (0, sound.TonePip)),
        ("pipmod", (0, sound.SAMTone)),
        ("noise", (0, sound.NoisePip)),
        ("noisemod", (0, sound.SAMNoise)),
        ("clicks", (0, sound.ClickTrain)),
        ("fmsweep", (0, sound.FMSweep)),
        ("dmr", (0, sound.DynamicRipple)),
        ("ssn", (0, sound.SpeechShapedNoise)),
        ("rss", (0, sound.RandomSpectrumShape)),
    ]
)


def play(args):
    plots = False

    stimarg = args.stimtype
    plots = args.showplot

    PS = pystim.PyStim()

    cf = 2e3
    Fs = PS.out_sampleFreq  # sample frequency
    level = 80.0
    seed = 34978
    fmod = 20.0
    dmod = 20.0

    # reduce to 2 plots only, not whole mess.

    # if plots:
    #     # waveforms
    #     win = pg.GraphicsWindow()
    #     pipwin = win.addPlot(title='sound pip', row=0, col=0)
    #     pipmodwin = win.addPlot(title='100 \% SAM modulated pip', row=1, col=0)
    #     noisewin = win.addPlot(title='WB noise', row=2, col=0)
    #     noisemodwin = win.addPlot(title='100 \% SAM Modulated WB Noise', row=3, col=0)
    #     clickwin = win.addPlot(title='clicks', row=4, col=0)
    #     fmwin = win.addPlot(title='fmsweep', row=5, col=0)
    #     dmrwin = win.addPlot(title='dmr', row=6, col=0)
    #     ssnwin = win.addPlot(title='ssn', row=7, col=0)
    #     rsswin = win.addPlot(title='rss', row=8, col=0)
    #
    #     # spectra
    #     pipwins = win.addPlot(title='sound pip Spec', row=0, col=1)
    #     pipmodwins = win.addPlot(title='100 \% SAM modulated pip', row=1, col=1)
    #     noisewins = win.addPlot(title='WB noise', row=2, col=1)
    #     noisemodwins = win.addPlot(title='100 \% SAM Modulated WB Noise', row=3, col=1)
    #     clickwins = win.addPlot(title='click spec', row=4, col=1)
    #     fmwins = win.addPlot(title='fmsweep spec', row=5, col=1)
    #     dmrwins = win.addPlot(title='dmr spec', row=6, col=1)
    #     ssnwins = win.addPlot(title='ssn', row=7, col=1)
    #     rsswins = win.addPlot(title='rss', row=8, col=1)
    #
    # else:
    #     pipwin = None
    #     pipmodwin = None
    #     noisewin = None
    #     noisemodwin = None
    #     clickwin = None
    #     fmwin = None
    #     dmrwin = None
    #     ssnwin = None
    #     rsswin = None
    #
    #     pipwins = None
    #     pipmodwins = None
    #     noisewins = None
    #     noisemodwins = None
    #     clickwins = None
    #     fmwins = None
    #     dmrwins = None
    #     ssnwins = None
    #     rsswins = None
    #
    # stims = OrderedDict([('pip', (pipwin, sound.TonePip)),
    #                      ('pipmod', (pipmodwin, sound.SAMTone)),
    #                      ('noise', (noisewin, sound.NoisePip)),
    #                      ('noisemod', (noisemodwin, sound.SAMNoise)),
    #                      ('clicks', (clickwin, sound.ClickTrain)),
    #                      ('fmsweep', (fmwin, sound.FMSweep)),
    #                      ('dmr', (dmrwin, sound.DynamicRipple)),
    #                      ('ssn', (ssnwin, sound.SpeechShapedNoise)),
    #                      ('rss', (rsswin, sound.RandomSpectrumShape)),
    #                  ])
    #
    # specs = OrderedDict([('pip', (pipwins, sound.TonePip)),
    #                      ('pipmod', (pipmodwins, sound.SAMTone)),
    #                      ('noise', (noisewins, sound.NoisePip)),
    #                      ('noisemod', (noisemodwins, sound.SAMNoise)),
    #                      ('clicks', (clickwins, sound.ClickTrain)),
    #                      ('fmsweep', (fmwins, sound.FMSweep)),
    #                      ('dmr', (dmrwins, sound.DynamicRipple)),
    #                      ('ssn', (ssnwins, sound.SpeechShapedNoise)),
    #                      ('rss', (rsswins, sound.RandomSpectrumShape)),
    #
    #                  ])
    if stimarg == "all":
        stimlist = stims
    else:
        if stimarg in list(stims.keys()):
            stimlist = [stimarg]
        else:
            raise ValueError("Stimulus %s not in known stimulus types" % stimarg)
    for stim in stimlist:
        print(stim)
        if stim in ["clicks"]:
            soundwave = stims[stim][1](
                rate=Fs,
                duration=1.0,
                dbspl=level,
                click_duration=1e-4,
                click_starts=1e-3 * np.linspace(10, 500, 10),
            )
        elif stim in ["fmsweep"]:
            soundwave = stims[stim][1](
                rate=Fs,
                duration=0.5,
                dbspl=level,
                start=0.0,
                ramp="linear",
                freqs=[16000, 200],
            )
        elif stim in ["pip", "pipmod", "noise", "noisemod"]:
            soundwave = stims[stim][1](
                rate=Fs,
                duration=2.0,
                f0=cf,
                dbspl=level,
                pip_duration=1.8,
                pip_start=[10e-3],
                ramp_duration=2.5e-3,
                fmod=fmod,
                dmod=dmod,
                seed=seed,
            )
        elif stim in ["dmr"]:
            soundwave = stims[stim][1](rate=Fs, duration=10.0)
        elif stim in ["ssn"]:  # speech shaped noise
            waveform = None
            wavfile = Path("wav/testsentence.wav")
            with wave.open(str(wavfile)) as wf:
                nframes = wf.getnframes()
                waveform = np.array(wf.readframes(nframes))
                sr = wf.getframerate()
                nchan = wf.getnchannels()
            Fs = sr
            if waveform is not None:
                audio16 = np.frombuffer(waveform, dtype=np.int16)
                audio = audio16.astype(np.float32)/4096
                soundwave = stims[stim][1](
                    rate=sr, duration=0, waveform=audio, samplingrate=sr
                )
        elif stim in ["rss"]:
            soundwave = stims[stim][1](
                rate=Fs,
                duration=0.5,
                dbspl=level,
                ramp="linear",
                ramp_duration=1e-2,
                f0=4000,
                pip_duration=0.4,
                pip_start=[50e-3],
                amp_group_size=8,
                amp_sd=12,
                spacing=64,
                octaves=3,
            )

        print(("Playing %s" % stim))
        print(np.max(soundwave.sound))
        PS.play_sound(soundwave.sound, soundwave.sound, Fs)

        if plots:  # make one graph for each waveform requested
            fig, ax = mpl.subplots(3, 1, figsize=(8, 10))
            fig.suptitle("Waveform and Spectrum")
            ax = ax.ravel()
            # print wave.time.shape
            # print wave.sound.shape
            ax[0].plot(wave.time, wave.sound)
            f, Pxx_spec = scipy.signal.periodogram(
                wave.sound, Fs
            )  # , window='flattop', nperseg=8192,
            # noverlap=512, scaling='spectrum')
            ax[1].semilogy(f[1:], np.sqrt(Pxx_spec)[1:])
            # ax[1].get_shared_x_axes().join(ax[1], ax[2])
            ax[1].set_xticklabels([])
            nfft = 256
            specfreqs, spectime, Sxx = scipy.signal.spectrogram(
                wave.sound, nperseg=int(0.01 * Fs), fs=Fs
            )
            thr = 1e-8
            Sxx[Sxx <= thr] = thr
            pcm = ax[2].pcolor(
                spectime,
                specfreqs,
                Sxx,
                norm=colors.LogNorm(vmin=Sxx.min(), vmax=Sxx.max()),
                cmap="PuBu_r",
            )
            fig.colorbar(pcm, ax=ax[2], extend="max")
            # Pxx, freqs, bins, im = mpl.specgram(wave.sound, NFFT=nfft, Fs=Fs, noverlap=nfft/4)
            mpl.show()

    # if plots and sys.flags.interactive == 0:
    #      pg.QtGui.QApplication.exec_()
def main():
    parser  = argparse.ArgumentParser(
        description="Play test sounds using the pysound library",
        argument_default=argparse.SUPPRESS,
        fromfile_prefix_chars="@",
    )
    known_stimtypes = list(stims.keys()) # get from the dictionary
    parser.add_argument(dest="stimtype", action="store", default="pip", 
        choices=known_stimtypes,
        help = f"Stimulus types: {str(known_stimtypes):s}",
    )
    parser.add_argument("-p", "--plot", dest="showplot", action="store_true", default=False,
        help="show plots of waveforms")
    args = parser.parse_args()
    play(args)

if __name__ == "__main__":
    main()
