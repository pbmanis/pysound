"""
pysounds: a python Class for interacting with hardware to produce sounds and
record signals.

Output hardware is either an National Instruments DAC card or a system sound
card. If the NI DAC is available, TDT system 3 hardware is assumed as well for
the# attenuators (PA5) and an RP2.1 to input the startle response. Second
channel of RP2.1 is collected as well. Use this for a microphone input to
monitor sound in the chamber. If the system sound card is used, stimuli are
generated and microphone input is collected, but they are not simultaneous. This
is used only for testing.


12/17/2008 Paul B. Manis, Ph.D. UNC Chapel Hill Department of
Otolaryngology/Head and Neck Surgery Supported by NIH Grants DC000425-22 and
DC004551-07 to PBM. Copyright Paul Manis, 2008, 2009

"""

import ctypes
import platform
import struct

import numpy as np
import pyaudio
import scipy
import scipy.signal

osname = platform.system()
nidaq_available = False
if osname == "Windows":
    import win32com.client
    import win32com

    try:
        import nidaq

        nidaq_avaiable = True
    except:
        pass
# The following are reference values for rough calibrations
# They do not correct for system frequency responses

REF_ES_dB = 86.0  # calibration info -  Assumes 10 dB padding with attenuator.
REF_ES_volt = 2.0  # output in volts to get refdb
REF_MAG_dB = 100.0  # right speaker is mag... different scaling.


class Pysounds:
    def __init__(self):
        ################################################################################
        # the first thing we must do is find out what hardware is available and
        # what system we are on.
        ################################################################################
        self.debugFlag = False
        self.hardware = []  # list of hardware found on this system
        self.find_hardware()

    def find_hardware(self):
        if (
            osname != "Windows" or not nidaq_available
        ):  # If not on a Windows system, just set up soundcard
            self.setup_soundcard()
            self.hardware.append("Soundcard")
        else:
            if self.setup_nidaq():
                self.hardware.append("NIDAQ")
            if self.setup_RP21():
                self.hardware.append("RP21")
            if self.setup_PA5():
                self.hardware.append("PA5")

    def setup_soundcard(self):
        if self.debugFlag:
            print("pysounds.init: OS or hardware only supports standard sound card")
        self.hardware.append("pyaudio")
        self.out_sampleFreq = 44100.0
        self.in_sampleFreq = 44100.0

    def setup_nidaq(self, devicename="Dev2"):
        # get the drivers and the activeX control (win32com)
        #            from nidaq import NIDAQ as n

        if self.debugFlag:
            print("pysounds.init: Attempt to Assert num devs > 0:", end=" ")
        if len(nidaq.NIDAQ.listDevices()) <= 0:
            return False

        device = nidaq.NIDAQ.getDevice(devicename)
        if self.debugFlag:
            print("pysounds.init: found nidaq devices.")
            print("devices: %s" % nidaq.NIDAQ.listDevices())
            print("getDevice:", end=" ")
            print("  ", devicename)

            print("\nAnalog Output Channels:", end=" ")
            # print "  AI: ", dev0.listAIChannels()
            print(" AO: ", device.listAOChannels())  # check output only
        return True

    def setup_PA5(self, usbdevices=[1, 2]):
        # active X connection to attenuators
        if osname != "Windows":
            return False
        PA5 = win32com.client.Dispatch("PA5.x")
        for devnum in usbdevices:
            a = PA5.ConnectPA5("USB", usbdevices[devnum])
            if a > 0 and self.debugFlag:
                print("pysounds.init: Connected to PA5 Attenuator 1")
            else:
                print("pysounds.init: Failed to connect to PA5 Attenuator 1")
                return False
            PA5.SetAtten(120.0)
        return True

    def setup_RP21(self):
        if osname != "Windows":
            return False
        self.RP21 = win32com.client.Dispatch("RPco.x")  # connect to RP2.1
        a = self.RP21.ConnectRP2("USB", 1)
        if a > 0 and self.debugFlag:
            print("pysounds.init: RP2.1 Connect is good: %d" % (a))
        else:
            print("pysounds.init: Failed to connect to PA5 Attenuator 1")
            return False
        self.RP21.ClearCOF()
        self.samp_cof_flag = 2  # 2 is for 24.4 kHz
        self.samp_flist = [
            6103.5256125,
            12210.703125,
            24414.0625,
            48828.125,
            97656.25,
            195312.5,
        ]
        if self.samp_cof_flag > 5:
            self.samp_cof_flag = 5
        a = self.RP21.LoadCOFsf(r"C:\pyStartle\startle2.rco", self.samp_cof_flag)
        if a > 0:
            print("pysounds.init: Connected to TDT RP2.1 and startle2.rco is loaded")
        else:
            print("pysounds.init: Error loading startle2.rco?, error = %d" % (a))
            return False
        self.hardware = "nidaq"
        self.out_sampleFreq = self.samp_flist[self.samp_cof_flag]
        self.in_sampleFreq = self.samp_flist[self.samp_cof_flag]
        return True

    def getHardware(self):
        return (self.hardware, self.out_sampleFreq, self.in_sampleFreq)

    # internal debug flag to control printing of intermediate messages
    def debugOn(self):
        self.debugFlag = True

    def debugOff(self):
        self.debugFlag = False

    ################################################################################
    # STIMULUS GENERATION ROUTINES
    #
    # transcribed from Matlab. P. Manis, Nov. 28-December 1 2008.
    ################################################################################

    def StimulusMaker(
        self,
        mode="tone",
        amp=1,
        freq=(1000, 3000, 4000),
        delay=0,
        duration=2000,
        rf=2.5,
        phase0=0,
        samplefreq=44100,
        ipi=20,
        np=1,
        alternate=1,
        level=70,
        playSignal=False,
        plotSignal=True,
        channel=0,
    ):
        # generate a tsound (tone, bb noise, bpnoise)  pip with amplitude (V), frequency (Hz) (or frequencies, using a tuple)
        # delay (msec), duration (msec).
        # if no rf (risefall) time is given (units, msec), cosine^2 shaping with 5 msec ramp duration is applied.
        # if no phase is given, phase starts on 0, with positive slope.
        # level is in dB SPL as given by the reference calibration data above...
        #
        clock = (
            1000.0 / samplefreq
        )  # calculate the sample clock rate, and convert to points per msec (khz)
        uclock = 1000.0 * clock  # microsecond clock
        phi = 2.0 * np.pi * phase0 / 360.0  # convert phase from degrees to radians...
        Fs = 1000 / clock
        phi = 0  # actually, always 0 phase for start
        w = []
        fil = self.rfShape(
            0, duration, samplefreq, rf
        )  # make the shape filter with 0 delay
        jd = int(np.floor(delay / clock))  # beginning of signal buildup (delay time)
        if jd < 0:
            jd = 0
        jpts = np.arange(0, len(fil))
        signal = np.zeros(len(jpts))
        siglen = len(signal)

        if mode == "tone":
            for i in range(0, len(freq)):
                signal = signal + fil * amp * np.sin(2.0 * np.pi * freq[i] * jpts / Fs)
                if self.debugFlag:
                    print("Generated Tone at %7.1fHz" % (freq[i]))

        if mode == "bbnoise":
            signal = signal + fil * amp * np.random.normal(0, 1, siglen)
            if self.debugFlag:
                print("BroadBand Noise ")

        if mode == "bpnoise":
            tsignal = fil * amp * np.random.normal(0, 1, siglen)
            # use freq[0] and freq[1] to set bandpass on the noise
            if self.debugFlag:
                print("freqs: HP: %6.1f    LP: %6.1f" % (freq[0], freq[1]))
            sf2 = samplefreq * 2.0  # nyquist limit
            if freq[0] > sf2 or freq[1] > sf2:
                print("freqs: ", freq)
                print("nyquist limit: ", sf2)
                print("sample frequ: ", samplefreq)
                print("coefficients not bounded [0, 1] for w... ")
                return np.array(signal)
            wp = [float(freq[0]) / sf2, float(freq[1]) / sf2]
            ws = [0.75 * float(freq[0]) / sf2, 1.25 * float(freq[1]) / sf2]
            (filter_b, filter_a) = scipy.signal.iirdesign(
                wp, ws, gpass=2.0, gstop=60.0, ftype="ellip"
            )
            if self.debugFlag:
                print("BandPass Noise %7.1f-%7.1f" % (freq[0], freq[1]))
            signal = scipy.signal.lfilter(filter_b, filter_a, tsignal)

        if mode == "notchnoise":
            return np.array(signal)

        if mode == "multitones":
            return np.array(signal)

        if mode == "silence":
            return np.array(signal)

        # now build the waveform from the components
        w = np.zeros(np.ceil(ipi * (np - 1) / clock) + jd + siglen)
        sign = np.ones(np)
        if alternate is True:
            sign[list(range(1, np, 2))] = -1
        id = int(floor(ipi / clock))
        for i in range(0, np):  # for each pulse in the waveform
            j0 = jd + i * id  # compute start time
            w[list(range(j0, j0 + siglen))] = sign[i] * signal

        w = w * self.dbconvert(
            spl=level, chan=channel
        )  # aftera all the shaping ane scaling, we convert to generate a signal of w dB
        if playSignal is True:
            self.playSound(w, w, samplefreq)

        if plotSignal is True:
            self.plotSignal(w, w, clock)
        return array(w)

    #
    # Rise-fall shaping of a waveform. This routine generates an envelope with
    # 1 as the signal max, and 0 as the baseline (off), with cosine^2 shaping of
    # duration rf starting at delay (msec). The duration of the signal includes the
    # rise and fall, so the duration of the signal at full amplitude is dur - 2*rf.
    # Note that since samplefreq is in Hz, delya, rf and duratio are converted to
    # seconds from the msec in the call.
    def rfShape(self, delay=0, duration=100, samplefreq=44100, rf=2.5):
        jd = int(
            np.floor((delay / 1000.0) * samplefreq)
        )  # beginning of signal buildup (delay time)
        if jd < 0:
            jd = 0
        je = int(
            np.floor(((delay + duration) / 1000.0) * samplefreq)
        )  # end of signal decay (duration + delay)
        #
        # build sin^2 filter from 0 to 90deg for shaping the waveform
        #
        nf = int(np.floor((rf / 1000.0) * samplefreq))  # number of points in the filter
        fo = 1.0 / (
            4.0 * (rf / 1000.0)
        )  # filter "frequency" in Hz - the 4 is because we use only 90deg for the rf component

        pts = np.arange(jd, jd + nf)
        fil = np.zeros(je)
        fil[list(range(jd, jd + nf))] = (
            np.sin(2 * np.pi * fo * pts * samplefreq) ** 2
        )  # filter
        fil[list(range(jd + nf, je - nf))] = 1
        pts = list(range(je - nf, je))
        kpts = list(range(jd + nf, jd, -1))
        fil[pts] = fil[kpts]
        return fil

    #
    # insertGap takes a waveform and inserts a shaped gap into it.
    # currently, gap is all the way off, i.e., 0 intensity.
    # a future change is to include relative gap level (-dB from current waveform)
    #
    def insertGap(self, wave, delay=20, duration=20, rf=2.5, samplefreq=44100):
        fil = self.rfShape(
            delay, duration, samplefreq, rf
        )  # make the shape filter with 0 delay
        lenf = len(fil)
        lenw = len(wave)
        if lenw > lenf:
            fil = np.append(fil, np.zeros(lenw - lenf))
        if lenf > lenw:
            fil = np.append(fil, np.zeros(lenf - lenw))
        return wave * (1.0 - fil)

    #
    # compute voltage from reference dB level
    # db = 20 * log10 (Vsignal/Vref)
    #
    def dbconvert(self, spl=0, chan=0):
        ref = REF_ES_dB
        if chan == 1:
            ref = REF_MAG_dB

        zeroref = REF_ES_volt / (10 ** (ref / 20.0))
        sf = zeroref * 10 ** (spl / 20.0)
        # actually, the voltage needed to get spl out...
        if self.debugFlag:
            print("pysounds.dbconvert: scale = %f for %f dB" % (sf, spl))
        return sf  # return a scale factor to multiply by a waveform normalized to 1

    ################################################################################
    # hardware interactions:
    #
    # set the attenuators on the PA5.
    # If no args are given, set to max attenuation

    def setAttens(self, atten_left=120.0, atten_right=120.0):
        if self.hardware == "nidaq":
            PA5.ConnectPA5("USB", 1)
            PA5.SetAtten(attenl)
            PA5.ConnectPA5("USB", 2)
            PA5.SetAtten(attenr)

    #
    # playSound sends the sound out to an audio device. In the absence of NI card
    # and TDT system, it will use the system audio device (sound card, etc)
    # The waveform is played in stereo.
    # Postduration is given in seconds...
    def play_sound(self, wavel, waver, samplefreq, postduration=0.35):
        if "pyaudio" in self.hardware:
            self.audio = pyaudio.PyAudio()
            chunk = 1024
            FORMAT = pyaudio.paFloat32
            CHANNELS = 2
            RATE = samplefreq
            if self.debugFlag:
                print("pysounds.play_sound: samplefreq: %f" % (RATE))
            self.stream = self.audio.open(
                format=FORMAT,
                channels=CHANNELS,
                rate=int(RATE),
                output=True,
                input=True,
                frames_per_buffer=chunk,
            )
            # play stream
            # print self.stream
            wave = np.zeros(2 * len(wavel))
            if len(wavel) != len(waver):
                print(
                    "pysounds.play_sound: waves not matched in length: %d vs. %d (L,R)"
                    % (len(wavel), len(waver))
                )
                return
            (waver, clipr) = self.clip(waver, 20.0)
            (wavel, clipl) = self.clip(wavel, 20.0)
            wave[0::2] = waver
            wave[
                1::2
            ] = wavel  # order chosen so matches entymotic earphones on my macbookpro.
            postdur = int(float(postduration * self.in_sampleFreq))
            rwave = self.read_array(len(wavel) + postdur, CHANNELS)
            self.write_array(wave)
            self.stream.stop_stream()
            self.stream.close()
            self.audio.terminate()
            self.ch1 = rwave[0::2]
            self.ch2 = rwave[1::2]
            return

        if "NIDAQ" in self.hardware:
            self.task = dev0.createTask()  # creat a task for the NI 6731 board.
            self.task.CreateAOVoltageChan(
                "/Dev2/ao0", "ao0", -10.0, 10.0, nidaq.Val_Volts, None
            )
            self.task.CreateAOVoltageChan(
                "/Dev2/ao1", "ao1", -10.0, 10.0, nidaq.Val_Volts, None
            )  # use 2 channels
            wlen = 2 * len(wavel)
            self.task.CfgSampClkTiming(
                None, samplefreq, nidaq.Val_Rising, nidaq.Val_FiniteSamps, len(wavel)
            )
            # DAQmxCfgDigEdgeStartTrig (taskHandle, "PFI0", DAQmx_Val_Rising);
            self.task.SetStartTrigType(nidaq.Val_DigEdge)
            self.task.CfgDigEdgeStartTrig("PFI0", nidaq.Val_Rising)
            daqwave = np.zeros(wlen)
            (wavel, clipl) = self.clip(wavel, 10.0)
            (waver, clipr) = self.clip(waver, 10.0)

            daqwave[0 : len(wavel)] = wavel
            daqwave[len(wavel) :] = waver
            # concatenate channels (using "groupbychannel" in writeanalogf64)
            dur = wlen / float(samplefreq)
            self.task.write(daqwave)
            if "RP21" in self.hardware:
                # now take in some acquisition...
                a = self.RP21.ClearCOF()
                if a <= 0:
                    print("pysounds.playSound: Unable to clear RP2.1")
                    return
                a = self.RP21.LoadCOFsf(
                    r"C:\pyStartle\startle2.rco", self.samp_cof_flag
                )
                if a > 0 and self.debugFlag:
                    print(
                        "pysounds.playSound: Connected to TDT RP2.1 and startle2.rco is loaded"
                    )
                else:
                    print(
                        "pysounds.playSound: Error loading startle2.rco?, error = %d"
                        % (a)
                    )
                    hwerr = 1
                    return
                self.trueFreq = self.RP21.GetSFreq()
                Ndata = ceil(0.5 * (dur + postduration) * self.trueFreq)
                self.RP21.SetTagVal(
                    "REC_Size", Ndata
                )  # old version using serbuf  -- with
                # new version using SerialBuf, can't set data size - it is fixed.
                # however, old version could not read the data size tag value, so
                # could not determine when buffer was full/acquisition was done.
            if "PA5" in self.hardware:
                self.setAttens(10.0, 10.0)  # set equal, but not at minimum...

            self.task.start()  # start the NI AO task
            if "RP21" in self.hardware:
                a = self.RP21.Run()  # start the RP2.1 processor...
                a = self.RP21.SoftTrg(
                    1
                )  # and trigger it. RP2.1 will in turn start the ni card
            while not self.task.isTaskDone():  # wait for AO to finish?
                if not self.PPGo:  # while waiting, check for stop.
                    if "RP21" in self.hardware:
                        self.RP21.Halt()
                    self.task.stop()
                    return
            self.task.stop()  # done, so stop the output.
            if "PA5" in self.hardware:
                self.setAttens()  # attenuators down (there is noise otherwise)
            # read the data...
            curindex1 = self.RP21.RP21.GetTagVal("Index1")
            curindex2 = self.RP21.GetTagVal("Index2")
            if "RP21" in self.hardware:
                while (
                    curindex1 < Ndata or curindex2 < Ndata
                ):  # wait for input data to be sampled
                    if not self.PPGo:  # while waiting, check for stop.
                        self.RP21.Halt()
                        return
                    curindex1 = self.RP21.GetTagVal("Index1")
                    curindex2 = self.RP21.GetTagVal("Index2")
            self.task.stop()
            if "RP21" in self.hardware:
                self.ch2 = self.RP21.ReadTagV("Data_out2", 0, Ndata)
                # ch2 = ch2 - mean(ch2[1:int(Ndata/20)]) # baseline: first 5% of trace
                self.ch1 = self.RP21.ReadTagV("Data_out1", 0, Ndata)
                self.RP21.Halt()

    def retrieveInputs(self):
        return (self.ch1, self.ch2)

    def HwOff(self):  # turn the hardware off if you can.
        if "Soundcard" in self.hardware:
            self.stream.stop_stream()
            self.stream.close()
            self.audio.terminate()

        if "nidaq" in self.hardware:
            self.task.stop()
            self.setAttens()
            if "RP21" in self.hardware:
                self.RP21.Halt()

    # clip data to max value (+/-) to avoid problems with daqs
    def clip(self, data, maxval):
        if self.debugFlag:
            print(
                "pysounds.clip: max(data) = %f, %f and maxval = %f"
                % (max(data), min(data), maxval)
            )
        clip = 0
        u = np.where(data >= maxval)
        ul = list(np.transpose(u).flat)
        if len(ul) > 0:
            data[ul] = maxval
            clip = 1  # set a flag in case we want to know
            if self.debugFlag:
                print("pysounds.clip: clipping %d positive points" % (len(ul)))
        minval = -maxval
        v = np.where(data <= minval)
        vl = list(np.transpose(v).flat)
        if len(vl) > 0:
            data[vl] = minval
            clip = 1
            if self.debugFlag:
                print("pysounds.clip: clipping %d negative points" % (len(vl)))
        if self.debugFlag:
            print(
                "pysounds.clip: clipped max(data) = %f, %f and maxval = %f"
                % (np.max(data), np.min(data), maxval)
            )
        return (data, clip)

    ################################################################################
    # the following was taken from #http://hlzr.net/docs/pyaudio.html
    # it is used for reading and writing to the system audio devie
    #
    ################################################################################
    def write_array(self, data):
        """
        Outputs a numpy array to the audio port, using PyAudio.
        """
        # Make Buffer
        buffer_size = struct.calcsize("@f") * len(data)
        output_buffer = ctypes.create_string_buffer(buffer_size)

        # Fill Up Buffer
        # struct needs @fffff, one f for each float
        format = "@" + "f" * len(data)
        struct.pack_into(format, output_buffer, 0, *data)

        # Shove contents of buffer out audio port
        self.stream.write(output_buffer)

    def read_array(self, size, channels=1):
        input_str_buffer = np.zeros((size, 1))  # self.stream.read(size)
        input_float_buffer = struct.unpack(
            "@" + "f" * size * channels, input_str_buffer
        )
        return np.array(input_float_buffer)


if __name__ == "__main__":
    P = Pysounds()
    print(P.hardware)
