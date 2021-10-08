#!/usr/bin/env python


"""
pystim: a Python Class for interacting with hardware to produce sounds and
record signals.

Output hardware is either an National Instruments DAC card or a system sound card
If the NI DAC is available, TDT system 3 hardware is assumed as well for the
attenuators (PA5) and an RP2.1.

If the system sound card is used, stimuli are generated. This is used only for testing.


12/17/2008 Paul B. Manis, Ph.D.
UNC Chapel Hill
Department of Otolaryngology/Head and Neck Surgery
Supported by NIH Grants DC000425-22 and DC004551-07 to PBM.
Copyright Paul Manis, 2008, 2009

Refactored and modified version, includes access to rz5d to help synchronize acquisition.
August, 2017
"""


"""
TDT manual: 
Sweep Control
To use the sweep control circuit constructs the following names are required:
zSwPeriod: The period of the sweep duration. This is set in OpenWorkbench
and can not be modified during block acquisition.
If it is necessary to change this value during the experiment, an
*asynchronous next sweep control circuit* construct should be used
 See Asynchronous Next Sweep Control, page 324 for more information.
317
OpenEx User's Guide
318
zSwCount: The maximum number of sweeps before the signal is terminated.
If this requires manual or external control, the value should be set to -1 through the OpenWorkbench protocol.

"""
import ctypes
import os
import platform
import struct
import time

import numpy as np

from pysound import sound as sound

opsys = platform.system()
nidaq_available = False
if opsys in ["nt", "Windows"]:
    try:
        import nidaqmx
        import tdt
        # import nidaq
        import win32com.client
        from nidaqmx.constants import AcquisitionType, Edge, VoltageUnits
        nidaq_available = True
    except:
        pass

if opsys in ["Darwin", "Linux"] or nidaq_available == False:
    import pyaudio

# The following are reference values for rough calibrations
# They do not correct for system frequency responses
# They are old.

REF_ES_dB = 86.0  # calibration info -  Assumes 10 dB padding with attenuator.
REF_ES_volt = 2.0  # output in volts to get refdb
REF_MAG_dB = 100.0  # right speaker is mag... different scaling.

RZ5D_Idle = 0
RZ5D_Preview = 2
RZ5D_Standby = 1
RZ5D_Run = 3


class PyStim:
    def __init__(self, hdw=["Soundcard"], devicename="dev1"):
        """ 
        During initialization, we identify what hardware is available.
        
        Parameters
        ----------
        hdw : list : (Default: ['Soundcard'])
            A list of the names of devices we expect to be able to use
            For example: ['PA5', 'NIDAQ', 'RZ5D'] for an attenuator, an NI
            card (for DAC output) and the TDT RZ5D DSP unit.
        devicename : str (Default: 'dev1')
            The device name for the NI device we will use.
        """
        self.index = 0
        self.debugFlag = False
        self.task = None  # NI Task
        self.required_hardware = hdw  # Require specific hardware
        self.hardware = []  # list of hardware actually found on this system
        self.find_hardware(
            device_info={"devicename": devicename}
        )  # population the self.hardware list
        self.TankName = []

    def find_hardware(self, device_info=None):
        """
        Find the hardware on the system.
        For non-windows systems, this just finds the system soundcard for testing
        Otherwise it looks for the requested hardware.
        Keeps track of available hardware in the self.hardware list
        
        Parameters
        ----------
        device_info : dict (default: None)
            User defined dictionary of parameters to pass to devices
        
        """
        if opsys in ["Darwin", "Linux"] or nidaq_available is False:  # If not on a Windows system, just set up soundcard
            self.setup_soundcard()
            self.hardware.append("Soundcard")
            self.out_samplefreq = 44100
        else:
            if "NIDAQ" in self.required_hardware and self.setup_nidaq(device_info):
                self.hardware.append("NIDAQ")
            if "RP21" in self.required_hardware and self.setup_RP21(
                "c:\\TDT\\OpenEx\\MyProjects\\Tetrode\\RCOCircuits\\tone_search.rcx"
            ):
                self.hardware.append("RP21")
            if "PA5" in self.required_hardware and self.setup_PA5():
                self.hardware.append("PA5")
            if "RZ5D" in self.required_hardware and self.setup_RZ5D():
                self.hardware.append("RZ5D")

    def setup_soundcard(self):
        if self.debugFlag:
            print("pysounds.init: OS or hardware only supports standard sound card")
        self.hardware.append("pyaudio")
        self.out_sampleFreq = 44100.0
        self.in_sampleFreq = 44100.0

    def setup_nidaq(self, device_info):
        # get the drivers and the activeX control (win32com)

        # if len(nidaq.NIDAQ.listDevices()) <=  0:
        #     if self.requireNIPA5:
        #         raise IOError('NIDAQ requirement requested, but device not found')
        #     else:
        #         return False
        self.NIDevice = nidaqmx.system.System.local()
        self.NIDevicename = self.NIDevice.devices.device_names
        # self.NIDevice = nidaq.NIDAQ.getDevice(device_info['devicename'])
        # self.NIDevicename = device_info['devicename']
        self.out_sampleFreq = 200000
        return True

    def show_nidaq(self):
        print("pysounds.init: found nidaq devices.")
        print("devices: %s" % self.NIDevice.devices.device_names)
        # print ("devices: %s" % nidaq.NIDAQ.listDevices())
        print("getDevice: " % self.NIDevice)
        print(
            "\nAnalog Output Channels: %d"
            % self.NIDevice.devices[self.NIDevicename].ao_physical_chans.channel_names
        )
        # print ("\nAnalog Output Channels: %d" %  self.NIDevice.listAOChannels())

    def setup_PA5(self, devnum=1):
        # active X connection to attenuators
        self.PA5 = win32com.client.Dispatch("PA5.x")
        a = self.PA5.ConnectPA5("USB", devnum)
        if a > 0:
            if self.debugFlag:
                print("pysounds.init: Connected to PA5 Attenuator %d" % devnum)
        else:
            if "PA5" in self.required_hardware:
                raise IOError("PA5 requirement requested, but device not found")
            else:
                return False
        self.PA5.SetAtten(120.0)
        return True

    def setup_RP21(self, rcofile):
        self.RP21_rcofile = rcofile
        self.RP21 = win32com.client.Dispatch("RPco.x")  # connect to RP2.1
        a = self.RP21.ConnectRP2("USB", 1)
        if a > 0 and self.debugFlag:
            print("pystim.setup_RP21: RP2.1 Connect is good: %d" % (a))
        else:
            print("pystim.setup_RP21: Failed to connect to RP2.1")
            return False
        self.RP21.ClearCOF()
        self.samp_cof_flag = 5  # 2 is for 24.4 kHz
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
        a = self.RP21.LoadCOFsf(rcofile, self.samp_cof_flag)
        if a > 0:
            print(
                "pystim.setup_RP21: File %s loaded\n      and sample rate set to %f"
                % (rcofile, self.samp_fllist[self.camp_cof_flag])
            )
        else:
            print("pysounds.init: Error loading RCO file %s, error = %d" % (rcofile, a))
            return False
        self.out_sampleFreq = self.samp_flist[self.samp_cof_flag]
        self.in_sampleFreq = self.samp_flist[self.samp_cof_flag]
        return True

    def show_RP21(self):
        """
        TODO: report RP2.1 info: cof rate, loaded circuit, sample freqs
        """
        pass

    def setup_RZ5D(self):
        self.RZ5D = tdt.SynapseAPI()
        # self.RZ5D = win32com.client.Dispatch('TDevAcc.X')
        # self.RZ5D.ConnectServer('Local')

        # # RZ5D Parameter tag definitions - mostly from the TDT CoreSweep macro
        # self.RZ5D_ParTags = {'SweepPeriod': 'ACQ_16ch.zSwPeriod',
        #                      'SweepTrigger': 'ACQ_16Ch.SweepTrigger',
        #                      'TotalSweepCount': 'ACQ_16ch.zSwCount',
        #                      'CurrentSweep': 'ACQ_16ch.zSwNum',
        #                      'SweepDone': 'ACQ_16ch.zSwDone',
        #                      }
        # self.getParams_RZ5D()

        # self.RZ5DParams['SampleFrequency'] = self.RZ5D.GetDeviceSF(self.RZ5DParams['device_name']) # get device sample frequency
        if self.RZ5D.getModeStr() != "Idle":
            self.RZ5D.setModeStr("Idle")
        # self.RZ5D.SetSysMode(RZ5D_Standby) # Standby needed to set up parameters....
        #     # time.sleep(2.0)
        # time.sleep(0.5)

        # self.RZ5D.SetTargetVal(self.RZ5D_ParTags['TotalSweepCount'], 3)

        # self.RZ5D.SetTargetVal(self.RZ5D_ParTags['SweepPeriod'], 1.0*self.RZ5D.GetDeviceSF(self.RZ5DParams['device_name'])) # initially set for one second
        # self.RZ5DParams['zSwPeriod'] = self.RZ5D.GetTargetVal(self.RZ5D_ParTags['SweepPeriod'])
        return True

    def getParams_RZ5D(self):
        self.RZ5DParams = {}  # keep a local copy of the parameters
        self.RZ5DParams["device_name"] = self.RZ5D.getGizmoNames()
        self.RZ5DParams["device status"] = self.RZ5D.getModeStr()
        # self.RZ5DParams['device_name'] = self.RZ5D.GetDeviceName(0)
        # self.RZ5DParams['RCO'] = self.RZ5D.GetDeviceRCO(self.RZ5DParams['device_name'])
        # self.RZ5DParams['device_status'] = self.RZ5D.GetDeviceStatus(self.RZ5DParams['device_name'])
        # self.RZ5DParams['zSwCount'] = self.RZ5D.GetTargetVal(self.RZ5D_ParTags['TotalSweepCount'])

    def show_RZ5D(self):
        # print('Device is using RCO/X file: {0:s}'.format(self.RZ5DParams['RCO']))
        print("Device Status: {0:d}".format(self.RZ5DParams["device_status"]))
        # print('RZ5D Sample Frequency: %f' % self.RZ5DParams['SampleFrequency'])
        # print('zSwCount: %d' % self.RZ5DParams['zSwCount'])
        # print('zSwPeriod points (N): {0:d}'.format(int(self.RZ5DParams['zSwPeriod'])))
        # for tag in [68, 73, 76, 80, 83, 65]:
        #     self.getTags(self.RZ5D, self.RZ5DParams['device_name'], tag)

    # def getTags(self, device, device_name, tagnum):
    #     tag =  device.GetNextTag(device_name, tagnum, 1)
    #     if len(tag) > 0:
    #         print('Type {0:2d}: Tag = {1:s}'.format(tagnum, tag))
    #     while len(tag) > 0:
    #         tag =  device.GetNextTag(device_name, tagnum, 0)
    #         if len(tag) == 0:
    #             return
    #         print('         Tag = {0:s}'.format(tag))

    def present_stim(
        self,
        waveforms,
        stimulus_period=1.0,
        reps=1,
        runmode="Record",
        protocol="Search",
    ):
        # sf = self.RZ5D.GetDeviceSF(self.RZ5DParams['device_name'])
        # # self.RZ5D.SetSysMode(RZ5D_Preview) # Preview needed to set up parameters....
        # self.RZ5D.setTargetVal(self.RZ5D_ParTags['SweepPeriod'], stimulus_period*sf)
        # self.RZ5D.setTargetVal(self.RZ5D_ParTags['TotalSweepCount'], reps+1)
        if self.RZ5D.getModeStr() != runmode:
            #TFR 20101007 removing the block checking/matching/storing stuff- important! This needs to be revised for record mode!
            # if runmode == "Record":
            #     protocol.replace(" ", "")
            #     subject = self.RZ5D.getCurrentSubject()
            #     self.TankName = self.RZ5D.getCurrentTank()
            #     # tankname=self.RZ5D.getCurrentTank()
            #     newblock = subject + protocol + "{:03d}".format(self.index)
            #     for checkBlocks in self.RZ5D.getKnownBlocks():
            #         if os.path.join(self.TankName, newblock) in checkBlocks:
            #             self.index = self.index + 1
            #             newblock = subject + protocol + "{:03d}".format(self.index)
            #     self.RZ5D.setCurrentBlock(newblock)
            #     self.index = self.index + 1
            self.RZ5D.setModeStr(runmode)

            time.sleep(1.0)
        self.prepare_NIDAQ(waveforms)  # load up NIDAQ to go
        # time.sleep(2.0) # just wait a few msec #added two seconds here 20180306

        # self.RZ5D.SetSysMode(runmode)

    def RZ5D_close(self):
        # self.RZ5D.SetSysMode(RZ5D_Idle) # Idle
        if self.RZ5D.getModeStr() != "Idle":
            self.RZ5D.setModeStr("Idle")
        time.sleep(1.0)
        # self.RZ5D.CloseConnection();

    def getHardware(self):
        return (self.hardware, self.out_sampleFreq, self.in_sampleFreq)

    # internal debug flag to control printing of intermediate messages
    def debugOn(self):
        self.debugFlag = True

    def debugOff(self):
        self.debugFlag = False

    def dbconvert(self, spl=0, chan=0):
        """
        compute voltage from reference dB level
        db = 20 * log10 (Vsignal/Vref)
        """
        ref = REF_ES_dB
        if chan == 1:
            ref = REF_MAG_dB

        zeroref = REF_ES_volt / (10 ** (ref / 20.0))
        sf = zeroref * 10 ** (spl / 20.0)
        # actually, the voltage needed to get spl out...
        if self.debugFlag:
            print("pystim.dbconvert: scale = %f for %f dB" % (sf, spl))
        return sf  # return a scale factor to multiply by a waveform normalized to 1

    def setAttens(self, atten_left=120.0, atten_right=None):
        if "PA5" in self.hardware:
            self.PA5.ConnectPA5("USB", 1)
            self.PA5.SetAtten(atten_left)
            if atten_right is not None:
                self.PA5.ConnectPA5("USB", 2)
                self.PA5.SetAtten(atten_right)

    def play_sound(
        self,
        wavel,
        waver=None,
        samplefreq=44100,
        postduration=0.05,
        attns=[20.0, 20.0],
        isi=1.0,
        reps=1,
        protocol="Search",
        storedata=True,
    ):
        """
        play_sound sends the sound out to an audio device.
        In the absence of NI card, and TDT system, we use the system audio device (sound card, etc)
        The waveform is played in both channels on sound cards, possibly on both channels
        for other devices if there are 2 channels.
        
        Parameters
        ----------
        wavel : numpy array of floats
            Left channel waveform
        waver : numpy of floats
            Right channel waveform
        samplefreq : float
            output sample frequency (Hz)
        postduration : float (default: 0.35)
            Time after end of stimulus, in seconds
        attns : 2x1 list (default: [20., 20.])
            Attenuator settings to use for this stimulus
        
        """
        if storedata:
            runmode = "Record"
        else:
            runmode = "Preview"
        # create an output waveform that has the stimulus repeated reps times with the selected ISI
        samplefreq = self.out_sampleFreq
        stimulus_duration = isi * reps  # len(wavel)*samplefreq + postduration
        pts_per_rep = int(float(isi) * samplefreq)
        if wavel.shape[0] < pts_per_rep:
            wavel = np.concatenate(
                (wavel, np.zeros(pts_per_rep - wavel.shape[0])), axis=0
            )
        wavel = np.tile(wavel, reps)
        if waver is not None:
            if waver.shape[0] < pts_per_rep:
                waver = np.concatenate(
                    (waver, np.zeros(pts_per_rep - waver.shape[0])), axis=0
                )
            waver = np.tile(waver, reps)

        # different approaches to playing out the sound for different hardware configuration:

        if "pyaudio" in self.hardware:
            self.audio = pyaudio.PyAudio()
            chunk = 1024
            FORMAT = pyaudio.paFloat32
            # CHANNELS = 2
            CHANNELS=1
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
            ] = wavel  # order chosen so matches etymotic earphones on my macbookpro.
            postdur = int(float(postduration * self.in_sampleFreq))
            # rwave = read_array(len(wavel)+postdur, CHANNELS)
            write_array(self.stream, wave)
            self.stream.stop_stream()
            self.stream.close()
            self.audio.terminate()
            # self.ch1 = rwave[0::2]
            # self.ch2 = rwave[1::2]
            return

        if "PA5" in self.hardware:
            self.setAttens(atten_left=attns)

        if "RZ5D" in self.hardware:
            # swcount = -1
            self.present_stim(
                wavel, isi, reps, runmode, protocol
            )  # this sets up the NI card as well.
            #             deadmantimer = isi*(reps+1)+0.5  # just in case it doesn't stop as it should
            #             start_time = time.time()  # deadman start time
            # #            print('done? ', self.RZ5D.GetTargetVal(self.RZ5D_ParTags['SweepDone']))
            #             while self.RZ5D.GetTargetVal(self.RZ5D_ParTags['SweepDone']) == 0:  # wait for zSwDone to be set
            #                 cs = self.RZ5D.GetTargetVal(self.RZ5D_ParTags['CurrentSweep'])
            #                 if cs > swcount:
            #                     print('   Sweep = %d' % cs)
            #                     swcount = swcount + 1
            #                 time.sleep(0.1)
            #                 elapsed_time = time.time() - start_time  # elapsed time is in seconds
            #                 if elapsed_time > deadmantimer:
            #                     print('DeadmanExit')
            #                     break
            if runmode == "Preview":
                return
            else:
                self.RZ5D.setModeStr("Idle")  # was (RZ5D_Standby)
                time.sleep(2.0)  # added a couple of seconds here 20180306
            # self.task.stop()
            self.setAttens(atten_left=120)
            #    self.present_stim(wavel, waver)

        if "RP21" in self.hardware:
            # now take in some acquisition...
            a = self.RP21.ClearCOF()
            if a <= 0:
                print("pystim.playSound: Unable to clear RP2.1")
                return
            a = self.RP21.LoadCOFsf("C:\pyStartle\startle2.rco", self.samp_cof_flag)
            if a > 0 and self.debugFlag:
                print(
                    "pystim.playSound: Connected to TDT RP2.1 and startle2.rco is loaded"
                )
            else:
                print("pystim.playSound: Error loading startle2.rco?, error = %d" % (a))
                return
            self.trueFreq = self.RP21.GetSFreq()
            Ndata = np.ceil(0.5 * (stimulus_duration) * self.trueFreq)
            self.RP21.SetTagVal("REC_Size", Ndata)  # old version using serbuf  -- with
            # new version using SerialBuf, can't set data size - it is fixed.
            # however, old version could not read the data size tag value, so
            # could not determine when buffer was full/acquisition was done.

            if "PA5" in self.hardware:
                self.setAttens(
                    atten_left=attns[0], atten_right=attns[1]
                )  # set equal, but not at minimum...

            self.task.start()  # start the NI AO task

            a = self.RP21.Run()  # start the RP2.1 processor...
            a = self.RP21.SoftTrg(
                1
            )  # and trigger it. RP2.1 will in turn start the ni card

            while not self.task.isTaskDone():  # wait for AO to finish?
                self.RP21.Halt()
                if "NIDAQ" in self.hardware:
                    self.task.stop()
                return

            if "PA5" in self.hardware:
                self.setAttens()  # attenuators down (there is noise otherwise)
            # read the data...
            curindex1 = self.RP21.GetTagVal("Index1")
            curindex2 = self.RP21.GetTagVal("Index2")

            while (
                curindex1 < Ndata or curindex2 < Ndata
            ):  # wait for input data to be sampled
                self.RP21.Halt()
                return
                curindex1 = self.RP21.GetTagVal("Index1")
                curindex2 = self.RP21.GetTagVal("Index2")

            self.ch2 = self.RP21.ReadTagV("Data_out2", 0, Ndata)
            # ch2 = ch2 - mean(ch2[1:int(Ndata/20)]) # baseline: first 5% of trace
            self.ch1 = self.RP21.ReadTagV("Data_out1", 0, Ndata)
            self.RP21.Halt()

    def prepare_NIDAQ(self, wavel, waver=None):
        samplefreq = self.out_sampleFreq
        print("samplefreq :", samplefreq)
        wlen = len(wavel)
        #        daqwave = np.zeros(wlen)
        (wavel, clipl) = self.clip(wavel, 10.0)
        with nidaqmx.Task() as task:
            task.ao_channels.add_ao_voltage_chan(
                "/Dev1/ao0", min_val=-10.0, max_val=10.0, units=VoltageUnits.VOLTS
            )
            task.timing.cfg_samp_clk_timing(
                samplefreq,
                source="",
                sample_mode=AcquisitionType.FINITE,
                samps_per_chan=wlen,
            )
            # print(task.write(wavel,auto_start=True))
            # task.triggers.start_trigger.trig_type.DIGITAL_EDGE
            task.triggers.start_trigger.cfg_dig_edge_start_trig(
                trigger_source="/Dev1/PFI0", trigger_edge=Edge.RISING
            )
            # task.write(data,auto_start=True)
            task.write(wavel)
            task.start()

            task.wait_until_done(timeout=10)
            task.stop()
            print("loaded wavel")

            # task.start()

        # self.task = self.NIDevice.createTask()  # creat a task for the NI 6731 board.
        # self.task.CreateAOVoltageChan("/%s/ao0" % self.NIDevicename, "ao0", -10., 10.,
        # nidaq.Val_Volts, None)

    #  self.task.CreateAOVoltageChan("/%s/ao1" % self.NIDevicename, "ao1", -10., 10.,
    #                                nidaq.Val_Volts, None) # use 2 channels

    # (waver, clipr) = self.clip(waver, 10.0)

    #        daqwave[0:len(wavel)] = wavel
    # daqwave[len(wavel):] = waver # concatenate channels (using "groupbychannel" in writeanalogf64)
    # self.task.CfgSampClkTiming(None, samplefreq, nidaq.Val_Rising,
    #                            nidaq.Val_FiniteSamps, len(wavel))

    # # set up triggering for the NI (hardware trigger)
    # self.task.CfgDigEdgeStartTrig('PFI0',  nidaq.Val_Rising)
    # self.task.SetStartTrigType(nidaq.Val_DigEdge)

    # Need to set up count of trigger events, but not seeing it in NIDAQmx.h
    # It appears to require programming both clocks in the card:
    # see http://www.ni.com/tutorial/5382/en/
    # Counter 1 is set as the sample rate clock. It goes both to the Counter 2
    # and to the AO.
    # A "Samples Per Channel" value is set both in Counter 2, and AO
    # The Counter 2 output becomes the cock source for the AO task
    # Counter 2 is set up to be triggered from the hardware source
    # The Counter 2 task is configured for finite retriggerable pulse generation
    # The AO task is configured for continuous operation (but this is controlled by the
    # counters)

    # self.task.SetStartTrigRetriggerable(True) # so we set up a short waveform and let NI get retriggered by Rz5d
    # self.task.write(wavel)
    # self.task.start()

    def retrieveRP21_inputs(self):
        return (self.ch1, self.ch2)

    def HwOff(self):  # turn the hardware off.

        if "Soundcard" in self.hardware:
            try:
                self.stream.stop_stream()
                self.stream.close()
                self.audio.terminate()
            except:
                pass  # possible we never created teh stream...

        # if 'NIDAQ' in self.hardware:
        # if self.task is not None:
        #     self.task.SetStartTrigRetriggerable = 0
        #     self.task.stop()

        if "RP21" in self.hardware:
            self.RP21.Halt()

        if "RZ5D" in self.hardware:
            self.RZ5D_close()

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


"""
the following was taken from #http://hlzr.net/docs/pyaudio.html
it is used for reading and writing to the system audio device

"""


def write_array(stream, data):
    """
    Outputs a numpy array to the audio port, using PyAudio.
    """
    # Make Buffer
    buffer_size = struct.calcsize("@f") * len(data)
    output_buffer = ctypes.create_string_buffer(buffer_size)

    # Fill Up Buffer
    # struct needs @fffff, one f for each float
    dataformat = "@" + "f" * len(data)
    struct.pack_into(dataformat, output_buffer, 0, *data)

    # Shove contents of buffer out audio port
    stream.write(output_buffer)


def read_array(stream, size, channels=1):
    input_str_buffer = np.zeros((size, 1))  # stream.read(size)
    input_float_buffer = struct.unpack("@" + "f" * size * channels, input_str_buffer)
    return np.array(input_float_buffer)


if __name__ == "__main__":

    p = PyStim(hdw=["PA5", "NIDAQ", "RZ5D"], devicename="dev1")
    ni_sf = 100000
    w = np.cos(2 * np.pi * 2000.0 * np.arange(0, 0.2, 1.0 / ni_sf))
    p.setAttens(atten_left=30)
    p.present_stim(w)
    print("here we are")
    time.sleep(2.0)
    p.RZ5D.setModeStr("Idle")
    p.task.stop()
    p.setAttens(atten_left=120)
