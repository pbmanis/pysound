#!/usr/bin/env python


"""
pystim: a Python Class for interacting with hardware to produce sounds and
record signals.

Output hardware is either an National Instruments DAC card or a system sound card
If the NI DAC is available, TDT system 3 hardware is assumed as well for the
attenuators (PA5) and an RP2.1. or RZ5D

Hardware on the Rig 5 (ABR) system includes:
RP2.1
RZ5D
NI6731 (high speed 4 channel dac)
2 x PA5

If the system sound card is used, stimuli are generated. This is used only for testing.


12/17/2008-2024 Paul B. Manis, Ph.D.
UNC Chapel Hill
Department of Otolaryngology/Head and Neck Surgery
Supported by NIH Grants DC000425, DC004551 and DC015093 to PBM.

Refactored and modified version, includes access to rz5d to help synchronize acquisition.
August, 2017 and later
"""


"""
Old:
    (TDT manual system 3): 
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

New:
    (Synapse):
    Sweep cycling is controlled by PulseGen1 clock, which sets the
    interstimulus interval.
    The NIDAQ system is using a callback to reload the card at the end of every
    output, and is retriggered by the PulseGen1 (on digital out 0).
    The NIDAQ system can be turned off or on independently of the RZ5D state,
    and can be reloaded in between stimuli as well.

"""


import ctypes
from dataclasses import dataclass, field
import os
import pyqtgraph as pg
from pathlib import Path
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


def defemptylist():
    return []


@dataclass
class Stimulus_Status:
    """
    Create data structure for the status of the stimulus generator
    """

    controller: object = None
    running: bool = False
    stimulus_count: int = 0
    done: bool = False
    index: int = 0
    debugFlag: bool = False
    NI_devicename: str = ""
    NI_task: object = None  # dac output task
    NI_digitalin: object=None  # digital in task (for ttl pulse measurements)
    required_hardware: list = field(default_factory=defemptylist)
    hardware: list = field(default_factory=defemptylist)
    max_repetitions: int = 10


class Stimulus_Parameters:
    """
    Create data structure for the stimulus parameters
    """

    out_sampleFreq: float = 44100.0
    in_sampleFreq: float = 44100.0
    atten_left: float = 30.0
    atten_right: float = 120.0


class PyStim:
    def __init__(
        self, required_hardware=["NIDAQ"], ni_devicename="dev1", controller=None
    ):
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
        controller : object
            The parent class that provides the controls.
        """

        print("Required Hardware: ", required_hardware)
        self.State = Stimulus_Status()  # create instance of each data structure (class)
        self.State.required_hardware = required_hardware
        self.State.NI_devicename = ni_devicename
        self.State.controller = controller
        self.Stimulus = Stimulus_Parameters()
        self.enable_digital_display=False
        self.new_data = False
        self.repetition_number = 0
        self.online_plot = None
        self.find_hardware()
        #            device_info={"devicename": devicename}
        #        )  # population the self.State.hardware list
        self.TankName = []

    def find_hardware(self):
        """
        Find the hardware on the system.
        For non-windows systems, this just finds the system soundcard for testing
        Otherwise it looks for the requested hardware.
        Keeps track of available hardware in the self.State.hardware list

        Parameters
        ----------
        None

        """
        print("Operating system: ", opsys)
        print("nidaq_available: ", nidaq_available)
        if (
            opsys in ["Darwin", "Linux"] or nidaq_available is False
        ):  # If not on a Windows system, just set up soundcard
            self.setup_soundcard()
            self.State.hardware.append("Soundcard")
            self.Stimulus.out_samplefreq = 44100
        else:
            if "NIDAQ" in self.State.required_hardware and self.setup_nidaq():
                self.State.hardware.append("NIDAQ")
                self.setup_nidaq()
            if "RP21" in self.State.required_hardware and self.setup_RP21(
                "c:\\TDT\\OpenEx\\MyProjects\\Tetrode\\RCOCircuits\\tone_search.rcx"
            ):
                self.State.hardware.append("RP21")
            if "RZ5D" in self.State.required_hardware and self.setup_RZ5D():
                self.State.hardware.append("RZ5D")
            if "PA5" in self.State.required_hardware and self.setup_PA5(devnum=1):
                self.State.hardware.append("PA5")

        print("State.hardware: ", self.State.hardware)


    def setup_soundcard(self):
        if self.State.debugFlag:
            print(
                "pysounds.init: OS or available hardware only supports a standard sound card"
            )
        self.State.hardware.append("pyaudio")
        self.Stimulus.out_sampleFreq = 44100.0
        self.Stimulus.in_sampleFreq = 44100.0

    def setup_nidaq(self):
        # get the drivers and the activeX control (win32com)

        self.NIDevice = nidaqmx.system.System.local()
        self.NIDevicename = self.NIDevice.devices.device_names
        self.Stimulus.out_sampleFreq = 250000  # output frequency, in Hz
        return True

    def show_nidaq(self):
        """
        Report some information regardign the nidaq setup
        """

        print("pystim.show_nidaq  Devices:")
        print("devices: %s" % self.NIDevice.devices.device_names)
        # print ("devices: %s" % nidaq.NIDAQ.listDevices())
        print("getDevice: " % self.NIDevice)
        print(
            "\nAnalog Output Channels: %d"
            % self.NIDevice.devices[self.NIDevicename].ao_physical_chans.channel_names
        )
        # print ("\nAnalog Output Channels: %d" %  self.NIDevice.listAOChannels())

    def setup_PA5(self, devnum=1):
        """
        active X connection to attenuators

        Parameters
        ----------
        devnum : int (default = 1)
            The device number to connect to for the attenuator
        """
        self.PA5 = win32com.client.Dispatch("PA5.x")
        a = self.PA5.ConnectPA5("USB", devnum)
        if a > 0:
            if self.State.debugFlag:
                print("pystim.setup_PA5: Connected to PA5 Attenuator %d" % devnum)
        else:
            if "PA5" in self.State.required_hardware:
                raise IOError(f"PA5 requirement requested, but device  {devnum:d} not found")
            else:
                return False
        self.PA5.SetAtten(120.0)
        return True

    def setup_RP21(self, rcofile: str = ""):
        """
        active x connection to the RP2.1 Real-Time Processor

        Parameters
        ----------
        rcofile : str (default : '')
            The RCO file to connect to. Must be full path.
        """

        self.RP21_rcofile = rcofile
        self.RP21 = win32com.client.Dispatch("RPco.x")  # connect to RP2.1
        a = self.RP21.ConnectRP2("USB", 1)
        if a > 0 and self.State.debugFlag:
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
        a = self.RP21.LoadCOFsf(self.RP21_rcofile, self.samp_cof_flag)
        if a > 0:
            print(
                "pystim.setup_RP21: File %s loaded\n      and sample rate set to %f"
                % (self.RP21_rcofile, self.samp_fllist[self.camp_cof_flag])
            )
        else:
            print(
                "pystim.setup_RP21: Error loading RCO file %s, error = %d"
                % (rcofile, a)
            )
            return False
        self.Stimulus.out_sampleFreq = self.samp_flist[self.samp_cof_flag]
        self.Stimulus.in_sampleFreq = self.samp_flist[self.samp_cof_flag]
        return True

    def show_RP21(self):
        """
        TODO: maybe report RP2.1 info: cof rate, loaded circuit, sample freqs
        """
        pass

    def setup_RZ5D(self):
        self.RZ5D = tdt.SynapseAPI()
        if self.RZ5D.getModeStr() != "Idle":
            self.RZ5D.setModeStr("Idle")
        # print(dir(self.RZ5D))
        self.get_RZ5D_Params()
        # print(self.RZ5DParams["device_names"])
        # exit()
        return True


    def get_RZ5D_Params(self):
        self.RZ5DParams = {}  # keep a local copy of the parameters
        self.RZ5DParams["device_names"] = self.RZ5D.getGizmoNames()
        self.RZ5DParams["device status"] = self.RZ5D.getModeStr()

    def show_RZ5D(self):
        print("RZ5D Device Status: {0:d}".format(self.RZ5DParams["device_status"]))

    def get_RZ5D_Mode(self):
        return self.RZ5D.getModeStr()

    def RZ5D_close(self):
        if self.RZ5D.getModeStr() != "Idle":
            self.RZ5D.setModeStr("Idle")

    def getHardware(self):
        return (
            self.State.hardware,
            self.Stimulus.out_sampleFreq,
            self.Stimulus.in_sampleFreq,
        )

    # internal debug flag to control printing of intermediate messages
    def debugOn(self):
        self.State.debugFlag = True

    def debugOff(self):
        self.State.debugFlag = False

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
        if self.State.debugFlag:
            print("pystim.dbconvert: scale = %f for %f dB" % (sf, spl))
        return sf  # return a scale factor to multiply by a waveform normalized to 1

    def setAttens(self, atten_left=120.0, atten_right=120.0):
        if "PA5" in self.State.hardware:
            self.PA5.ConnectPA5("USB", 1)
            self.PA5.SetAtten(atten_left)
            if atten_right is not None:
                self.PA5.ConnectPA5("USB", 2)
                self.PA5.SetAtten(atten_right)

    def play_sound(
        self,
        wavel,
        waver=None,
        samplefreq:int=44100,  # default is for sound card
        postduration:float=0.005,
        attns:list=[20.0, 20.0],
        interstimulus_interval:float=1.0, # seconds
        repetitions:int=1,
        protocol:str="Search",
        storedata:bool=True,
    ):
        """
        play_sound sends the sound out to an audio device. In the absence of NI
        card or a usable TDT system, we (try to) use the system audio device (sound
        card, etc).
        The waveform is played in both channels on sound cards,
        possibly on both channels for other devices if there are 2 channels.

        Parameters
        ----------
        wavel : numpy array of floats
            Left channel waveform
        waver : numpy of floats
            Right channel waveform
        samplefreq : float
            output sample frequency (Hz)
        postduration : float (default: 0.05)
            Time after end of stimulus, in seconds, to set values to 0
        attns : 2x1 list (default: [20., 20.])
            Attenuator settings to use for this stimulus
        interstimulus_interval : float (default 1.0)
            InterStimulus interval: time from start to start of a repeated stimulus
        repetitions : int (default 1)
            Number of repetitions before returning.
        protocol: str (default "Search")
            protocol mode to use.
        storedata : bool (default: True)
            flag to force storage of data at end of run

        """
        if storedata:
            runmode = "Record"
        else:
            runmode = "Preview"
        # add end points (zeros) to the array for the specified time. This keeps the array
        # from having random final values.
        append_pts = int(postduration*samplefreq)
        wavel = np.append(wavel, np.zeros(append_pts))
        waver = np.append(waver, np.zeros(append_pts))
        # print("pystim hardware: ", self.State.hardware)
        self.stimstarttime = time.time()
        if "pyaudio" in self.State.hardware:
            print("pystim: Playing through pyaudio, system sound card")
            self.audio = pyaudio.PyAudio()
            chunk = 1024
            FORMAT = pyaudio.paFloat32
            # CHANNELS = 2
            CHANNELS = 1
            if self.State.debugFlag:
                print("pystim.play_sound: samplefreq: %f" % (RATE))
            self.stream = self.audio.open(
                format=FORMAT,
                channels=CHANNELS,
                rate=int(self.Stimulus.out_samplefreq),
                output=True,
                input=True,
                frames_per_buffer=chunk,
            )
            wave = np.zeros(2 * len(wavel))
            if len(wavel) != len(waver):
                print(
                    "pystim.play_sound: waves not matched in length: %d vs. %d (L,R)"
                    % (len(wavel), len(waver))
                )
                return
            (waver, clipr) = self.clip(waver, 20.0)
            (wavel, clipl) = self.clip(wavel, 20.0)
            wave[0::2] = waver
            wave[
                1::2
            ] = wavel  # order chosen so matches etymotic earphones on my macbookpro.
            postdur = int(float(postduration * self.Stimulus.in_sampleFreq))

            write_array(self.stream, wave)
            self.stream.stop_stream()
            self.stream.close()
            self.audio.terminate()
            return

        if "PA5" in self.State.hardware:
            # print("setting PA5")``
            self.setAttens(atten_left=attns, atten_right=attns)

        if "RZ5D" in self.State.hardware:
            # print("Accessing RZ5D")
            swcount = -1
            timeout = interstimulus_interval * repetitions + 1
            # Start and run the stim/recording for specified # sweeps/time.
            # self.RZ5D.setModeStr(runmode)
            self._present_stim(
                waveforms=wavel,
                interstimulus_interval=interstimulus_interval,
                repetitions=repetitions,
                runmode=runmode,
                protocol=protocol,
                timeout=timeout,
            )  # this sets up the NI card.
           

    def _present_stim(
        self,
        waveforms,
        interstimulus_interval: float = 1.0,
        repetitions: int = 1,
        runmode: str = "Preview",
        protocol: str = "Search",
        timeout: float = 10.0,
    ):
        """
        Set up and enable the stimulus presentation.
        We:
            1. start the RZ5D
            2. set up the pulse generator
            3. load and arm the NI card
            Everything after that is controlled by synapse and hardware
            triggers, so we just return.
        """
        print("   present_stim: RZ5D")

        if (
            self.RZ5D.getModeStr() != runmode
        ):  # make sure the rz5d is in the requested mode first
            self.RZ5D.setModeStr(runmode)

        ##################################################################################
        # Set up the stimulus timing
        # We use the PulseGen0 to write to digital line out 0
        # This bit controls/triggers the timing of the stimuli (interstimulus interval)
        # by trigginering the NI card.
        pgen = "PulseGen1"
        # print("   rz5d interstimulus interval: ", interstimulus_interval)
        params = self.RZ5D.getParameterNames(pgen)
        # print(f"{pgen:s} params: {params!s}")
        print("PulseGen1:pulse period ", self.RZ5D.getParameterValue(pgen, "PulsePeriod"))
        self.RZ5D.setParameterValue(pgen, "PulsePeriod", interstimulus_interval)
        self.RZ5D.setParameterValue(pgen, "DutyCycle", 2.0)  # 1 msec pulse
        self.RZ5D.setParameterValue(pgen, "Enable", 1)
        # print("after: ", self.RZ5D.getParameterValue(pgen, "PulsePeriod"))
        ##################################################################################

        for nr in range(repetitions):
            self.State.done = False
            # load up NIDAQ to go. This takes about 50 msec depending on the waveform
            # ni_timer = time.time()
            self.prepare_NIDAQ(waveforms, repetitions=repetitions)
            self.repetition_number = nr


    def stop_nidaq(self):
        """
        Only stop the DAC, not the RZ5D
        This is used when reloading a new stimulus.
        """
        self.State.running = False
        # if self.State.NI_task is not None:
            # self.State.NI_task.wait_until_done(timeout=2.0)
            # self.State.NI_task.close()  # release resources
            # self.State.NI_task = None  # need to destroy value
            # self.State.running = False
        # if self.State.NI_digitalin is not None:
        #     self.State.NI_digitalin.close()
        #     self.State.NI_digitalin = None


    def stop_recording(self):
        """
        Stop the entire system (DAC and RZ5D)
        """

        self.RZ5D.setModeStr("Idle")
        i = 0
        while(self.RZ5D.getModeStr() != "Idle"):
            time.sleep(0.2)
            i += 1
            if i > 25:
                break
        self.stop_nidaq()
        self.setAttens(atten_left=120)

    def arm_NIDAQ(self):
        """
        Load up the NI card output buffer, and set the triggers
        This gets the card ready to put out the buffer with the
        next trigger pulse
        """
        # print("    Arm: ", time.time()-self.stimstarttime)
        self.State.NI_task.write(self.waveout, auto_start=False)
        #  self.State.NI_task.triggers.start_trigger.trig_type.DIGITAL_EDGE
        self.State.NI_task.triggers.start_trigger.cfg_dig_edge_start_trig(
            trigger_source="/Dev1/PFI0",
            trigger_edge=Edge.RISING,
        )
        # print("daq armed")
        self.State.running = True

    def re_arm_NIDAQ(self, task_handle, status, callback_data):
        """
        Callback for when the daq is done...
        Re arm the dac card and start the task again
        """
        print("    re-arm: ", time.time() - self.stimstarttime)
        if status != 0:
            self.stop_recording()  # nidaq failure?
            return False

        if self.State.NI_task.is_task_done():
            self.State.NI_task.stop()
            self.arm_NIDAQ()  # reload and re-arm the trigger
            self.State.NI_task.start()
            self.State.stimulus_count += 1

            # counter_elapsed = self.State.stimulus_count > self.repetitions
            #   controller_running = self.State.controller.running
            timeout = False  # (time.time() - self.start_time) > self.timeout
            # if counter_elapsed or (not self.State.running) or timeout:
            if (not self.State.running) or timeout:
                self.stop_nidaq()
                self.State.done = True
                return False
        return True

    def load_and_arm_NIDAQ(self):
        """
        Setup of NI card for AO.
        Creates a task for the card, sets parameters, clock rate, and loads the waveform.
        Wrapping in a with statement tears down the tasks when we are done.
        This does not block the GUI.
        """
        # print("    Load and run: ")
        # print("Stimulus duration: ", len(self.waveout)/self.Stimulus.out_sampleFreq)
        this_starttime = time.time()
        failed = False
        self.new_data = False # flag to alert us to new data.
        self.digital_data = None
        with nidaqmx.task.Task("NI_DAC_out") as self.State.NI_task, nidaqmx.task.Task("NI_DI_In") as self.State.NI_digitalin:
            channel_name = f"/{self.State.NI_devicename:s}/ao0"
            self.State.NI_task.ao_channels.add_ao_voltage_chan(  # can only do this once...
                channel_name, min_val=-10.0, max_val=10.0, units=VoltageUnits.VOLTS
            )
            # self.State.NI_task.register_done_event(self.re_arm_NIDAQ)

            # print("Load and arm: SFreq = ", self.Stimulus.out_sampleFreq, "Wave len: ", len(self.waveout))
            # print("     stim dur: ", len(self.waveout)/self.Stimulus.out_sampleFreq)
            self.State.NI_task.timing.cfg_samp_clk_timing(
                self.Stimulus.out_sampleFreq,
                source="",
                sample_mode=AcquisitionType.FINITE,
                samps_per_chan=len(self.waveout),
            )
            self.State.NI_task.write(self.waveout, auto_start=False)

            self.State.NI_task.triggers.start_trigger.cfg_dig_edge_start_trig(
                    trigger_source="/Dev1/PFI0",
                    trigger_edge=Edge.RISING,
                )
            
            for line in [0,1,2,3]:
                self.State.NI_digitalin.di_channels.add_di_chan(
                f"/Dev1/port0/line{line:d}",
                )

            self.State.NI_digitalin.timing.cfg_samp_clk_timing(
                self.Stimulus.out_sampleFreq,
                source="/Dev1/ao/SampleClock",
                # active_edge=Edge.RISING,
                sample_mode=AcquisitionType.FINITE,
                samps_per_chan=len(self.waveout)
            )

            # if not self.State.running:
            #     self.State.NI_task.stop()
            #     return False
            # self.arm_NIDAQ()
        #  self.State.NI_task.triggers.start_trigger.trig_type.DIGITAL_EDGE

        # print("daq armed")
            self.State.NI_digitalin.start()
            self.State.NI_task.start()  # and start it
            # print("    restarting at ", time.time()-self.start_time, "secs since last stim")
            # print("    ", time.time() - self.stimstarttime, "secs since start of sequence")
            while not self.State.NI_task.is_task_done():
                now_time = time.time()
                if now_time - this_starttime > 5.0:
                    failed = True
                    raise ValueError("arming nidaq/task execution FAILED")

            self.digital_data = np.array(self.State.NI_digitalin.read(len(self.waveout)))
            print("data shape: ", self.digital_data.shape)
            self.new_data = True
        self.State.running = False
        print("prepping plot", self.enable_digital_display)
        if self.new_data and self.enable_digital_display:
            print("ok to plot")
            t = np.linspace(0, len(self.waveout)/self.Stimulus.out_sampleFreq, len(self.waveout))
            n_erase = 10
            thresh = 0.5
            # if (self.repetition_number % n_erase )== 0:
            #     print("clearing", n_erase, self.repetition_number)
            #     self.online_plot.clear()
            print("data shape: ", self.digital_data.shape)
            tr_x = []
            tr_y = []

            for i in range(self.digital_data.shape[0]):

                tcross = np.diff(self.digital_data[i,:] > thresh, prepend=False)
                tindex = np.argwhere(tcross)[::2,0]
                # print("tindex: ", tindex)
                # print("outfreq: ", self.Stimulus.out_sampleFreq)
                tindex = tindex/self.Stimulus.out_sampleFreq
                # if len(tindex) == 1:
                #     continue

                tv = np.ones_like(tindex)+(i % n_erase)/n_erase
                # mpl.plot(t, n + data[i,:])
                # print(tv, tindex)
                # print(len(tindex), len(tv), np.max(tindex), np.min(tindex))
                tr_x.extend(tindex)
                tr_y.extend(tv)
                self.online_plot.scatterPlot(tindex, tv, 
                                                 symbolBrush=pg.mkBrush("r"), # (pg.intColor(i, hues=n_erase)), 
                                                 symbol='o', symbolSize=3)
            self.online_plot.setXRange(0, 0.5)
            self.online_plot.setYRange(0, 5)
            # pg.QtGui.QGuiApplication.processEvents()
            self.repetition_number += 1

            self.new_data = False
        if not failed:
            return True
        else:
            return False

    def prepare_NIDAQ(
        self, wavel, waver=None, repetitions: int = 1, timeout: float = 1200.0
    ):
        """
        Set up and initialize the NIDAQ card for output,
        then let it run and keep up with each task completion
        so it can be retriggered on the next trigger pulse.
        Configured so that if we are currently running, the run is immediately stopped
        so we can setup right away.
        """
        # print("\nPrepare NIDAQ")
        # self.stop_nidaq()  # stop the DAC if it is running
        # update the waveform and rep counter
        self.waveout = wavel
        self.repetitions = repetitions
        self.State.stimulus_count = 0
        (self.waveout, clipl) = self.clip(
            self.waveout, 10.0
        )  # clip the wave if it's >10V
        self.start_time = time.time()
        self.timeout = timeout
        self.load_and_arm_NIDAQ()

    def retrieveRP21_inputs(self):
        return (self.ch1, self.ch2)

    def HwOff(self):  # turn the hardware off.

        if "Soundcard" in self.State.hardware:
            try:
                self.stream.stop_stream()
                self.stream.close()
                self.audio.terminate()
            except:
                pass  # possible we never created teh stream...

        if "NIDAQ" in self.State.hardware:
            self.stop_nidaq()

        if "RP21" in self.State.hardware:
            self.RP21.Halt()

        if "RZ5D" in self.State.hardware:
            self.RZ5D_close()

    # clip data to max value (+/-) to avoid problems with daqs
    def clip(self, data, maxval):
        # t0 = time.time()
        if self.State.debugFlag:
            print(
                "pystim.clip: max(data) = %f, %f and maxval = %f"
                % (max(data), min(data), maxval)
            )
        clip = 0
        u = np.where(data >= maxval)
        ul = list(np.transpose(u).flat)
        if len(ul) > 0:
            data[ul] = maxval
            clip = 1  # set a flag in case we want to know
            if self.State.debugFlag:
                print("pystim.clip: clipping %d positive points" % (len(ul)))
        minval = -maxval
        v = np.where(data <= minval)
        vl = list(np.transpose(v).flat)
        if len(vl) > 0:
            data[vl] = minval
            clip = 1
            if self.State.debugFlag:
                print("pystim.clip: clipping %d negative points" % (len(vl)))
        if self.State.debugFlag:
            print(
                "pystims.clip: clipped max(data) = %f, %f and maxval = %f"
                % (np.max(data), np.min(data), maxval)
            )
        # print("clipping took: ", time.time() - t0)
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
    ni_sampld_frequency = 100000
    w = np.cos(2 * np.pi * 2000.0 * np.arange(0, 0.2, 1.0 / ni_sample_frequency))
    p.setAttens(atten_left=30)
    p._present_stim(w)
    time.sleep(2.0)
    p.RZ5D.setModeStr("Idle")
    p.task.stop()
    p.setAttens(atten_left=120)
