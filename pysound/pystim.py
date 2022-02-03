#!/usr/bin/env python


"""
pystim: a Python Class for interacting with hardware to produce sounds and
record signals.

Output hardware is either an National Instruments DAC card or a system sound card
If the NI DAC is available, TDT system 3 hardware is assumed as well for the
attenuators (PA5) and an RP2.1.

Hardware on the Rig 5 (ABR) system includes:
RP2.1
RZ5D
NI6731 (high speed 4 channel dac)
2 x PA5

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
from dataclasses import dataclass, field
import os
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
    controller : object = None
    running : bool = False
    index: int = 0
    debugFlag: bool = False
    NI_devicename: str = ''
    NI_task: object = None
    required_hardware: list = field(default_factory=defemptylist)
    hardware : list = field(default_factory=defemptylist)
    max_repetitions: int = 10

class Stimulus_Parameters:
    """
    Create data structure for the stimulus parameters
    """
    out_sampleFreq: float = 44100.
    in_sampleFreq: float = 44100.
    atten_left: float = 30.0
    atten_right: float = 120.0


class PyStim:
    def __init__(self, required_hardware=["Soundcard"], ni_devicename="dev1", controller=None):
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
        
        self.State = Stimulus_Status()  # create instance of each data structure (class)
        self.State.required_hardware = required_hardware
        self.State.NI_devicename = ni_devicename
        self.State.controller = controller
        self.Stimulus = Stimulus_Parameters() 
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
            if "PA5" in self.State.required_hardware and self.setup_PA5():
                self.State.hardware.append("PA5")
            if "RZ5D" in self.State.required_hardware and self.setup_RZ5D():
                self.State.hardware.append("RZ5D")

    def setup_soundcard(self):
        if self.State.debugFlag:
            print("pysounds.init: OS or available hardware only supports a standard sound card")
        self.State.hardware.append("pyaudio")
        self.Stimulus.out_sampleFreq = 44100.0
        self.Stimulus.in_sampleFreq = 44100.0

    def setup_nidaq(self):
        # get the drivers and the activeX control (win32com)

        self.NIDevice = nidaqmx.system.System.local()
        self.NIDevicename = self.NIDevice.devices.device_names
        self.Stimulus.out_sampleFreq = 200000  # output frequency, in Hz
        return True

    def show_nidaq(self):
        """
        Report some information regardign the nidaq setup
        """
        
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
                print("pysounds.init: Connected to PA5 Attenuator %d" % devnum)
        else:
            if "PA5" in self.State.required_hardware:
                raise IOError("PA5 requirement requested, but device not found")
            else:
                return False
        self.PA5.SetAtten(120.0)
        return True

    def setup_RP21(self, rcofile:str = ''):
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
            print("pysounds.init: Error loading RCO file %s, error = %d" % (rcofile, a))
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
        return True

    def getParams_RZ5D(self):
        self.RZ5DParams = {}  # keep a local copy of the parameters
        self.RZ5DParams["device_name"] = self.RZ5D.getGizmoNames()
        self.RZ5DParams["device status"] = self.RZ5D.getModeStr()

    def show_RZ5D(self):
        print("Device Status: {0:d}".format(self.RZ5DParams["device_status"]))

    def _present_stim(
        self,
        waveforms,
        stimulus_period: float = 1.0,
        repetitions: int = 1,
        runmode: str = "Record",
        protocol: str = "Search",
        timeout: float = 10.0,
    ):
        ##################################################################################
        # We use the PulseGen1 to write to digital line out 0
        # This bit controls/triggers the timing of the stimuli (interstimulus interval)

        params = self.RZ5D.getParameterNames('PulseGen1')
        self.RZ5D.setParameterValue('PulseGen1', 'PulsePeriod', stimulus_period)
        self.RZ5D.setParameterValue('PulseGen1', 'DutyCycle', 1.0) # 1 msec pulse
        self.RZ5D.setParameterValue('PulseGen1', 'Enable', 1.0)
        # for param in params:
        #     info = self.RZ5D.getParameterInfo('PulseGen1', param)
        #     print(f"Param: {param:s}, Info: {str(info):s}")
        ##################################################################################


        self.prepare_NIDAQ(waveforms, repetitions=repetitions)  # load up NIDAQ to go

        if self.RZ5D.getModeStr() != runmode:
            # TFR 20101007 removing the block checking/matching/storing stuff- important! 
            # This needs to be revised for record mode!
            # if runmode == "Record":
            #     protocol.replace(" ", "")
            #     subject = self.RZ5D.getCurrentSubject()
            #     self.TankName = self.RZ5D.getCurrentTank()
            #     # tankname=self.RZ5D.getCurrentTank()
            #     newblock = subject + protocol + "{:03d}".format(self.State.index)
            #     for checkBlocks in self.RZ5D.getKnownBlocks():
            #         if os.path.join(self.TankName, newblock) in checkBlocks:
            #             self.State.index = self.State.index + 1
            #             newblock = subject + protocol + "{:03d}".format(self.State.index)
            #     self.RZ5D.setCurrentBlock(newblock)
            #     self.State.index = self.State.index + 1
            self.RZ5D.setModeStr(runmode)

    def RZ5D_close(self):
        if self.RZ5D.getModeStr() != "Idle":
            self.RZ5D.setModeStr("Idle")
        # time.sleep(1.0)

    def getHardware(self):
        return (self.State.hardware, self.Stimulus.out_sampleFreq, self.Stimulus.in_sampleFreq)

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
        In the absence of NI card, and TDT system, we (try to) use the system audio device (sound card, etc)
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
        isi : float (default 1.0)
            Interstimulus interval
        reps : int (default 1)
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
        if "pyaudio" in self.State.hardware:
            self.audio = pyaudio.PyAudio()
            chunk = 1024
            FORMAT = pyaudio.paFloat32
            # CHANNELS = 2
            CHANNELS = 1
            if self.State.debugFlag:
                print("pysounds.play_sound: samplefreq: %f" % (RATE))
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
            postdur = int(float(postduration * self.Stimulus.in_sampleFreq))

            write_array(self.stream, wave)
            self.stream.stop_stream()
            self.stream.close()
            self.audio.terminate()
            return

        if "PA5" in self.State.hardware:
            self.setAttens(atten_left=attns, atten_right=attns)

        if "RZ5D" in self.State.hardware:
            swcount = -1
            timeout = isi * reps + 1
             # Start and run the stim/recording for specified # sweeps/time.
           # self.RZ5D.setModeStr(runmode)
            self._present_stim(
                wavel,
                stimulus_period=isi,
                repetitions=reps,
                runmode=runmode,
                protocol=protocol,
                timeout=timeout,
            )  # this sets up the NI card as well.
 
            
            # while time.time()-start_time < deadmantimer:
            #     time.sleep(0.01)
            # sweeps_start = self.RZ5D.getSystemStatus()['recordSecs']
            # currTime = 0
            # prevTime = 0
            # print(f"ISI: {isi:.3f}  reps: {reps:d}, runmode: {runmode:s}")
            # print(f"Running for maximum of: {deadmantimer:.2f} seconds")
            # while currTime < deadmantimer:
            #     currTime = self.RZ5D.getSystemStatus()['recordSecs']-sweeps_start
            #     if prevTime != currTime:
            #         print(f"Running, sweeps time elapsed: {currTime:d} sec")
            #     prevTime = currTime
            # TFR end comment 10/12/21
            # print("nidaq has stopped")
            # if runmode == "Preview":
            #     return
            # else:
            #     self.RZ5D.setModeStr("Idle")  # was (RZ5D_Standby)

            # self.setAttens(atten_left=120)

    def stop_stim(self):
        if self.State.NI_task is not None:
            self.State.NI_task.close()  # release resources
            self.State.NI_task = None  # need to destroy value
        self.RZ5D.setModeStr("Idle")
        self.setAttens(atten_left=120)

    def arm_NIDAQ(self):
        """
        Load up the NI card output buffer, and set the triggers
        This gets the card ready to put out the buffer with the
        next trigger pulse
        """
        self.State.NI_task.write(self.waveout, auto_start=False)
      #  self.State.NI_task.triggers.start_trigger.trig_type.DIGITAL_EDGE
        self.State.NI_task.triggers.start_trigger.cfg_dig_edge_start_trig(
            trigger_source="/Dev1/PFI0", trigger_edge=Edge.RISING,
        )

    def re_arm_NIDAQ(self, task_handle, status, callback_data):
        """
        Callback for when the daq is done... 
        Re arm the dac card and start the task again
        """

        if status != 0:
            self.stop_stim()
            return False

        if self.State.NI_task.is_task_done():
            self.State.NI_task.stop()
            self.arm_NIDAQ()  # reload
            self.State.NI_task.start()
            self.stim_counter += 1
     
            counter_elapsed = self.stim_counter > self.repetitions
            controller_running = self.State.controller.running
            timeout = False # (time.time() - self.start_time) > self.timeout
            if (
                counter_elapsed
                or (not controller_running)
                or timeout
            ):
                # print(counter_elapsed, controller_running, timeout)
                # print("Stopping NI task... (in re_arm_NIDAQ)")
                self.stop_stim()
                return False
        return True

    def load_and_arm_NIDAQ(self):
        """
        Initial setup of NI card for AO. 
        Creates a task for the card, sets parameters, clock rate,
        and does setup if needed. 
        A callback is registered so that when the task is done, the
        board is re-armed for the next trigger. 
        This does not block the GUI.
        """
        self.State.NI_task = nidaqmx.task.Task("NI_DAC_out")
        channel_name = f"/{self.State.NI_devicename:s}/ao0"
        self.State.NI_task.ao_channels.add_ao_voltage_chan(  # can only do this once...
            channel_name, min_val=-10.0, max_val=10.0, units=VoltageUnits.VOLTS
        )
        self.State.NI_task.register_done_event(self.re_arm_NIDAQ)

        self.State.NI_task.timing.cfg_samp_clk_timing(
            self.Stimulus.out_sampleFreq,
            source="",
            sample_mode=AcquisitionType.FINITE,
            samps_per_chan=len(self.waveout),
        )

        if not self.State.controller.running:
            self.State.NI_task.stop()
            return False
        self.arm_NIDAQ()   # setup the DAC card
        self.State.NI_task.start()  # and start it
        return True

    def prepare_NIDAQ(
        self, wavel, waver=None, repetitions: int = 1, timeout: float = 1200.0
    ):
        """
        Set up and initialize the NIDAQ card for output,
        then let it run and keep up with each task completion
        so it can be retriggered on the next trigger pulse.
        """

        self.stop_flag = False
        self.State.NI_task = None
        self.waveout = wavel
        self.repetitions = repetitions
        (self.waveout, clipl) = self.clip(
            self.waveout, 10.0
        )  # clip the wave if it's >10V
        self.stim_counter = 0
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
            if self.State.NI_task is not None:
                self.State.NI_task.stop()

        if "RP21" in self.State.hardware:
            self.RP21.Halt()

        if "RZ5D" in self.State.hardware:
            self.RZ5D_close()

    # clip data to max value (+/-) to avoid problems with daqs
    def clip(self, data, maxval):
        if self.State.debugFlag:
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
            if self.State.debugFlag:
                print("pysounds.clip: clipping %d positive points" % (len(ul)))
        minval = -maxval
        v = np.where(data <= minval)
        vl = list(np.transpose(v).flat)
        if len(vl) > 0:
            data[vl] = minval
            clip = 1
            if self.State.debugFlag:
                print("pysounds.clip: clipping %d negative points" % (len(vl)))
        if self.State.debugFlag:
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
    ni_sampld_frequency = 100000
    w = np.cos(2 * np.pi * 2000.0 * np.arange(0, 0.2, 1.0 / ni_sample_frequency))
    p.setAttens(atten_left=30)
    p._present_stim(w)
    time.sleep(2.0)
    p.RZ5D.setModeStr("Idle")
    p.task.stop()
    p.setAttens(atten_left=120)
