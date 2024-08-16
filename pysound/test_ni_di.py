"""
test reading timing of digital transistions on NI6731

aftger:
https://nspyre.readthedocs.io/en/latest/guides/ni-daqmx.html

"""

import nidaqmx
import numpy as np
import matplotlib.pyplot as mpl

from nidaqmx.constants import (
    AcquisitionType,
    CountDirection,
    Edge,
    READ_ALL_AVAILABLE,
    TaskMode,
    TriggerType,
    VoltageUnits,
)
from nidaqmx.stream_readers import CounterReader


required_hardware = ["NIDAQ"]
ni_devicename = "Dev1"

NI_Device = nidaqmx.system.System.local()
NI_Devicename = NI_Device.devices.device_names
print(NI_Devicename[0])
sampling_rate = 1000  # input frequency, in Hz (10 usec resolution)

# Let's load up the NI-DAQmx system that is visible in the
# Measurement & Automation Explorer (MAX) software of NI-DAQmx for
# the local machine.
system = nidaqmx.system.System.local()
# We know on our current system that our DAQ is named 'Dev1'
DAQ_device = system.devices["Dev1"]
# create a list of all the counters available on 'DAQ1'
counter_names = [ci.name for ci in DAQ_device.ci_physical_chans]
print("Counter names: ", counter_names)
# note that using the counter output channels instead of the inputs
# includes the '[device]/freqout' output, which is not a counter
print(
    "counter out Physical channels: ", [co.name for co in DAQ_device.co_physical_chans]
)
print(
    "analog output Physical channels: ",
    [co.name for co in DAQ_device.ao_physical_chans],
)
print("digital input Physical channels: ", [co.name for co in DAQ_device.di_ports])
print("digital input Physical lines: ", [co.name for co in DAQ_device.di_lines])
print(
    "digital trig_usage, trig_supported, max rate: ",
    DAQ_device.di_trig_usage,
    DAQ_device.dig_trig_supported,
    DAQ_device.di_max_rate,
)

print("Physical AO channels: ", [ao.name for ao in DAQ_device.ao_physical_chans])
# print(dir(DAQ_device))
duration = 1.0
sample_rate = 100000
npoints = int(duration*sample_rate)
print(npoints)
t = np.linspace(0, npoints/sample_rate, npoints)
wave = np.sin(2*np.pi*1000*t*sample_rate)

f, ax = mpl.subplots(1,1)

def plot_data(data, n, ax):
    t = np.linspace(0, npoints/sample_rate, npoints)
    cols = ['r', 'g', 'b', 'y']
    for i in range(data.shape[0]):
        thresh = 0.5
        tcross = np.diff(data[i,:] > thresh, prepend=False)
        tindex = np.argwhere(tcross)[::2,0]
        tindex = tindex/sample_rate
        tv = 0.5*np.ones_like(tindex)+0.1*i + n
        # mpl.plot(t, n + data[i,:])
        ax.plot(tindex, tv, 'o', color=cols[i], markersize=2)
    mpl.draw()

for i in range(10):
    print("Rep: ", i)
    with nidaqmx.Task("test_di") as di_task, nidaqmx.Task("ao clock") as ao_task:
        ao_task.ao_channels.add_ao_voltage_chan(  # can only do this once...
            "Dev1/ao0", min_val=-10.0, max_val=10.0, units=VoltageUnits.VOLTS,
        )
        ao_task.timing.cfg_samp_clk_timing(
            sample_rate, source="", sample_mode=AcquisitionType.FINITE,
                samps_per_chan=npoints
        )
        ao_task.triggers.start_trigger.cfg_dig_edge_start_trig(
            trigger_source="/Dev1/PFI0",
            trigger_edge=Edge.RISING,
        )
        ao_task.write(wave, auto_start=False)
        ao_task.control(TaskMode.TASK_COMMIT)

        for line in [0,1,2,3]:
            di_task.di_channels.add_di_chan(
            f"/Dev1/port0/line{line:d}",
            # edge=Edge.RISING,
            )

        di_task.timing.cfg_samp_clk_timing(
            sampling_rate,
            source="/Dev1/ao/SampleClock",
            # active_edge=Edge.RISING,
            sample_mode=AcquisitionType.FINITE,
            samps_per_chan=npoints
        )

        di_task.start()
        ao_task.start()

        while not ao_task.is_task_done():
            pass
        data = np.array(di_task.read(npoints))
    plot_data(data, n=i, ax=ax)

mpl.show()

    