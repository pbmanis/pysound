
import sys, time, numpy
import nidaqmx    
from nidaqmx.constants import AcquisitionType, Edge, VoltageUnits


NIsys=nidaqmx.system.System.local()    
for device in NIsys.devices:
  print(device)
  NIDevice=device
print("  ", NIDevice)
print(NIDevice.ao_trig_usage[1])
NIDevice.ao_trig_usage[1]
print(dir(NIDevice))
print('DAQ connected, info above')

## Output task

def outputTest():
  
  with nidaqmx.Task() as task:
    task.ao_channels.add_ao_voltage_chan('/Dev1/ao0',min_val=-10.,max_val=10.,units=VoltageUnits.VOLTS)
    task.timing.samp_clk_rate=200000.
    task.timing.samp_timing_type.HANDSHAKE
    print('sample rate: ',task.timing.samp_clk_rate)
    data = numpy.zeros((1000,), dtype=numpy.float64)
    data[20:40] = 5.0
    data[60:80] = 5.0

    task.triggers.start_trigger.trig_type.DIGITAL_EDGE
    task.timing.cfg_samp_clk_timing(200000,source=u'',sample_mode=AcquisitionType.FINITE,samps_per_chan=1000)

    task.triggers.start_trigger.trig_type.DIGITAL_EDGE
    task.triggers.start_trigger.cfg_dig_edge_start_trig(trigger_source=u'/Dev1/PFI0',trigger_edge=Edge.RISING)
    # task.write(data,auto_start=True)
    task.write(data)
    task.start()

    task.wait_until_done(timeout=1)
    task.stop()


outputTest()
