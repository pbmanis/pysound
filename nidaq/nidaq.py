from ctypes import *
import sys
from . import cheader
import numpy
import re
import types
import ctypes
import os


int8 = c_byte
uInt8 = c_ubyte
int16 = c_short
uInt16 = c_ushort
int32 = c_long
uInt32 = c_ulong
float32 = c_float
float64 = c_double
int64 = c_longlong
uInt64 = c_ulonglong
bool32 = uInt32
TaskHandle = uInt32
DAQmxEveryNSamplesEventCallbackPtr = CFUNCTYPE(int32, c_ulong, c_long, c_ulong, c_void_p)
DAQmxDoneEventCallbackPtr = CFUNCTYPE(int32, c_ulong, c_long, c_void_p)
DAQmxSignalEventCallbackPtr = CFUNCTYPE(int32, c_ulong, c_long, c_void_p)

def init():
  ## System-specific code
  headerFiles = [os.path.join(os.path.dirname(__file__), "NIDAQmx.h")]
  xmlFiles = [os.path.join(os.path.dirname(__file__), "NIDAQmx.xml")]
  defs = cheader.getDefs(headerFiles)
  global NIDAQ
  NIDAQ = _NIDAQ()
  for k in defs:
    if k is not None:
      setattr(sys.modules[__name__], re.sub('^DAQmx_?', '', k), defs[k])
  NIDAQ.functions = cheader.getFuncs(xmlFiles)


class NIDAQError(Exception):
  pass
class NIDAQWarning(Exception):
  pass

class _NIDAQ:
  NIDAQ_CREATED = False
  def __init__(self):
    if _NIDAQ.NIDAQ_CREATED:
      raise Exception("Will not create another nidaq instance--use the pre-existing NIDAQ object.")
    self.nidaq = windll.nicaiu
    self.devices = {}
    # :TODO: initialize the driver
    _NIDAQ.NIDAQ_CREATED = True
  #  print('Created NIDAQ instance _NIDAQ')

  def __repr__(self):
    return "<niDAQmx driver wrapper>"

  def listDevices(self):
    return self.GetSysDevNames().split(b", ")

  def getDevice(self, dev):
    return Device(dev, self)
  
  def __getattr__(self, attr):
    if attr[0] != "_" and hasattr(self.nidaq, 'DAQmx' + attr):
      return lambda *args: self.call(attr, *args)
    else:
      raise NameError(attr)

  # import codecs
  # def makebytes(x):
  #     return codecs.latin_1_encode(x)[0]

  def call(self, func, *args):
    func = 'DAQmx' + func
    ret = None
    retType, argSig = self.functions[func]
    #print "CALL: ", func, args, argSig
    
    if func[:8] == "DAQmxGet":
      if argSig[-1][0] == 'data':
        ret = getattr(ctypes, argSig[-1][1])()
        args += (byref(ret),)
      elif argSig[-2][1:] == ('c_char', 1) and argSig[-1][1:] in [('c_ulong', 0), ('c_long', 0)]:
        #print "correct for buffer return"
        tmpargs = args + (getattr(ctypes, argSig[-2][1])(), getattr(ctypes, argSig[-1][1])())
        buffSize = self._call(func, *tmpargs)
        ret = create_string_buffer(b'\0' * buffSize)
        args += (ret, buffSize)
    
    cArgs = []
    for i in range(0, len(args)):
      arg = args[i]
      # print('ARG: ',arg)
      if type(args[i]) in [float, int, int, bool] and argSig[i][2] == 0:
        #print func, i, argSig[i][0], argSig[i][1], type(arg)
        arg = getattr(ctypes, argSig[i][1])(arg)
      # if type(args[i]) in [str]:
      #   print('arg: ', str(args[i]))
      #   s = bytearray(arg, 'utf-8')
      #   print('   ', s)
      #   arg = getattr(ctypes, argSig[i][1])(makebytes(args[i]))
      cArgs.append(arg)
    
    #print "  FINAL CALL: ", cArgs
    errCode = self._call(func, *cArgs)
    if errCode < 0:
      raise NIDAQError(errCode, "Function '%s%s'" % (func, str(args)), *self.error(errCode))
    elif errCode > 0:
      raise NIDAQWarning(errCode, "Function '%s%s'" % (func, str(args)), *self.error(errCode))
    
    if ret is None:
      return True
    else:
      return ret.value
    
  def _call(self, func, *args):
    try:
      return getattr(self.nidaq, func)(*args)
    except:
      print((func, args))
      raise
    
  def error(self, errCode):
    return (self.GetErrorString(errCode),
           self.GetExtendedErrorInfo())

  def __del__(self):
    self.__class__.NIDAQ_CREATED = False

class Device:
  def __init__(self, dev, nidaq):
    self.dev = dev
    self.nidaq = nidaq

  def createTask(self, taskName=""):
    return Task(self, taskName)

  def getType(self):
    return self.nidaq.GetDevProductType(self.dev)
    
  def getSerialNumber(self):
    sn = uInt32(0)
    self.nidaq.GetDevSerialNum(self.dev, byref(sn))
    return sn.value

  def listAIChannels(self):
    return self.nidaq.GetDevAIPhysicalChans(self.dev).split(b", ")

  def listAOChannels(self):
    return self.nidaq.GetDevAOPhysicalChans(self.dev).split(b", ")

  def listDILines(self):
    return self.nidaq.GetDevDILines(self.dev).split(b", ")

  def listDIPorts(self):
    return self.nidaq.GetDevDIPorts(self.dev).split(b", ")

  def listDOLines(self):
    return self.nidaq.GetDevDOLines(self.dev).split(b", ")

  def listDOPorts(self):
    return self.nidaq.GetDevDOPorts(self.dev).split(b", ")

class Task:
  TaskHandle = uInt32
  
  def __init__(self, device, taskName=""):
    self.handle = Task.TaskHandle(0)
    self.device = device
    self.nidaq = device.nidaq
    self.nidaq.CreateTask(taskName,byref(self.handle))

  def __del__(self):
    self.nidaq.ClearTask(self.handle)

  def __getattr__(self, attr):
    func = getattr(self.nidaq, attr)
    return lambda *args: func(self.handle, *args)

  #def addChannel(self, channelType, *args):
    #"""channelType must be named to fit the API for their channel creation functions."""
    #self.nidaq.call("DAQmxCreate%sChan" % channelType, self.handle, *args)

  #def addTiming(self, timingType, *args):
    #self.nidaq.call("DAQmxCfg%sTiming" % timingType, self.handle, *args)

  def start(self):
    self.nidaq.StartTask(self.handle)

  def stop(self):
    self.nidaq.StopTask(self.handle)

  def isDone(self):
    return self.nidaq.IsTaskDone(self.handle)

  def read(self, samples=None, timeout=10.):
    if samples is None:
      samples = self.GetSampQuantSampPerChan()
    numChans = self.GetTaskNumChans()
    
    shape = (numChans, samples)
    #print "Shape: ", shape
    buf = numpy.empty(shape, dtype=float64)
    samplesRead = int32()
    self.ReadAnalogF64(samples, timeout, Val_GroupByChannel, buf.ctypes.data, buf.size, byref(samplesRead), None)
    return (buf, samplesRead.value)

  def write(self, data, timeout=10.):
    numChans = self.GetTaskNumChans()
    samplesWritten = int32()
    self.WriteAnalogF64(data.size / numChans, False, timeout, Val_GroupByChannel, data.ctypes.data, byref(samplesWritten), None)
    return samplesWritten.value

init()
