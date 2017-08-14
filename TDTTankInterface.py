from __future__ import print_function

import sys
import os
import datetime
import numpy as np
from PyQt4 import QtGui, QtCore 
import time
from collections import OrderedDict
import pprint

if os.name == 'nt':
    import win32com.client

class TDTTankInterface():
    def __init__(self):
        if os.name == 'nt':
            self.available = True
        else:
            self.available = False
        self.tank_directory = ''
        print('setup tank interface')
    
    def show_tank_path(self):
        print('Tank path: ', self.tank_directory)

    def open_tank(self):
        self.TDT_Tank = win32com.client.Dispatch('TTank.X')
        self.TDT_Tank.ConnectServer('Local','Me')
        self.TDT_Tank.OpenTank(self.tank_directory,'R')  # set with path to the tank

        
    def find_last_block(self):
        if not self.available:
            return 0
        blocknum = 0
        while self.TDT_Tank.QueryBlockName(blocknum) != '':
            blocknum = blocknum + 1
        #blocknum = blocknum-1
        return blocknum
    
    def close_tank(self):
        self.TDT_Tank.CloseTank()
    
    def release_server(self):
        self.TDT_Tank.ReleaseServer()
        
    def add_tank(self, tankname):
        TT.AddTank(tankname, self.tank_directory)
        