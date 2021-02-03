__author__ = "Tessa Ropp, Paul B. Manis"
__version__ = "0.50"

import pysound.pysounds
import pysound.sound


import platform
# only import these under windows
nidaq_available = False
if platform.system() == "Windows":
    try:
    	from nidaq import *
    	from cheader import *
    	nidaq_available = True
    except:
    	pass
