__author__ = "Tessa Ropp, Paul B. Manis"
__version__ = "0.50"

import pysound.pysounds
import pysound.sound


import platform
# only import these under windows
if platform.system() == "Windows":
    from nidaq import *
    from cheader import *

