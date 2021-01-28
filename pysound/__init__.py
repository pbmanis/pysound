

import platform
# only import these under windows
if platform.system() == "Windows":
    from nidaq import *
    from cheader import *

