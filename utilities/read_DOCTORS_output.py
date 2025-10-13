# -*- coding: utf-8 -*-
"""
Read DOCTORS output binary file

DOCTORS output binary file structure:
    number of energy groups 'int'
    number of x elements    'int'
    number of y elements    'int'
    number of z elements    'int'
    x node positions        'float'
    y node positions        'float'
    z node positions        'float'
    flux value              'double'

@auther: Xin Liu
"""

import numpy as np
import matplotlib.pyplot as plt
import struct

# The following code is to read DOCTORS collided flux output
with open("collided_flux_iso_gpu_0.dat", "rb") as bin_file:
    data = bin_file.read()
    gcount, xcount, ycount, zcount = struct.unpack("iiii", data[:16])
    xnode = struct.unpack("f" * xcount, data[16:4*xcount + 16])
    ynode = struct.unpack("f" * ycount, data[4*xcount+16:4*xcount+4*ycount + 16])
    znode = struct.unpack("f" * zcount, data[4*xcount+4*ycount+16:4*xcount+4*ycount+4*zcount + 16])
    cflux = struct.unpack("d" * xcount * ycount * zcount * gcount, data[4*xcount+4*ycount+4*zcount+16:])
    
cflux = np.asarray(cflux)

cflux_reshaped = cflux.reshape((gcount, xcount, ycount, zcount))

# Plot 2D slice at specified energy group
plt.figure()
plt.title("Collided flux center slice at energy group 10")
plt.imshow(np.squeeze(cflux_reshaped[10,:,:,7]))
plt.colorbar()

    
    

