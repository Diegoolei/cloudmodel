import struct
import numpy as np

# Python example to read a cloud file
f = open('outputdata1/nube3101.sal', 'rb')
recl = struct.unpack('i', f.read(4))[0]
numval = recl/np.dtype('float32').itemsize
data = np.fromfile(f, dtype='float64', count=-1)
endrec = struct.unpack('i', f.read(4))[0]
if endrec is not recl:
    print("error unexpected end rec")
f.close()
