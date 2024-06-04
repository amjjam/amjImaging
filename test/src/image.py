import numpy

def read(filename):
    with open(filename,"rb") as f:
        nR=int.from_bytes(f.read(4),byteorder='little')
        nC=int.from_bytes(f.read(4),byteorder='little')
        return numpy.fromfile(f,(numpy.float64,(nR,nC)),-1)[0]
