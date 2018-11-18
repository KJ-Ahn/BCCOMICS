## Converts bccomics_setup-generated binaries into enzo ICs.
## Author: Britton Smith

import h5py
import numpy as np
import struct

def read_binary_file(filename):
    block_size = 32768
    dsize = 8
    f = open(filename, "rb")
    f.seek(0, 2)
    file_size = f.tell()
    f.seek(0)
    nblocks = np.ceil(float(file_size) /
                      block_size).astype(np.int64)
    array_size = int(file_size / dsize)
    data = np.empty(array_size, dtype=np.float64)
    i = 0
    offset = 0
    for ib in range(nblocks):
        my_block = min(block_size, file_size - offset)
        buff = f.read(my_block)
        for ibuff in range(int(my_block / dsize)):
            data[i] = struct.unpack("d", buff[ibuff*dsize:(ibuff+1)*dsize])[0]
            i += 1
        offset += my_block
    f.close()
    return data

def convert_file(input_filename, output_filename,
                 grid=True):
    print ("Converting %s to %s." % (input_filename, output_filename))
    dtype = np.int64
    if not isinstance(input_filename, tuple):
        input_filename = (input_filename,)
    data = [read_binary_file(filename)
            for filename in input_filename]

    crank = len(data)
    if grid:
        rank = 3
        csize = data[0].size
        dim = np.int(np.round(csize**(1./rank)))
        shape = dim * np.ones(rank, dtype=dtype)
        topgriddims = shape
        end = shape
        for i in range(len(data)):
            data[i] = np.reshape(data[i], shape, order="C")
    else:
        rank = 1
        csize = data[0].size
        shape = np.array([csize], dtype=dtype)
        topgriddims = -99999 * np.ones(3, dtype=dtype)
        end = topgriddims - 1
#    start = np.zeros(rank, dtype=dtype)
    start = np.zeros(3, dtype=dtype)
    data = np.array(data)

    fh = h5py.File(output_filename, "w")
    # dataset should have the same name as the file
    dataset = fh.create_dataset(output_filename, data=data)
    # write the attributes
    dataset.attrs["Component_Rank"] = crank
    dataset.attrs["Component_Size"] = csize
    dataset.attrs["Rank"]           = rank
    dataset.attrs["Dimensions"]     = getattr(shape, "tolist", lambda x=shape: x)()
    dataset.attrs["TopGridDims"]    = getattr(topgriddims, "tolist", lambda x=topgriddims: x)()
    dataset.attrs["TopGridEnd"]     = getattr(end, "tolist", lambda x=end: x)()
    dataset.attrs["TopGridStart"]   = getattr(start, "tolist", lambda x=start: x)()
    fh.close()

if __name__ == "__main__":
    convert_file(("cpos1", "cpos2", "cpos3"), "ParticlePositions", grid=False)
    convert_file(("vc1", "vc2", "vc3"), "ParticleVelocities", grid=False)
    convert_file(("vb1", "vb2", "vb3"), "GridVelocities", grid=True)
    convert_file("db", "GridDensity", grid=True)
    convert_file("etherm", "GasThermalSpecEnergy", grid=True)
    convert_file("etot", "GasTotalSpecEnergy", grid=True)
