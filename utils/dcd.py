#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@author: Naoto Hori
'''

import os
import struct
from lop.elements.error import MyError

class DcdHeader(object):
    def __init__(self):
        self.format = 'cafemol'
        self.block1 = None
        self.nset = None
        self.istart = None
        self.nstep_save = None
        self.nstep = None
        self.nunit_real = None
        self.delta = None
        self.title = None
        self.tempk = None
        self.lunit2mp = None
        self.nmp_real = None
        self.with_unit_cell = False

    def show(self):
        print('format', self.format)
        print('The 1st block:')
        for i, data in enumerate(self.block1):
            print('   {:2d}:'.format(i+1), self.block1[i])
        if self.title is not None:
            for line in self.title :
                print(line.decode('utf-8'))
        print('nset', self.nset)
        print('istart', self.istart)
        print('nstep_save', self.nstep_save)
        print('nstep', self.nstep)
        print('nunit_real', self.nunit_real)
        print('delta', self.delta)
        print('tempk', self.tempk)
        if self.nunit_real is not None:
            for i in range(self.nunit_real) :
                print('lunit2mp[', i, ']', self.lunit2mp[i])
        print('nmp_real', self.nmp_real)

class DcdFile :
    def __init__(self, filename) :
        self._filename = filename
        self._header = DcdHeader()
        self._seek_data = None
        self._seek_mark = None

    def open_to_read(self):
        self._file = open(self._filename, 'rb')

    def open_to_write(self):
        self._file = open(self._filename, 'wb')

    def close(self):
        self._file.close()

    def flush(self):
        self._file.flush()

    def read_header(self):
        if not self._file :
            raise MyError('DcdFile', 'read_header', 'Logical: _file is None')

        ###########################################
        # Check the format reading title lines
        ###########################################

        self._file.seek(0)

        # Skip the first block
        b = self._pick_data()

        # title block 
        b = self._pick_data()
        size = len(b)
        nline = (size - 4) // 80
        bdata = struct.unpack(('i' + '80s' * nline), b)
        try:
            if bdata[1].decode("utf-8").find('Git commit') != -1 or bdata[2].decode("utf-8").find('CafeMol') != -1:
                self._header.format = 'cafemol'
            else:
                self._header.format = 'charmm'
        except:
            self._header.format = 'charmm'


        ###########################################
        # Load header infomation
        ###########################################

        self._file.seek(0)

        if self._header.format == 'cafemol':
            # first block 
            b = self._pick_data()
            bdata = struct.unpack('4siii5iifi9i', b)
            #bdata = struct.unpack('4siii5iid9i',b)

            self._header.block1 = bdata

            if bdata[0].decode("utf-8") != 'CORD' :
                raise MyError('DcdFile', 'read_header', '%s is not coordinate file' % (self._filename,))
            self._header.nset = bdata[1]
            self._header.istart = bdata[2]
            self._header.nstep_save = bdata[3]
            self._header.nstep = bdata[4]
            self._header.nunit_real = bdata[5]
            self._header.delta = bdata[10]

            # title block (number of line can be changed)
            b = self._pick_data()
            bdata = struct.unpack(('i' + '80s' * (3 + self._header.nunit_real)), b)
            self._header.title = (bdata[1], bdata[2])
            #self._header.tempk = float(bdata[3].strip('\0 '))
            self._header.tempk = float(bdata[3])
            self._header.lunit2mp = []
            for i in range(self._header.nunit_real) :
    #            self._header.lunit2mp.append(int(bdata[i + 4].strip('\0 ')))
                self._header.lunit2mp.append(int(bdata[i + 4]))
    
            # nmp_real
            b = self._pick_data()
            self._header.nmp_real = struct.unpack('i', b)[0]

        elif self._header.format == 'charmm':
            # first block 
            b = self._pick_data()
            #bdata = struct.unpack('4s20i', b)
            bdata = struct.unpack('4s9if10i', b)

            self._header.block1 = bdata

            if (bdata[11] == 1):
                self._header.with_unit_cell = True
            #bdata = struct.unpack('4siii5iid9i',b)
            #if bdata[0] != 'CORD' :
            #    raise MyError('DcdFile', 'read_header', '%s is not coordinate file' % (self._filename,))
            #self._header.nset = bdata[1]
            #self._header.istart = bdata[2]
            #self._header.nstep_save = bdata[3]
            #self._header.nstep = bdata[4]
            #self._header.nunit_real = bdata[5]
            #self._header.delta = bdata[10]

            # title block (number of line can be changed)
            b = self._pick_data()
            size = len(b)
            nline = (size - 4) // 80
            bdata = struct.unpack(('i' + '80s' * nline), b)
            self._header.title = []
            for i in range(nline):
                self._header.title.append(bdata[i+1])
    
            # nmp_real
            b = self._pick_data()
            self._header.nmp_real = struct.unpack('i', b)[0]

        self._seek_data = self._file.tell()


    def write_header(self):
        if not self._header :
            raise MyError('DcdFile', 'write_header', 'Logical: _header is None')

        self._file.seek(0)

        if self._header.nset is None:
            self._header.nset = 0
        if self._header.istart is None:
            self._header.istart = 0
        if self._header.nstep_save is None:
            self._header.nstep_save = 0
        if self._header.nstep is None:
            self._header.nstep = 0
        if self._header.nunit_real is None:
            self._header.nunit_real = 0
        if self._header.delta is None:
            self._header.delta = 0.0
        if self._header.tempk is None:
            self._header.tempk = 0.0
        if self._header.lunit2mp is None:
            self._header.lunit2mp = []

        #first block
        form = '4siii5iifi9i'
        binary = struct.pack(form,
                             b'CORD',
                             self._header.nset,
                             self._header.istart,
                             self._header.nstep_save,
                             self._header.nstep,
                             self._header.nunit_real,
                             0, 0, 0, 0,
                             self._header.delta,
                             0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        self._put_data(binary, struct.calcsize(form))

        #title block
        if self._header.format == 'cafemol':
            binary = struct.pack('i' + '80s' * 2,
                                 (3 + len(self._header.lunit2mp)),
                                 self._header.title[0],
                                 self._header.title[1])

            import re
            re_null = re.compile(b'\0')

            p = struct.pack('80s', bytes('%f' % self._header.tempk, 'utf-8'))
            binary += re_null.sub(b' ', p)
            for i in range(self._header.nunit_real) :
                p = struct.pack('80s', bytes('%i' % self._header.lunit2mp[i], 'utf-8'))
                binary += re_null.sub(b' ', p)

            form = 'i' + '80s' * (3 + self._header.nunit_real)
            self._put_data(binary, struct.calcsize(form))

        else:
            binary = struct.pack('i', len(self._header.title))
            for t in self._header.title:
                binary += struct.pack('80s', t)
            form = 'i' + '80s' * len(self._header.title)
            self._put_data(binary, struct.calcsize(form))

        # nmp_real
        binary = struct.pack('i', self._header.nmp_real)
        self._put_data(binary, 4)


    def show_header(self):
        """print header information"""
        self._header.show()

    def set_header(self, header):
        import copy
        self._header = copy.deepcopy(header)

    def get_header(self):
        return self._header

    def read_onestep(self):

        if (self._header.with_unit_cell):
            num = struct.unpack('i', self._file.read(4))[0]
            self._file.seek(4+num, os.SEEK_CUR)

        """return 2-dimensional lists"""
        coord_matrix = []
        b = self._pick_data()
        x = struct.unpack('f' * self._header.nmp_real, b)
        b = self._pick_data()
        y = struct.unpack('f' * self._header.nmp_real, b)
        b = self._pick_data()
        z = struct.unpack('f' * self._header.nmp_real, b)

        for i in range(self._header.nmp_real) :
            xyz = [x[i], y[i], z[i]]
            coord_matrix.append(xyz)

        return coord_matrix

    def read_onestep_np(self):
        """return ndarray"""
        import numpy as np
        data = np.empty((self._header.nmp_real, 3))
        b = self._pick_data()
        data[:,0] = struct.unpack('f' * self._header.nmp_real, b)
        b = self._pick_data()
        data[:,1] = struct.unpack('f' * self._header.nmp_real, b)
        b = self._pick_data()
        data[:,2] = struct.unpack('f' * self._header.nmp_real, b)

        return data

    def read_onestep_np_solute(self, nsolute):
        """return ndarray"""
        import numpy as np
        data = np.empty((self._header.nmp_real, 3))
        b = self._pick_data()
        data[:,0] = struct.unpack('f' * self._header.nmp_real, b)
        b = self._pick_data()
        data[:,1] = struct.unpack('f' * self._header.nmp_real, b)
        b = self._pick_data()
        data[:,2] = struct.unpack('f' * self._header.nmp_real, b)

        return data[:nsolute, 0:3]

    def read_onestep_npF(self):
        """return ndarray in Fortran format"""
        import numpy as np
        data = np.empty((3,self._header.nmp_real), order='F')
        b = self._pick_data()
        data[0,:] = struct.unpack('f' * self._header.nmp_real, b)
        b = self._pick_data()
        data[1,:] = struct.unpack('f' * self._header.nmp_real, b)
        b = self._pick_data()
        data[2,:] = struct.unpack('f' * self._header.nmp_real, b)

        return data

    def skip_onestep(self):
        try:
            if (self._header.with_unit_cell):
                num = struct.unpack('i', self._file.read(4))[0]
                self._file.seek(4+num, os.SEEK_CUR)
            num = struct.unpack('i', self._file.read(4))[0]
            self._file.seek(4+num, os.SEEK_CUR)
            num = struct.unpack('i', self._file.read(4))[0]
            self._file.seek(4+num, os.SEEK_CUR)
            num = struct.unpack('i', self._file.read(4))[0]
            self._file.seek(4+num, os.SEEK_CUR)
        except:
            raise EOFError('EOF in a middle of dcd.skip_onestep().')

    ''' This will throw EOFError exception if there are not enough frames.'''
    def skip(self, num):
        for _ in range(num):
            self.skip_onestep()

    ''' This does NOT throw EOFError exception if there are not enough frames.'''
    ''' Instead, this function will return the number of frames skipped successfully.'''
    def skip_as_many_as_possible_upto(self, num):
        for i in range(num):
            try:
                self.skip_onestep()
            except EOFError:
                return i
        return num

    def write_onestep(self, coord_matrix):
        # for X
        binary = b''
        for xyz in coord_matrix :
            binary += struct.pack('f', xyz[0])
        self._put_data(binary, 4 * self._header.nmp_real)
        # for Y
        binary = b''
        for xyz in coord_matrix :
            binary += struct.pack('f', xyz[1])
        self._put_data(binary, 4 * self._header.nmp_real)
        # for Z
        binary = b''
        for xyz in coord_matrix :
            binary += struct.pack('f', xyz[2])
        self._put_data(binary, 4 * self._header.nmp_real)

    def has_more_data(self):
        """return True or False"""
        char = self._file.read(4)
        if not char :
            return False
        else :
            self._file.seek(-4, os.SEEK_CUR)
            return True

    def count_frame(self):
        self.set_mark()

        if self._seek_data is None:
            self.read_header()

        self.rewind()

        n = 0
        while self.has_more_data():
            try:
                self.skip_onestep()
            except EOFError:
                break
            n += 1

        self.go_mark()

        return n

    def rewind(self):
        self._file.seek( self._seek_data, os.SEEK_SET)

    def set_mark(self):
        self._seek_mark = self._file.tell()

    def go_mark(self):
        self._file.seek( self._seek_mark, os.SEEK_SET)

    def _pick_data(self):
        """return binary data between 'integer' and 'integer'. 'integer' indicates the number of bytes"""
        num = struct.unpack('i', self._file.read(4))[0]
        b = self._file.read(num)
        self._file.seek(4, os.SEEK_CUR)
        return b

    def _put_data(self, binary, size):
        self._file.write(struct.pack('<L', size))
        self._file.write(binary)
        self._file.write(struct.pack('<L', size))
        # '<' indicates little endian, 'L' stands 4 byte data

    def _read_at(self, num):
        self._file.seek(0)
        self.read_header()
        for i in range(num - 1) :
            self.read_onestep()
        return self.read_onestep()

