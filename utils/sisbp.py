#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@author: Naoto Hori
'''

import os
import struct
from lop.elements.error import MyError

class SisbpFile :
    def __init__(self, filename) :
        self._filename = filename
        self._seek_data = None
        self._seek_mark = None

        self._int_kind = None
        self._real_kind = None
        self._kind_str = None
        self._step_byte = None

    def open_to_read(self):
        self._file = open(self._filename, 'rb')
        self.read_header()

    def open_to_write(self):
        self._file = open(self._filename, 'wb')

    def close(self):
        self._file.close()

    def flush(self):
        self._file.flush()

    def read_header(self):
        if not self._file :
            raise MyError('SisbpFile', 'read_header', 'Logical: _file is None')

        self._file.seek(0)

        # Read kind
        num1, num2 = struct.unpack('ii', self._file.read(8))
        if num1 == 2:
            self._int_kind = 2
            int_kind_str = 'h'
        elif num1 == 4:
            self._int_kind = 4
            int_kind_str = 'i'
        elif num1 == 8:
            self._int_kind = 8
            int_kind_str = 'q'
        else:
            raise MyError('SisbpFile', 'read_header', 'unknown int kind type: '+'%i'%num1)

        if num2 == 4:
            self._real_kind = 4
            real_kind_str = 'f'
        elif num2 == 8:
            self._real_kind = 8
            real_kind_str = 'd'
        else:
            raise MyError('SisbpFile', 'read_header', 'unknown real kind type: '+'%i'%num2)

        self._step_byte = 2 * self._int_kind + self._real_kind
        self._kind_str = int_kind_str + int_kind_str + real_kind_str
        self._seek_data = self._file.tell()

    def set_header(self, header):
        (int_kind, real_kind, kind_str, step_byte) = header
        self._int_kind = int_kind
        self._real_kind = real_kind
        self._kind_str = kind_str
        self._step_byte = step_byte

    def get_header(self):
        return (self._int_kind,
                self._real_kind,
                self._kind_str,
                self._step_byte)

    def copy_header(self, src):
        self._int_kind = src._int_kind
        self._real_kind = src._real_kind
        self._kind_str = src._kind_str
        self._step_byte = src._step_byte

    def write_header(self):
        if not self._file :
            raise MyError('SisbpFile', 'write_header', 'Logical: _file is None')

        self._file.write(struct.pack('ii', self._int_kind, self._real_kind))

    def read_onestep(self):

        pairs = []
        energies = []

        while True:

            i, j, e = self._pick_data()
            if i == 0:
                break

            pairs.append((i,j))
            energies.append(e)

        return pairs, energies

    def write_onestep(self, pairs, energies):

        form = self._kind_str  # 'hhf'
        for (i, j), energy in zip(pairs, energies):
            self._file.write(struct.pack(form, i, j, energy))

        self._file.write(struct.pack(form, 0, 0, 0.0))

    def skip_onestep(self):
        try:
            while(True):
                if (self._pick_data()[0] == 0):
                    break
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

    def has_more_data(self):
        """return True or False"""
        char = self._file.read(self._int_kind)
        if not char :
            return False
        else :
            self._file.seek(-self._int_kind, os.SEEK_CUR)
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
        return struct.unpack(self._kind_str, self._file.read(self._step_byte))

    def _put_data(self, binary, size):
        self._file.write(struct.pack('<L', size))
        self._file.write(binary)
        self._file.write(struct.pack('<L', size))
        # '<' indicates little endian, 'L' stands 4 byte data

    def _read_at(self, num):
        self._file.seek(0)
        self.read_header()
        for i in range(num - 1) :
            self.skip_onestep()
        return self.read_onestep()

