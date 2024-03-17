'''
@author: Naoto Hori
'''

class SisoutHeader(object):
    def __init__(self) :
        self.step = None
        self.label = None
        self.temp = None

        self.e_kin = None
        self.e_pot = None
        self.e_bond = None
        self.e_angl = None
        self.e_dih = None
        self.e_bp = None
        self.e_exv = None
        self.e_ele = None
        self.e_stage = None
        self.e_twz = None

    def show(self):
        print('step', self.step)
        print('label', self.label)
        print('temp', self.temp)

        print('e_kin', self.e_kin)
        print('e_pot', self.e_pot)
        print('e_bond', self.e_bond)
        print('e_angl', self.e_angl)
        print('e_dih', self.e_dih)
        print('e_bp', self.e_bp)
        print('e_exv', self.e_exv)
        print('e_ele', self.e_el)
        print('e_stage', self.e_stage)
        print('e_twz', self.e_twz)

class SisoutFile(object):
    def __init__(self, filename) :
        self._filename = filename
        self.header_lines = []
        self.head_str = None
        self.head_col = None
        self.flg_u_u = False
        self.flg_all = False
        self.num_unit = 0

    def open_to_read(self):
        self._file = open(self._filename, 'r')

    def open_to_write(self):
        self._file = open(self._filename, 'w')

    def close(self):
        self._file.close()

    def flush(self):
        self._file.flush()

    def read_header(self):
        ''' store lines to self.header_lines
            store information to self.head_str and self.head_col
        '''
        if not self._file :
            raise Exception('_file is None in sisout.py')

        self._file.seek(0)

        line = self._file.readline()

        self.header_lines.append(line)
        self.head_str = line.split()
        self.head_col = self._head_str2col(self.head_str)

    def write_header(self):
        for l in self.header_lines:
            self._file.write(l)

    def copy_header(self, source):
        import copy
        self.header_lines = copy.deepcopy(source.header_lines)
        self.head_str = copy.deepcopy(source.head_str)
        self.head_col = copy.deepcopy(source.head_col)
        self.flg_u_u = copy.deepcopy(source.flg_u_u)
        self.num_unit = copy.deepcopy(source.num_unit)

    def read_onestep(self):
        lines = []
        lines.append( self._file.readline() )

        sisout_list = []

        #if self.flg_all:
        #    nlines = self.num_unit + 1
        #else:
        #    nlines = self.num_unit
        nlines = 1

        #if self.flg_u_u:
        #    #for i in xrange(self.num_unit+1):
        #    #for i in xrange(self.num_unit):  # There is no "#all" line
        #    for i in range(nlines):
        #        l = self._file.readline()
        #        lines.append(l)
        #        sisout_list.append(l.split()[1:])

        #    for i in range(self.num_unit*(self.num_unit-1)/2):
        #        l = self._file.readline()
        #        lines.append(l)
        #        sisout_list.append(l.split()[1:])
        #else:
        #    sisout_list.append( lines[-1].split() )
            #for i in xrange(self.num_unit+1):
            #    l = self._file.readline()
            #    lines.append(l)
            #    sisout_list.append(l.split())
        #    for i in range(nlines):
        #        l = self._file.readline()
        #        lines.append(l)
        #        sisout_list.append(l.split()[1:])

        sisout_list.append( lines[-1].split() )
        #for i in range(nlines):
        #    l = self._file.readline()
        #    lines.append(l)
        #    sisout_list.append(l.split()[1:])

        return (sisout_list, lines)

    def write_onestep(self, lines):
        for l in lines:
            self._file.write(l)

    def skip_onestep(self):
        for i in range(self.num_unit +2):
            self._file.readline()

    def skip(self, num):
        for i in range((self.num_unit +2)*num):
            self._file.readline()

    #def write_onestep(self, coord_matrix):

    def has_more_data(self):
        """return True or False"""
        s = self._file.tell()
        try:
            if self._file.readline() == '':
                self._file.seek(s)
                return False
            else:
                self._file.seek(s)
                return True
        except:
            self._file.seek(s)
            return False

    def _head_str2col(self, strlist):
        h = SisoutHeader()
        for i, s in enumerate(strlist):
            if s.endswith(')nframe'):
                h.step = i
            elif s.endswith(')R'):
                h.label = i
            elif s.endswith(')T'):
                h.temp = i
            elif s.endswith(')Ekin'):
                h.e_kin = i
            elif s.endswith(')Epot'):
                h.e_pot = i
            elif s.endswith(')Ebond'):
                h.e_bond = i
            elif s.endswith(')Eangl'):
                h.e_angl = i
            elif s.endswith(')Edih'):
                h.e_dih = i
            elif s.endswith(')Ebp'):
                h.e_bp = i
            elif s.endswith(')Eexv'):
                h.e_exv = i
            elif s.endswith(')Eele'):
                h.e_ele = i
            elif s.endswith(')Estage'):
                h.e_stage = i
            elif s.endswith(')Etweezers'):
                h.e_stage = i
            else:
                raise Exception('unknown header in sisout.py')
        return h
