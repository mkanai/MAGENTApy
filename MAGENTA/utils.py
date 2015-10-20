import functools
import numpy as np
import os
import re
import sys

BUFLINENUM = 65536
MAXIDLEN = 256
MAXGENENUM_PER_GENESET = 16384


class Writer(object):
    def __init__(self, filename, bufsize=-1):
        self.file = open(filename, 'w', bufsize)

    def write(self, msg, *args, **kwargs):
        msg = msg.format(*args, **kwargs)
        self.file.write(msg)


class Logger(Writer):
    def __init__(self, filename, bufsize=1):
        Writer.__init__(self, filename, bufsize)

    def log(self, msg, *args, **kwargs):
        msg = msg.format(*args, **kwargs)
        self.write(msg)
        sys.stdout.write(msg)


def logical_and(*args):
    return functools.reduce(np.logical_and, args)


def logical_or(*args):
    return functools.reduce(np.logical_or, args)


def format_NaN(format_string, value, sep='\t'):
    return 'NaN' + sep if np.isnan(value) else format_string.format(value)


def get_valid_filename(filename, dir=''):
    return os.path.join(dir, filename.replace('/', '_'))


def inHLAregion(HumanGeneChrPos, HLA_start, HLA_end):
    chrom = HumanGeneChrPos[:, 0]
    start = HumanGeneChrPos[:, 1]
    end = HumanGeneChrPos[:, 2]

    return np.logical_and(np.equal(chrom, 6),
                          np.greater(np.fmin(end, HLA_end)
                                     - np.fmax(start, HLA_start), 0))

