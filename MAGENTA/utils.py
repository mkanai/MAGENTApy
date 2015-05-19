import functools
import numpy as np
import sys

BUFLINENUM = 65536
MAXIDLEN = 128
MAXGENESETNUM = 256


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


def inHLAregion(HumanGeneChrPos, st, en):
    chrom = HumanGeneChrPos[:, 0]
    start = HumanGeneChrPos[:, 1]
    end = HumanGeneChrPos[:, 2]

    return np.logical_and(np.equal(chrom, 6),
                          logical_or(logical_and(np.less(start, st), np.greater(end, st), np.less_equal(end, en)),
                                     np.logical_and(np.greater_equal(start, st), np.less_equal(end, en)),
                                     logical_and(np.greater_equal(start, st), np.less(start, en), np.greater(end, en)),
                                     np.logical_and(np.less_equal(start, st), np.greater_equal(end, en))))
