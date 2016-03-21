#! /usr/bin/env python
from inspect import getargspec
import numpy as np
import subprocess
from multiprocessing import pool
import re


def thenumbers(values, numsplit):
    """takes in a list of values then splits them into {numsplit} lists."""
    to_return = [[] for _ in xrange(numsplit)]
    flor = len(values) // numsplit
    remain = len(values) % numsplit
    remains_done = 0
    index = 0
    for i, val in enumerate(values):
        to_return[index].append(val)
        # (num persplit + a remainder)*index
        if i == (index + 1) * flor + (index < remain) + remains_done - 1:
            if index < remain:
                remains_done += 1
            index += 1
    return to_return


def uni(size):
    # a function to wrap numpy uniform so make_spectra will understand it.
    return np.random.uniform(size=size)


def make_spectra(numspec, init_file, params, dist, exposure, bkgfile, rspfile, ancrfile="", numxcm=1, bkgscale=0.25, prefix=""):
    """This function makes the spectra for the simulation.

    numspec(int)          : The number of spectra you wish to simulate this will be used as the value for size in dist.(number of data points * number of lightcurves)
    init_file(str)        : The file which is the intal setup for the model.
    params(dict)          : A dictionary the key is the  model parameter to be changed the value is a (min,max) tuple which defines the range.
    dist(func)            : A function which returns a value from 0->1 which defines the distibution it must take a minimum of the number of values to return
                           (if 1 input parameter with no defualts -- e.g.foo(size))
                           (if 2 input parameters with no defualts -- foo(size,range))
                           (if 3 input parameters with no defualts foo(size,min,max))
                           (more input parameters are allowed but they must have defualts).
                           Any seed should be taken care of in dist.
                           Also numpy.random fuctions and lambdas wont work as the inspect module doesn't play well with them, so wrap them like this:
                           def uni(size):
                             return np.random.uniform(size=size)

    exposure(float)       : The exposure time for each simulation (note: if there are two detectors coadded you may need to double it).
    bkgfile(str)          : The input background file.
    rspfile(str)          : The input responce file - if no ancrfile is given this should be a rsp not and rmf.
    ancrfile(str)('')     : The arf file.
    numxcm(int)           : The number of xcm files to produce.
    bkgscale(float)(0.25) : The background scaling factor (effectively the ratios of the effective area of the source and background regions)
    prefix(str)('')       : The prefix prepend the lc_{number}.fak file.

    Returns:(xcms,spectra)
    xcms              : a list of xcm files made
    spectra           : a list of spectra and the input parameters used to make them [(spectra,(parameters,)),...]
    """
    simulated_values = {}
    inspec = getargspec(dist)
    numargs = len(inspec.args) - len(inspec.defaults)
    if numargs == 1:
        for parameter, val_range in params:
                # dist should only return values between 0 and 1 so we mulitply it by the range of values required
                simulated_values[parameter] = min(val_range) + np.array(dist(numspec)) * (max(val_range) - min(val_range))
    elif numargs == 2:
        for parameter, val_range in params:
            # dist takes a range so it sould return values within that range
            simulated_values[parameter] = np.array(dist(numspec, (min(val_range), max(val_range))))
    elif numargs == 3:
        for parameter, val_range in params:
            # aren't I nice giving you an option of how to get out your function
            # dist takes a range so it sould return values within that range
            simulated_values[parameter] = np.array(dist(numspec, min(val_range), max(val_range)))
    else:
        raise Exception("Dist has toooooo meny arguments!")

    the_val_indexs = thenumbers(xrange(numspec), numxcm)
    spectra = []
    xcms = []
    for xcm in xrange(numxcm):
        with open(prefix + "mc_%i.xcm" % xcm, "w") as fw:
            xcms.append("mc_%i.xcm" % xcm)
            # put the init_files contents into new xcm
            fw.write(open(init_file).read() + "\n")
            for val_index in the_val_indexs[xcm]:
                for param, values in simulated_values.iteritem():
                    fw.write("new %i %f\n" % (param, values[val_index]))
                fw.write("fakeit %(bkg)s\n%(rsp)s\n%(arf)s\ny\n\n%(pre)slc%(val)i.fak\n%{exp}f,%{scale}f,%{exp}f\n\n\n" %
                         {"bkg": bkgfile, "rsp": rspfile, "arfs": ancrfile, "pre": prefix, "val": val_index, "exp": exposure, "scale": bkgscale})
                spectra.append(("%(pre)slc%(val)i.fak" % {"pre": prefix, "val": val_index}, ((param, values[val_index]) for param, values in simulated_values.iteritem())))
                fw.write("data none\n")  # must clear spectra or fake it goes a bit mental
            fw.write("exit\n\n\n")  # probably only need two \n but more don't hurt :)
    return xcms, spectra


def make_count_analysis(energy_bins, numxcm, spectra="", numspec=0, prefix=""):
    """Makes the files which analyse the simulated spectra. Each log file created should be a single lc.

    energy_bins(list) : The energy_bins in keV of the format (0.5,1,2,3,4) which are then made into ((0.5,1),(1,2),(2,3),(3,4)).
    numxcm(int)      : The number of files you wish to create number. The number of spectra should be divisable by this number otherwise you have wasted spectra.
    spectra(list)     : The output of make_spectra -or atleast just filenames of the spectra - if not set you must set numspec,prefix

    Required if spectra not set:
    numspec(int)      : The max value of spectra you wish to use(non inclusive so if 100 is given the spectra used are 0-99).
    prefix(str)('')   : The prefix prepend to the lc_{number}.fak file.

    Returns:
    xcms              : a list of xcm files made
    """
    assert (spectra is not "" and (isinstance(spectra, list) or isinstance(spectra, tuple)) and len(spectra) > 0) \
        or (numspec is not 0 and prefix is not ""), \
        "spectra or (numspec and prefix must be set: How am I ment to know what files to use?)"

    if spectra is not "":
        if len(spectra[0]) > 1:
            # _spectra = spectra #save spectra incase I need it later
            spectra = [spec[0] for spec in spectra]
    else:
        spectra = ["%(pre)slc%(val)i.fak" % {"pre": prefix, "val": i} for i in xrange(numspec)]

    spec_per_file = len(spectra) // numxcm
    enbins = zip(energy_bins[:-1], energy_bins[1:])
    xcms = []
    for n in xrange(numxcm):
        fw = open("lc_%i.xcm" % n, 'w')
        xcms.append("lc_%i.xcm" % n)
        fw.write("log lc_%i.log\n" % n)
        fw.write("data " + " ".join(["%(dgroup)i:%(dgroup)i %(spec)s" % {"dgroup": m + 1, "spec": spec} for m, spec in enumerate(spectra[n * spec_per_file:(n + 1) * spec_per_file])]) + "\n")
        for elow, ehigh in enbins:
            fw.write("not *:0.0-**\n")  # notice all channels
            fw.write("ign *:0.0-%.2f %.2f-**\n" % (elow, ehigh))  # filter energies
            fw.write("show data\n")  # this is where all the interesting info is!
        fw.write("exit\n\n\n")
    return xcms


def run(numpara, xcms, heainit="heainit", shell='tcsh', seperator=';'):
    """Runs the xcm files created by the make* functions.

    numpara(int)        : The number of processes which are spawned.
    xcm(list)           : The xcm files which are ran.
    shell(str)('tcsh')  : The shell you wish to run on.
    heainit(str)        : The command used to initate heasoft
    seperator(str)(';') : The seperator used to end a command
    """
    def _run(cmd):
        """runs cmds in shell."""
        return subprocess.Popen(cmd, shell=True, executable='/bin/' + shell, stdin=subprocess.PIPE, stdout=subprocess.PIPE)

    p = pool.Pool(40)
    p.map(run, [seperator.join([heainit, xcm]) for xcm in xcms])


def make_lcs(logfile):
    """Take a logfile and returns a lighcurve/count rates for each energy.

    Returns
    lcs : A dict of energy bin, lightcuve pairs {(elow,ehigh):[[count rate],[errors]]}.
    """
    count_regex = re.compile("Net count rate \(cts/s\) for Spectrum:\d+\s+(-?(?:\d*\.)?\d+e?[-|+]?\d*)\s+\+/-\s+(-?(?:\d*\.)?\d+e?[-|+]?\d*).*")
    ign_regex = re.compile("\d+>ign \*\:0.0-((?:\d*\.)?\d+e?-?\d*) ((?:\d*\.)?\d+e?-?\d*)-\*\*;")
    lc = {}  # this houses the output lightcurves lc[energy]=[counts],[error]
    # this is voodoo which collects only count rates and ignore commands from the log and seperates them in to energies effectively
    for en_counts in "".join([l for l in open(logfile).read().split("\n") if "Net" in l or "ign " in l]).split("!XSPEC")[1:]:
        energy, counts = en_counts.split("#", 1)
        counts = np.array([map(float, count_regex.match(c).groups()) for c in counts.split("#")])
        energy = tuple(map(float, ign_regex.match(energy).groups()))
        lc[energy] = counts.T
    return lc
