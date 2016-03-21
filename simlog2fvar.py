#! /usr/bin/env python
import os
import numpy as np
from collections import defaultdict
from uncertainties import ufloat, nominal_value, std_dev
from uncertainties import unumpy as unp
from simulations import make_lcs, dofvar
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit


def logs_to_fvar(prefix, dir=".", dofunc=dofvar):
    """Take all the log files in a dir with {prefix} in the name and returns a list of (energies,error) and a numpy array of count rates.

    prefix(str)  : log files should have 'lc', 'log' and '{prefix}' in it
    dir(str)     : the directory to look for the log files
    dofunc(func) : the function to apply to the lightcurves

    Return:
    fvar_ens,fvar_values
    """
    fvars = defaultdict(list)
    for logfile in [l for l in os.listdir(dir) if "log" in l and "lc" in l and prefix in l]:
        for key, value in make_lcs(logfile).iteritems():
            fvars[key].append(ufloat(*dofunc(*value)))
    fvar_ens = [(sum(key) / 2, round(abs(key[1] - key[0]), 4) / 2) for key in fvars.keys()]
    fvar_values = np.dstack(fvars.values())[0].T
    return fvar_ens, fvar_values


def logs_to_fv_multi(prefix, dir=".", dofunc=dofvar):
    """Take all the log files in a dir with {prefix} in the name and returns a list of (energies,error) and a numpy array of count rates.

    prefix(list): log files should have 'lc', 'log' and one of '{prefix}' in it this will concatinate lcs then input them into dofunc
    dir(str)   : the directory to look for the log files
    dofunc(func) : the function to apply to the lightcurves

    Return:
    fvar_ens,fvar_values
    """
    fvars = defaultdict(list)
    for logfile in zip(*[sorted([l for l in os.listdir(dir) if "log" in l and "lc" in l and p in l]) for p in prefix]):
        lcs = [make_lcs(log) for log in logfile]
        for key, value in lcs[0].iteritems():
            value = np.hstack([l[key] for l in lcs])
            fvars[key].append(ufloat(*dofunc(*value)))
    fvar_ens = [(sum(key) / 2, round(abs(key[1] - key[0]), 4) / 2) for key in fvars.keys()]
    fvar_values = np.dstack(fvars.values())[0].T
    return fvar_ens, fvar_values


def unp_to_tuple(fvals):
    """Takes a list/array of ufloats and returns a list of (value,error)."""
    fvar_vals = []
    for fv in fvals:
        fvar_vals.append((nominal_value(fv), std_dev(fv)))
    return fvar_vals


def nanmean(fvals):
    """Does nanmean ... numpy was giving a strange result.

    fvals(list): the input values -> rows will be averaged over
    """
    fvar_vals = []
    for fv in fvals:
        fvar_vals.append(np.mean(fv[unp.isnan(fv) is False]))
    return fvar_vals


def nanmed(fvals):
    """Does median of non-nan values.

    fvals(list): the input values -> rows will be averaged over
    """
    fvar_vals = []
    for fv in fvals:
        fvar_vals.append(np.median(fv[unp.isnan(fv) is False]))
    return fvar_vals


# These probably aren't useful [
def get_spec_mean(fv):
    """Returns spectrum closest to the mean."""
    return fv.T[np.argmin(abs(np.array(nanmean(fv.T)) - np.mean(nanmean(fv.T))))]


def get_spec_max(fv):
    return fv.T[np.argmax(nanmean(fv.T))]


def get_spec_min(fv):
    return fv.T[np.argmin(nanmean(fv.T))]
# ]


def nonan(fv):
    """Removes spectrums with nans in them."""
    ret = []
    for f in fv.T:
        if np.any(unp.isnan(f)):
            continue
        ret.append(f)
    return np.array(ret)


def write_fvar(fvar_ens, fvar_values, prefix=""):
    """Takes a list of energies (value,error) and fvars (values,error) and writes a file.

    fvar_ens(list)    : list of (energies,error) will make up column 1 and 2
    fvar_values(list) : list of (fvar,error) will make up column 3 and 4
    prefix(str)       : a prefix for output file
    """
    if prefix == "":
        prefix = os.getcwd().split("/")[-1]
    assert len(fvar_values) == len(fvar_ens)
    with open(prefix + "_fvar.dat", "w") as fw:
        fw.writelines(sorted(["%f %f " % fvar_ens[n] + "%f %f\n" % tuple(fvar_values[n]) for n in xrange(len(fvar_values))]))


def calc_dist(fvars, nbins):
    """Calulates and plots the ditributions of the simulated data. This is then fit with a gaussian so we can derive errors."""
    hist_grams = []
    for n, f in enumerate(fvars):
        fv = nonan(f).T
        fv = zip(*unp_to_tuple(fv))[0]
        dist, bins, patches = plt.hist(fv, nbins, normed=True)
        bins = (bins[1:] + bins[:-1]) / 2.0
        hist_grams.append((bins, dist, fit_dist(bins, dist)))
    return hist_grams


def gauss_func(x, *p):
    """just a guassian - not exciting."""
    A, mu, sigma = p
    return A * np.exp(-(x - mu) ** 2 / (2. * sigma ** 2))


def fit_dist(bins, dist):
    """Lets fit the distributions so we can figure out the errorbars."""
    p0 = [1, bins.mean(), bins.std()]
    coeff, var_matrix = curve_fit(gauss_func, bins, dist, p0=p0)
    params = dict(zip(["norm", "mu", "sigma"], coeff))
    params['sigma'] = abs(params['sigma'])
    return gauss_func(bins, *coeff), params


def plot_dist_fit(dists, i1, i2):
    """Assumes interactive plt."""
    indexs = np.vstack([[(x, y) for y in xrange(i2)] for x in xrange(i1)])
    fig, axs = plt.subplots(i1, i2)
    axs = np.array(axs)
    for n, (bins, dist, (gaus, params)) in enumerate(dists):
        ax = axs[tuple(indexs[n])]
        ax.plot(bins, gaus)
        ax.bar(bins, dist, np.mean((bins[1:] - bins[:-1]) / 2))


def get_dist_sigma(dists):
    """Assumes return from calc_dist is input and returns the sigma values from the dist coeff's."""
    sigmas = []
    for bins, dist, (gaus, params) in dists:
        sigmas.append(params['sigma'])
    return sigmas


def main(inp=False, prefix="", remove_nans=True):
    """This does a general treatment of the output logs and then writes a mean fvar spectrum to file.

    Prameters:
    ----------

    inp* :  energies and fvars in the form of ([engies],[[fvars for en 1],[...],...]), if False they're calculated from logs_to_fvar
    prefix* : prefix of the input files used if inp is False
    remove_nans*: if True spectrums with nans within them are ignored

    * optional
    """
    if inp is False:
        ens, fv = logs_to_fvar(prefix)
    else:
        ens, fv = inp
    if remove_nans:
        fv = nonan(fv).T

    dists = calc_dist(fv, 40)
    plot_dist_fit(dists, 3, 6)
    sigs = get_dist_sigma(dists)

    # the reason unp_to_tuple and nanmean are seperate
    # is cuz you might not want to take the mean
    np.save(os.getcwd().split("/")[-1] + "_allfvars", zip(ens, fv))
    fv2 = unp_to_tuple(nanmean(fv))
    fv2 = zip(zip(*fv2)[0], sigs)
    write_fvar(ens, fv2, prefix=os.getcwd().split("/")[-1])

    '''These where just exploratory but I might need them again some day
    fv2 = unp_to_tuple(nanmed(fv))
    write_fvar(ens, fv2, prefix=os.getcwd().split("/")[-1] + "_median")

    fv2 = unp_to_tuple(get_spec_mean(fv))
    write_fvar(ens, fv2, prefix=os.getcwd().split("/")[-1] + "_meanspec")

    fv2 = unp_to_tuple(get_spec_min(fv))
    write_fvar(ens, fv2, prefix=os.getcwd().split("/")[-1] + "_minspec")

    fv2 = unp_to_tuple(get_spec_max(fv))
    write_fvar(ens, fv2, prefix=os.getcwd().split("/")[-1] + "_maxspec")

    write_fvar(ens, unp_to_tuple(np.min(fv, axis=1)), prefix=os.getcwd().split("/")[-1] + "_min")
    write_fvar(ens, unp_to_tuple(np.max(fv, axis=1)), prefix=os.getcwd().split("/")[-1] + "_max")
    '''


if __name__ == '__main__':
    main()
