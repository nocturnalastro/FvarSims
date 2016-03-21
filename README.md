# FvarSims

##TODO
* make use of pyxspec so we can remove the need to make loads of files.
* Add prefix to lc log files 
* remove uncertanites as its not needed - then test it still works

## Install
I haven't packaged this!

##Example

Firstly an xspec commands file must be made which will load the data and model into xspec (init.xcm in exmaple).
You must note down the parameters you wish to vary by parameter number and the ranges you wish to allow them to vary between.

OK so import the required functions

`
from simulations import make_spectra, make_count_analysis, run
`

To produce the fake spectra we run:

```
sxcm,specs=make_spectra(numspec=1600,init_file="init.xcm",params={2:(0.2,0.6),3:(1.9,2.3),4:(0.1,0.8)},dist=uni,exposure=11520, bkgfile="bkg.pi", rspfile="rmf.rmf",ancrfile="arf.arf",numxcm=40, bkgscale=0.25, prefix="big_")
run(40,sxcm)
```

This will sample the distribution uni(a wrapper for numpy.uniform) between the ranges defined in params for each parameter this produce 40 xcm files. Each xcm will produce a number of fake spectrum big_lc*.fak each with a single realisation of the parameters. `specs` is a list of filenames,(param values,...).  Then `run` will set running 40 instances of xspec which will run one of the produced xcm files.

We can then make xcm files which will then sample the `count rate` of the fake spectrum in given energy bands to with the aim to produce lightcurves from these `count rates`

```
lxcm=make_count_analysis((0.5,2,5,10),80, specs)
run(40,sxcm)
```

This will produce the files `lc0.log->lc79.log` (I havent put a prefix for lc*.log files).
To convert the log files into something useful like calculating the fractional variabilty (this is mostly written in the `main` for a vary basic case)

```
from simlog2fvar import logs_to_fvar, dofvar, nonan, unp_to_tuple, nanmean, write_fvar
ens, fv = logs_to_fvar(prefix="lc", dir=".", dofunc=dofvar) #not needed in this case but as an example
```
**N.B** if you just want the something different to the lightcurves you can change dofunc to something like `lambda x: x` so the lightcurve is returned (there is a function `make_lcs` in simulatons.py but it only does one file - bad naming I know).

So now we have alot of realistions of Fvar spectra so lets take the mean. however some of the lightcurves can produced `NaN`s due to the nature of Fvar, so lets throw them away and take the mean.

```
fv = nonan(fv).T
fmean,_ unp_to_tuple(nanmean(fv))
```

Error bars are nice! 
So lets calculate them from the width of the distribution. Also it's probably a good idea to plot them to make sure they don't look odd.
```
dists = calc_dist(fv, 8) # we only have 80 fvars so not really all tha usefull
plt.ion()
plot_dist_fit(dists,1,4) #1 row 4 columns  or the other way round I can't remember
sigmas = get_dist_sigma(dists)
```
So lets write this to a handy ascii file!
```
write_fvar(ens, zip(fmean, sigs), prefix="big_mean")
```
