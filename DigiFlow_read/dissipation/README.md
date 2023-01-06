# Internal Waves Energetics Calculations

These scripts are designed to calculate TKE dissipation from PIV fields, Written by Peter Sutherland for the Hydralab Internal Waves project. The primary objective was to study internal solitary wave propagation in ice covered waters. The scripts have now been modified by Sam Hartharn-Evans for extended use in the Novak Laboratory, and for advanced energetics analysis to study internal solitary wave propagation in ice-covered waters. 

## Script usage ##

TODO Stuff: sort stuff

## Find Boundaries
The script [`find_boundaries`](./find_boundaries.m) finds the location of the pycnocline and/or ice. In the case of ice, it can be used to mask that part of the flow. In both cases, a brightness threshold is used, although this is much more fiddly to find for the ice. 

## Calculate APE
The script [`dfi_ts2APE`](./dfi_ts2APE.m) calculates APE density (J/m) from a column timeseries image using the method of Boegman et al. (2005). JFM, 531, 159-180. .
Pycnocline location calculated from find_boundaries, pycnocline location tool, then smoothed. 
![Pycnocline detection](./dfi_ts2APE.png)

APE calculated as:

$$APE = cg * \Delta\rho * \int_{t_0}^{t_1}(\eta^2) dt$$

## Dissipation
### Direct
We can calculate dissipation directly, using the Direct estimate version of [Doron et al. (2001)](https://doi.org/10.1175/1520-0485(2001)0312108:TCADEI%3E2.0.CO;2), assuming isotropic turbulence:

$$    \epsilon_d = 3 \nu \Biggl[ \Biggl< \left(\dfrac{\delta u}{\delta x}\right) ^2 \Biggr>  +  \Bigg \langle \left(\dfrac{\delta w}{\delta z}\right) ^2 \Bigg \rangle  + \Bigg \langle \left(\dfrac{\delta u}{\delta z}\right) ^2 \Bigg \rangle + \Bigg \langle \left(\dfrac{\delta w}{\delta x}\right) ^2 \Bigg \rangle + \\
    2 \Bigg \langle \left(\dfrac{\delta u}{\delta z} \dfrac{\delta w}{\delta x}\right) \Bigg \rangle + 4/3 \Bigg \langle \left(\dfrac{\delta u}{\delta x} \dfrac{\delta w}{\delta z}\right) \Bigg \rangle \Biggr]
$$

### Spectral
Various forms of these exist in 1D and 2D. Essentially, these calculate the power spectra of the velocity data [`spec2_ps_nopad.m`](./spec2_ps_nopad.m)/ [`spec_ps_nopad.m`](./spec_ps_nopad.m) and calculates the value of $\epsilon$ which closest fits the observed spectrum:
![Nasmyth Spectra](./NasmythSpectra.png)


## Horizontally averaged dissipation profiles ##

Horizontal averaging of TKE dissipation is performed in `plot_dissipation_profile.m` as follows:

- Data from within 0.1m from each side of the plotted FOV (figure 5 in article) was discarded to remove potential edge effects from the PIV algorithm.  To visualize the generally low sensitivity of the results to the amount of edge cut off, produce figure 1111 with various different values of xr_min and xr_max (the amount of edge trimmed at the lower and upper ends of the figures).

- At each depth level, all remaining data were averaged to provide a single average dissipaiton \overbar{\epsilon} at that level.

