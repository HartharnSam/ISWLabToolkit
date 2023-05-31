# DigiFlow Read MATLAB Package

## What the project does 
MATLAB based package/workflow for processing output of DigiFlow

## Why the project is useful
Allows the user to carry out more advanced plot production/data analysis/combination with other data streams in a familiar language. Increasingly, more complex analysis codes are included as directories within this package, currently including the [`energetics`](./energetics) and [`particletracking`](./particletracking) subpackages. 

## How to get started with the project
Primary I/O (data reading in) scripts are labelled `df#read.m` or similar. Main script is the [`dfireadvel`](./dfireadvel.m) for reading in UNCOMPACTED, UNCOMPRESSED .dfi (**D**igi**F**low **I**mage) files into a structure `im`, containing the data, grid info and metadata. [`dfi_grid_read.m`](./dfi_grid_read.m) additionally reads this grid in all valuable forms into another structure (x, y vectors, matrices, Nx, Ny, dx, dy, xlim, ylim). 

### To read in .dfi images to MATLAB 
- To read in only the data run `im = dfireadvel('filename'); [x, y, dx, dy, nx, ny, x0, y0, xi, yi] = dfi_grid_read(im);`
- To Make a simple plot run plot_dfi;
- To make more complex plots or movies, run dfi_image.m with all parameters set accordingly

### To make initial plots
To make plots from any kind of dfi output (a timeseries, PIV, or raw image) use the [`plot_dfi.m`](./plot_dfi.m) function.
A very basic starting code could also be:
```
im = dfireadvel('filename'); [x, y, dx, dy, nx, ny, x0, y0, xi, yi] = dfi_grid_read(im);
imagesc(x, y, im.cdata(:, :, 1)); singlecycle; 
``` 
`dfi_image.m` is a now-redundant function to call `plot_dfi` and format those images, but it does help with example scripts. 

For sequential reading in of files: 
```
im = dfireadvel('output_0000.dfi'); [x, y, dx, dy, nx, ny, x0, y0, xi, yi] = dfi_grid_read(im);
for ii = 1:100
	fnm_in = strrep('output_####.dfi', '####', sprintf('%04d', ii));
	im = dfireadvel(fnm_in); 
	imagesc(x, z, im.cdata(:, :, 1));
end
``` 
is useful. 

## Contents
### I/O (Data reading) functions
-  [`dfireadvel.m`](./dfireadvel.m) is for reading .dfi (**D**igi**F**low **I**mage files). It is the only read script widely used and tested. 
-  [`dfi_grid_read.m`](./dfi_grid_read.m) is for turning the output of dfireadvel into a grid structure. Well used and tested.
-  [`dfdread.m`](./dfdread.m) is for reading .dfd (**D**igi**F**low **D**rawing files), but it has barely been tested, and was designed for the dfd files output from Particle Tracking Velocimetry. 
- [`dfsread.m`](./dfsread.m) is for the .dfs (**D**igi**F**low **S**tatus files), but is even less tested, and I don't think it works. 
- [`ptv2mat.m`](./ptv2mat.m) is for reading in the output of Particle Tracking Velocimetry, via [`dtf2txt.dfc`](../Cameras/dft2txt.dfc) and outputting as .mat

### Combining Scripts
Consider there are three ways to stitch images from different cameras together:
- [`combine_dfi_images.m`](./combine_dfi_images.m) sequentially reads in .dfi images and plots them over each other, such that the final image has the most recently plotted image on top in case of overlap. Can handle videos \& saving. Computationally less expensive, but graphically intensive due to overlap. 
- [`combine_im.m`](./combine_im.m) loads all the images and stitches them together using, and as a result is smoother (it incorporates all data in case of overlap), computationally quite expensive, but graphically not too bad. Output is a data matrix, which can be re-analysed in the same was as the output of dfireadvel. 
- [`combine_digiflow_avi.m`](./combine_digiflow_avi.m) the "easiest" fix to combine images, output each as .avi, and then re-stitch together. It kind of works, is quick to run, but it's not a "nice" solution. 
- [`convert_video.m`](./convert_video.m) - doesn't truly fit in this category, but a useful tool for data format conversion in MATLAB, slower than an online thing, but no risk of dodgy downloads, or limits. It doesn't work for sound (yet). 

### Plotting
- [`plot_dfi`](./plot_dfi.m) is the baseline plotter, it "should" fairly easily produce some nice looking plots. 
- [`dfi_image.m`](./dfi_image.m) is a now-redundant function to call `plot_dfi` and format those images, but it does help with example scripts. 
- [`mstreamline`](./mstreamline.m) produces streamlines from the velocity field, works in a slightly different way to `streamline`, which is sometimes beneficial to users. See plot_dfi.m for example usage. 
- [`vekplot2.m`](./veckplot2.m) plots quiver arrows, but unlike the matlab quiver function, they are scaled properly!
- [`vekLeg.m`](./vekLeg.m) produces the legend (scale arrow) for vekplot2.m

### Subpackages
- [`/colormaps`](./colormaps) contains all of DigiFlow's inbuilt colormaps, plus newbluewhitered (perceptually uniform) and bluewhitered (non-perceptually uniform) - which are similar to the [cmocean](https://uk.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps) 'delta' colormap. 
- [`energetics`](./energetics/) contains various tools for analysing energetics measures in the experimental output. Primarily the dissipation calculation code from Peter Sutherland, and the APE estimation code. Readme within, but future work should include introducing a Thorpe Scale code, and wider discussions around energetics in the lab. 
- [`particletracking`](./particletracking) codes for processing particle tracking, see individual readme for this. 




## Users will also need
- [subaxis.m package](https://www.mathworks.com/matlabcentral/fileexchange/3696-subaxis-subplot) - NOTE: this is in the extras part of the ISWLabToolkit.



