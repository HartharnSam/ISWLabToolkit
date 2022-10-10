# DigiFlow Read MATLAB Package

## What the project does 
MATLAB based package/workflow for processing output of DigiFlow

## Why the project is useful
Allows the user to carry out more advanced plot production/data analysis/combination with other data streams in a familiar language

## How to get started with the project
Primary I/O (data reading in) scripts are labelled `df#read.m` or similar. Main script is the `dfireadvel` for reading in UNCOMPACTED, UNCOMPRESSED .dfi (**D**igi**F**low **I**mage) files into a structure `im`, containing the data, grid info and metadata. 

To make plots from any kind of dfi output (a timeseries, PIV, or raw image) use the `plot_dfi.m` function. 

`dfi_image.m` is an overly complex function to call `plot_dfi` and format those images. 

`/colormaps` contains all of DigiFlow's inbuilt colormaps, plus newbluewhitered (perceptually uniform) and bluewhitered (non-perceptually uniform) - which are similar to the [cmocean](https://uk.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps) 'delta' colormap. 

`combine_digiflow_avi` might be helpful for you to combine videos from the two cameras into an .avi movie. But its not the best function, check out what you need (e.g. a dfi image with WCS set). `combine_digiflow_images.m` is not perfect, but it does a decent job of directly combining .dfi images and plots together. 

### To read in .dfi images to MATLAB 
- To read in only the data run `im = dfireadvel('filename'); [x, y, dx, dy, nx, ny, x0, y0, xi, yi] = dfi_grid_read(im);`
- To Make a simple plot run plot_dfi;
- To make more complex plots or movies, run dfi_image.m with all parameters set accordingly

## Help, Authorship & Licence
contact current maintainer via github.io/HartharnSam , email: s.hartharn-evans2@ncl.ac.uk , twitter @HartharnSam

Project created by Sam Hartharn-Evans, but most code is contributed from other authors, including J.K.Sveen@damtp.cam.ac.uk (digiflow_read), Peter Sutherland (dissipation toolbox), Marek Stastna, Magda Carr, and the DJLES toolbox is Dunphy et al. (2011).

This package is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This package is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

See <https://www.gnu.org/licenses/>.

## Users will also need
- [subaxis.m package](https://www.mathworks.com/matlabcentral/fileexchange/3696-subaxis-subplot)

\/ a lot of this may be out of date, but left for any value it has in lieu of the time to update it

Files:
    PROGRAMMES: 
	- PIV_image.m :         [SCRIPT] Makes an image from single .dfi image/frame of PIV output showing vorticity (background) and velocity (quivers or streamline)

    FUNCTIONS (called by above programmes):
	- dfireadvel : 		[FUNCTION] Reads in .dfi image to MATLAB, outputs a structure with all "tags" from the DF file
	- vekplot2.m : 		[FUNCTION] Equivalent to QUIVER but with known scale for comparisons


PIV_movie.m
    Outputs an .avi file to DigiFlow files directory and blanks off defined slope from PIV output
    Shows (by default) vorticity as background image and velocity as scaled quiver arrows
    Parameters to set:
    - DigiFlow Directory - location of digiflow .dfi files
    - Slope equation and limits
    - Files to read (general file name and number to read)
    - File output
    - Scale for colours (vorticity), x axis limits and velocity arrows
    Input .dfi must be uncompressed and uncompacted

    Calls upon dfireadvel; bluewhitered; vekplot2.m

streamline_movie.m
    Outputs an .avi file to DigiFlow files directory and blanks off defined slope from PIV output
    Shows (by default) vorticity as background image and streamlines as red contour lines
    Parameters to set:
    - DigiFlow Directory - location of digiflow .dfi files
    - Slope equation and limits
    - Files to read (general file name and number to read)
    - File output
    - Scale for colours (vorticity), x axis limits
    Input .dfi must be uncompressed and uncompacted

    Calls upon dfireadvel; bluewhitered; streamline (BUILT IN MATLAB FUNCTION)

PIV_image.m 
    Works the same as PIV_movie.m, outputs a single image instead of movie.
    Shows (by default) vorticity as background image and velocity as scaled quiver arrows
    Parameters to set:
    - DigiFlow Directory - location of digiflow .dfi files
    - Slope equation and limits
    - File to read (file name)
    - File output (file name and format code)
    - Scale for colours (vorticity), x axis limits and velocity arrows
    Input .dfi must be uncompressed and uncompacted

    Calls upon dfireadvel; bluewhitered; vekplot2.m

