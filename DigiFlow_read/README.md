**What the project does** 
MATLAB based package/workflow for processing output of DigiFlow

**Why the project is useful**
Allows the user to carry out more advanced plot production/data analysis/combination with other data streams in a familiar language

**How to get started with the project**
Primary I/O scripts are labelled df#read.m or similar. Main script is the dfireadvel for reading in UNCOMPACTED, UNCOMPRESSED .dfi files into a structure im, containing the data, grid info and metadata. 

To make plots from any kind of dfi output (a timeseries, PIV, or raw image) use the plot_dfi.m function. 

dfi_image.m is an overly complex function to call plot_dfi and format those images. 

/colormaps contains all of DigiFlow's inbuilt colormaps, plus newbluewhitered (perceptually uniform) and bluewhitered (non-perceptually uniform) - which are similar to cmocean('delta');

combine_digiflow_avi might be helpful for you to combine videos from the two cameras, but its not the best function, check out what you need (e.g. a dfi image with WCS set). In another script of mine I have something that combines directly the .dfi images and plots that, at some stage I'll write it as its own function and share (which would be far easier on many levels).

_To read in .dfi images to MATLAB_
- To read in only the data run im = dfireadvel('filename'); [x, y, dx, dy, nx, ny, x0, y0, xi, yi] = dfi_grid_read(im);  
- To Make a simple plot run plot_dfi;
- To make more complex plots or movies, run dfi_image.m with all parameters set accordingly

_To read in _

Where users can get help with your project
Who maintains and contributes to the project

**Users will also need**
subaxis.m package from https://www.mathworks.com/matlabcentral/fileexchange/3696-subaxis-subplot


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

