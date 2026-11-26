# plot_timeslices

Plots a panelled figure comparing lab (DigiFlow `.dfi`) and/or SPINS model output
across one or more time points, with one tiled column per data source.

```matlab
ax = plot_timeslices(ColumnProperties, c_range, labelled_columns, ylims)
```

## Inputs

| Argument | Required | Description |
|---|---|---|
| `ColumnProperties` | yes | Cell array, one struct per column. See schema below. |
| `c_range` | no (default `6*[-1 1]`) | 2-element colour axis limits, applied to every panel. |
| `labelled_columns` | no (default: last column only) | Indices of columns that get a colourbar. |
| `ylims` | no (default: `'tight'` per axis) | 2-element y-axis limits applied to every panel. |

## `ColumnProperties{jj}` schema

Each cell is a struct describing one column of the figure:

| Field | Used when | Description |
|---|---|---|
| `DirectoryName` | always | Folder to `cd` into for this column (cell array of folders, for `'lab'`). |
| `Type` | always | One of `'model'`, `'model_diff'`, `'lab'`. |
| `Times` | always | Vector of output numbers / frame numbers to plot for this column. **Must be the same length for every column** — `n_times` is taken from column 1 and reused for all columns. |
| `Field` | always | Which field to plot, e.g. `'rho'`, `'u'`, `'vorty'`, `'image'`. |
| `xlimits` | always | 2-element x-axis limits for this column. |
| `FGType` | optional | Foreground overlay: `'quiver'` (not yet functional, see review) or a field name to contour over the top. |
| `contourlevels` | if `FGType` set and not `'quiver'` | Vector of contour levels to draw. |
| `ReferenceDirectory`, `RefTimes` | `Type == 'model_diff'` only | Baseline directory/times to difference against. |

    
## Output

    
`ax` — cell array of axes handles, one cell per column, each containing the
 per-time axes handles for that column (`s` internally), for further
        customisation after plotting.
   
## Dependencies

`newbluewhitered`, `cmocean`, `figure_print_format`, `spinsgrid2d`,
        `spins_params`, `subplot_labels`, `find_contour`, `get_spins_data`
        plus SPINSmatlab and DigiFlow readers
        (`get_spins_data`, `dfiread`, `dfi_grid_read`) on the MATLAB path.

## Known issues (open)

- Multi-camera lab support and the quiver overlay are stubbed out / unreachable.
- Directory is `cd`'d into per-column but never restored afterwards.
- Trailing dead code block (slope outline, time-stamp label, grid cropping) still to be cleaned up or reimplemented.

## Example

See `example_plot_timeslices.m` for a worked example constructing
`ColumnProperties` for a two-column (model vs lab) comparison figure.