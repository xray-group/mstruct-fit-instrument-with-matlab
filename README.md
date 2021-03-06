# mstruct-fit-instrument-with-matlab

## Download

Download [mstruct-fit-instrument-with-matlab-0.1-beta.zip](https://github.com/xray-group/mstruct-fit-instrument-with-matlab/archive/v0.1-beta.tar.gz)

### Requirements:

- measured data in 'xy' format
- table with positions and intensities of diffraction lines (exported e.g. from PDF2)

#### note:
- adjustments to other data formats are straightforward
- the hkl table can be generated by any software, a screenshot below shows export from Pearson's Crystal Data

![Pearson's crystal data - hkl table export](/doc/PCD-export-LaB6-table.png)

## Outline

The workflow to receive the instrumental function parameters for MStruct is

1. semi-automatic fitting of pattern segments in Matlab
2. import fitting result into an Excel worksheet
3. review the instrument polynomial approximation in Excel

## Step 1: Matlab fitting

This is a semi-automatic procedure. You may want to adjust x-ray wavelength and spectral components ratio,
input patter file name or the whole data import procedure. Two methods for importing the standard peak table
are included. The first one could be used for the hkl table exported from the old version of ICDD PDF2
database. Authors have no idea how to achive that with PDF4, however any similar table is acceptable,
so an input for a table exported from Pearson's Crystal Data is shown. The hkl table should contain
reflections x-positions in the first column and some sort of relative intensities in the third column.
If wavelength or standard material is changed you will need to adjust bounds for pattern segments
where diffraction pattern background can be considerd semi-linear.

When you start working in Matlab you may need to add directory with some helper fitting and data
import functions to Matlab path.

```Matlab
cd mstruct-fit-instrument-with-matlab
addpath('helper_functions/')
```

As indicated, the procedure is semi-automatic. The Matlab script is paused after fitting each segment
leaving the user an opportunity to inspect the fit in detail. One needs to press 'space' to contine.
The script creates simple plots with the peak profile parameters in the end.

Profile parameters are saved to **params_1.txt**.

## Step 2: Excel import

Results from `params_1.txt` must be imported into the **lab6bb-pixcel-ns.xls** template excel workbook.

In Excel use:
```
File -> Open -> params_1.txt -> choose Deliminated -> choose Space -> ...
```
and delete columns except: 2Theta, intensity, hwhm, k, Asym. The format must be compatible with the
template table in the 'params_1' worksheet.

One can see the parmaeters interpolation in the next Excel worksheet: mstruct. If there is more or
less peaks than in the template workbook the bounds for Excel functions must be extended or rows deleted.

One can see also the Cohen-Wagner plot for displacement correction in the last worksheet. Note that this
lattice parameter evaluation method is valid only for cubic materials and Bragg-Brentano geometry.
A small discrepancy in the reulting lattice parameter in the example template indicates a possible
inzorrect setting for 2Theta-Zero. However this worksheet is unralated to the instrumental profile
determination.

## Step 3: Review the instrument polynomial approximation in Excel

The mstruct Excel worksheet gives directly the profile coefficients that can be used in `MStruct`
program. There is a possibility to compare with an additional set of parameters that can be set
in the worksheet. The secondary parameters configuration in the example worksheet comes from the
whole powder pattern fitting of the given experimental pattern in MStruct. One can see the results
are slightly different. Anyway one can set e.g. older parameters set for comparison etc.

The determined profile parameters should be also critically reviewd. A limitted angular range
is measured for the instrumental function determination (20-150 deg here) but often the same
instrumental function is used later at lower or higher angles. So it is extrapolated. It is
important the profile parameters are valid in the whole x-range required. In particular
check if the profile shape parameter `k <= 1` at high angles.
