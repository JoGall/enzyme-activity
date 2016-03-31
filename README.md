#Quantify the maximal rate of enzymatic reaction, Vmax, from kinematic optical density data

For enzymatic reactions following Michaelis–Menten kinetics, the maximal rate of activity (Vmax) can be measured as the slope of the reaction curve during the linear phase of the reaction, when the enzyme becomes saturated with substrate.

These [R](https://www.r-project.org/) scripts have been designed to quantify phenoloxidase / pro-phenoloxidase activity in spectrophotometric readings taken from insect haemolymph (e.g. [Moret & Siva-Jothy, 2003](rspb.royalsocietypublishing.org/content/270/1532/2475)), but can also be used for other turbidimetric assays such as those quantifying lysozyme-like activity.

The script 'VmaxFunctions.R' contains functions to calculate Vmax from a dataframe of optical density (OD) values. In both methods, OD values can be subsetted within set time thresholds and smoothed using a running median with adjustable window size, `smoothing`. `calcVmaxRegression()` calculates Vmax as a simple linear regression through the smoothed OD data, and `plotVmaxRegression()` visualises the smoothed data and fitted Vmax by returning a plot for each sample in a pdf. `calcVmaxIntervals()` bootstraps the smoothed OD data into overlapping windows of size `windowSize`, and finds the window with the greatest rate of change (as measured by linear regression); `plotVmaxIntervals()` prints a pdf to visualise this process.

The script 'SoftmaxToR9.R' transforms a [SoftMax®](http://www.moleculardevices.com/systems/microplate-readers/softmax-pro-data-acquisition-and-analysis-software) (<=v5.1) exported text file into a dataframe using the function, `SofToR()`. The directory contains a sample SoftMax® output file, 'example.text'.
