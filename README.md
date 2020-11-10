# ICQSplitter

**Readme**

The International Comet Quarterly Splitter (ICQSplitter) is a Python based open-source software which will take data from the ICQ, Comet OBServation Database (COBS), and JPL HORIZONS to produce lightcurves of a specified target. The pipeline can be run on Unix or Windows-based operating systems. ICQSplitter was used in this text to produce lightcurves of visual magnitude data from amateur astronomers, but it is capable of taking in any measurements, including those from charge-coupled devices (CCD) that are reported in ICQ's standard 80-column format. The user has the options to apply any combination of corrections discussed in the main body of this paper. For example, a user with observational magnitudes from a relatively non-dusty comet may wish to forgo the application of a phase correction.

At its base level (i.e., without any command line arguments), this program will read in the ICQ or COBS 80 column format (available from ICQ or COBS) and convert it to a .csv file that is more accessible to most people. As these data are from citizen astronomers and ICQ and COBS reports all data reported from an observer, problematic entries will exist in the data. This program filters from the data entries that do not meet a field standard set of criteria (such as removing observations made under reported poor weather conditions, only using one observation per observer per night, removing observations that were made with telescopes when the comet was too bright, etc...). A list of all criterion for why an observation is 'kept' or 'removed' can be found on 'reasons_data_were_removed.txt'. This code is set up such that if there is a reason included here for why a point is removed that you do not agree with then you can comment out that section in main().

This document describes the functionality of ICQSplitter Version 3.0 as of 28 January 2020. Also refer to the documentation for installation guides and additional support.

**Citing this Repository**

As per AAS's style guidelines (availble here: http://journals.aas.org/authors/references.html), please cite the DOI for the most recent software release as:

Olivia Curtis, David Rabson, Nathan Lastra, Sharlene Cruz-Gonzalez, & Maria Womack. (2020, January 26). ICQSplitter: Python Tool for Handling Visual or CCD Magnitudes of Comets (Version v3.0). Zenodo. http://doi.org/10.5281/zenodo.3628044

**1.1 Methods and Implementation**

This software is a standalone Python 3.6.4 script. It makes use of Python packages that are freely available and easy to install through the Python Package Index; required packages include NumPy, SciPy, matplotlib, and CALLHORIZONS. At its base level (i.e., without any optional command line arguments input) this program will read in ICQ or COBS 80 column format and filter out data that do not meet the criteria discussed in \ref{sec:removed}. The data must be from a a single small body comprised of either entirely CCD or visual magnitude data and comprised of observations from a single orbit around the sun spanning a date range no larger than five years.

**1.2 ICQSplitter Arguments**

ICQSplitter has a few optional command line arguments. Refer to the documentation provided online for more details on how to use them.

**1.2.1 --heliocentric**

This command applies heliocentric corrections to the raw magnitudes. In doing so, ICQSplitter will use the CALLHORIZONS package to query JPL HORIZONS to extract the heliocentric distance, geocentric distance, and phase angle of the target. This function will perform a single query of JPL HORIZONS over the range of dates provided in increments inputted by the user. The default time interval is 30 minutes increments. For instance,  if your first date is 1996:01:19 00:00, final date is 1996:01:19 01:00, and your increment size is every 30 minutes then it will query JPL Horizons for the ephemerides of your object at 1996:01:19 00:00, 1996:01:19 00:30, and 1996:01:19 01:00. As JPL only allows users to pull 100,000 in a single query, if one wishes to analyze all of their data in a single compilation of ICQSplitter then small time increments should be avoided. ICQSplitter uses the ephermides data to perform heliocentric corrections. 

**1.2.2 --phase**

Applies phase corrections to the given magnitudes. If --heliocentric and --phase are called at the same time then the phase angle corrections will be applied onto the heliocentric corrected magnitudes, else JPL will be queried for the first time and phase corrections will be applied to the raw magnitudes. The phase angles are cross referenced to Dave Schleicher's Composite Dust Function for Comets.

**1.2.3 --stats**

Performs the statistical analysis. The program will automatically split any dataset into pre- and post-perihelion and perform the statistics on each set separately. ICQSplitter follows procedures for regression analysis through the methods of singular value decomposition using NumPy's Linear Algebra package. After a polynomial fit has been taken to convergence, Python's Statistics package is used to perform the Students t and probability tests on each observer's data. If one observer is found to fail the stationarity test in either epoch, then that observer is removed from the dataset and the procedure is repeated. The --stats command is always issued after --heliocentric and --phase (if those commands have also been given). 

**1.2.4 --plot**

This command will produce plots with Pythons matplotlib package. For instance, if a user runs her code on data with the --phase and --stats commands then --plots will produce individual graphs of m_{tot}, m_{helio}, m_{phase}, and m_{shift}. 
