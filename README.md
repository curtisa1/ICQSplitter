# ICQSplitter
https://zenodo.org/badge/latestdoi/122151399


ICQSplitter is a Python package that will take data from the International Comet Quarterly (ICQ), 
Comet OBServation database (COBS), and JPL Horizons and produce lightcurves of a comet of interest 
with adjustable corrections available.  These data are usually provided in an 80 column text file format. 
See 'input_columns_meaning.txt' for a description of what each column represents. 

The program cofirmed works on Python 3.6.3. Remember to have the file 'Schleicher_Composite_Phase_Function.txt' in your 
working directory if you wish to perform phase angle corrections.

In analyzing data from multiple observers, questions arise of justifying which data to reject without biasing results. Moreover,
instrument calibration is not readily available to amateurs nor are such calibrations without flaw.
The purpose of this is to provide a tool to quickly sort data from amateur astronomers who regularly report data to ICQ or COBS. 
These data are widely available and have been shown to be comparable to magnitudes derived from CCD or photoelectric means, especially
for observations from experienced observers. Although these amateur data are used throughout the field you often see
people neglecting to perform heliocentric distance or phase angle corrections to the magnitude, and you 
never see people provide any statistical consideration on how to add observations from multiple observers
(see https://dps2017-aas.ipostersessions.com/default.aspx?s=3A-C3-A6-3A-8A-F3-36-EF-AE-42-79-03-2B-BE-1D-6D ). Although the 
latter corrections are not as important as the former two, this tool will allow a researcher to quickly sort or apply only the 
corrections they want to apply to any data given in ICQ or COBS format.

At its base level this program will read in the ICQ or COBS 80 column format (available from ICQ or COBS) and convert it to a .csv file 
that is more accessible to most people. As these data are from citizen astronomers and ICQ and COBS reports
all data reported from an observer, problematic entries will exist in the data. This program filters from the data
entries that do not meet a field standard set of criterion  (such as removing observations made under reported 
poor weather conditions, only using one observation per observer per night, removing 
observation that were made with telescopes when the comet was too bright, etc...). A list of all criterion for why
a observation is 'kept' or 'removed' can be found on 'reasons_data_were_removed.txt'. This code is set up such that
if there is a reason included here for why a point is removed that you do not agree with then you can comment out that section
in main().

There exists command line arguments --heliocentric and --phase that will pull ephemerides from JPL HORIZONS with
Michael Mommert's CALLHORIZONS package as well as Dave Schleicher's Composite Dust Phase Function for Comets available 
here: http://asteroid.lowell.edu/comet/dustphase.html that will perform heliocentric corrections, phase angle corrections,
or both to the kept points.

The command line argument --heliocentric will calculate a heliocentric corrected magnitude mhelio = mapp -5log(Δ). 
Where mapp is the apparent magnitude reported by the observer, delta is the distance from the observer to the target in au taken 
from JPL HORIZONS via Michael Mommert's CALLHORIZONS package written in this code in queryJPL(). 

The function queryJPL() will query JPL HORIZONS and pull the ephemerides of the object inputted at the beginning of the code.
The epoch range will be from the first time in the 'kept' observations (i.e., points remaining after the previous sorting) 
to the last time at time increments inputted by the user (default is every 30 minutes).  That is if your first date is 1996:01:19 00:00, 
final date is 1996:01:19 01:00 and your increment size is every 30 minutes then it will query JPL Horizons for the 
ephemerides of your object at 1996:01:19 00:00, 1996:01:19 00:30, and 1996:01:19 01:00 It then stores the required 
information from this query in global lists to be used later

There were two choices for this block, query JPL HORIZONS at each point in the data or query once over the entire date/time range.
I went with the latter as each individual query to JPL HORIZONS takes quite a bit of time, although this way
does require much more sorting later on. Additionally, for date ranges over long periods of time you can only
query JPL at 30 minute increments before running into their 100,000 epoch query limit. However, the change
in delta and phase angle during a 30 minute (or even 1 hour) period is insignificant to the fact that amateurs report
these magnitudes to one decimal place.

The command line argument --phase will calculate a phse angle corrected magnitude mph = mapp+2.5*log(Φ(ϴ)) where mapp
is the raw magnitude and Φ(ϴ) is Schleicher's composite phase function.

--heliocentric and --phase can be combined such that mhelio = mapp -5log(Δ) is calculated followed by mph = mhelio+2.5*loglog(Φ(ϴ)). All magnitudes along with heliocentric distance, r, are reported in the 'keepers.csv' file.

curtisa1 (at) mail.usf.edu, latest version: v1.0, 2018-19-19

*	v1.0: Sorts problematic entries from data, performs heliocentric and phase corrections.

Sources of error (albeit miniscule) which exist in this code include:
-Each observation in the data is only compared to JPL HORIZONS to the nearest half-hour, so the extracted phase 
 angles and deltas are only good to the nearest half hour.
-Schleicher provides his phase function to whole number degrees. So for instance if the ephemerides says the object had
 a phase angle of 30.67 at a certain date / time, then the composite phase function normalized to 0 degrees (column 2 in 
 'Schleicher_Composite_Phase_Function.txt') that is used in calculations would be 0.4142. Ideally this could be extrapolated so 
 that no rounding is performed for the phase angle, but that is deemed not necessary at thistime.
 -The observatory code used in JPL queries is 500, the center of the Earth, as it is unknown where the observer is just
  based off of their name. If that information is available then we could do a dynamic search based off of each observer's
  individual location.

Although error exists, it should be noted that not performing heliocentric distance or phase angle corrections
causes much more error than what is present in this code. As a comet's geocentric distance only varies by about 0.001 - 0.01 au over
a thirty minute interval and phase angle by even smaller amounts, this error is deemed reasonable seeing as these 
apparent magnitudes are only reported to the first decimal place by observers anyway.

TODO:
-Add option to graph r vs mapp, r vs mhelio, and r vs mph.
-Remove necessity of having 'Schleicher_Composite_Phase_Function.txt' in working directory.
-Add ability to increment the JPLHorizon query in other intervals besides 1-60 minutes.
-Add checks for all reference star catalogs that have special requirements as per ICQ's website.
-Add statistical analysis techniques to account for observer differences.
