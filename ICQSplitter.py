"""



ICQSplitter is a Python package that will take data from the International Comet Quarterly (ICQ), Comet OBServation database (COBS), and
JPL Horizons and produce lightcurves of the comet with adjustable corrections. These data are usually provided in an 80
column text file format. See 'input_columns_meaning.txt' for a description of what each column represents.

At its base level this program will read in this 80 column format (available from ICQ or COBS) and convert it to a .csv file 
that is more accessible to most people. As these data are from citizen astronomers and ICQ and COBS reports
all data reported from an observer, this program also filters from the data entries that do not meet a field standard set of criterion 
(such as removing observations made under reported poor weather conditions, only using one observation per observer per night, removing 
observation that were made with telescopes when the comet was too bright, etc...). A list of all criterion for why
a observation is 'kept' or 'removed' can be found on 'reasons_data_were_removed.txt'. This code is set up such that
if there is a reason included here for why a point is removed that you do not agree with then you can comment out that section
in main().

There exists command line arguments --heliocentric and --phase that will pull ephemerides from JPL HORIZONS with
Michael Mommert's CALLHORIZONS package as well as Dave Schleicher's Composite Dust Phase Function for Comets available 
here: http://asteroid.lowell.edu/comet/dustphase.html that will perform heliocentric corrections, phase angle corrections,
or both to the kept points.

Users may also perform statistical corrections with the --stats command line argument. Please see the file
'statistics_method_appendix.txt' for an understanding of what this command does. (note at least one of --heliocentric
and --phase must be used along with --stats to query JPL). WARNING: --stats should only be used for one
dates corresponding to one revolution around the sun for that object (i.e., if a comet has data with perihelions occuring on dates
x, y, and z then this function will work for dates between x and z exclusive, separating it into pre-perihelion data in range (x,y] and
post-perihelion data in range (y,z) ).

The command line argument --plot will also plot any available data (i.e., any combination of raw magnitudes, mehlio, mphase, and 
mshift Vs. heliocentric distance). (note at least one of --heliocentric and --phase must be used along with --plot to query JPL).



curtisa1 (at) mail.usf.edu, latest version: v1.0, 2018-04-04

Available command line arguments (type these into terminal when compiling program):
--heliocentric
--phase
--stats
--plot

*	v1.0: Sorts problematic entries from data, performs heliocentric distance and phase angle corrections.
*	v1.1: Added Input Argument CCD_Bool for people using only CCD Measurements.
*	v2.0: Added statistical correction and plotting command line arguments!
*   v2.1: Fixed issues with statistical analysis. Added option to get full detailed stats analysis. Uncomment 509 - 515 and 578 - 584 and 650 - 653 to see the output files!.


"""

###############################
####### Input Arguments #######
###############################

input_file = 'GreeneWithBiver.txt'			#Name of your input file
target_nickname = 'HB'						#Nickname of target for output file organization (for example, HB = Hale-Bopp)
small_body_designation = '902014;'			#Name of your small body ex) 'ceres' or 'eris'
JPL_Time_Increment = 30 					#How much to increment JPL queries in minutes up to 60.
ouput_file_kept_points = 'keepers.csv'		#Name of output file for points that meet all sorting criterion
output_file_rejected_points = 'removed.csv'	#Name of output file for points that were removed from the data
perihelion = '1997/04/01'					#Datetime of perihelion format YYYY/MM/DD
CCD_Bool = 1								#If 0 then user only has CCD measurements, if 1 then user has visual magnitude measurements

###############################
####### Input Arguments #######
###############################

try:
	import callhorizons
except:
	print('Please install the callhorizons python package with pip install callhorizons')
import os
import numpy as np
import math
import csv
import sys
import time
import matplotlib
import pylab as plt
import statistics
from scipy import stats
from scipy import linalg
from matplotlib.ticker import MultipleLocator
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

#Deletes the i-th observation in the dataset after it failed a sorting criterion
#Defined later, metalist is a list whose elements are each column (as lists) of the input file
#see mainblock to know what each element of metalist is (e.g., metalist[23] == list of observer for each observation)
#read as: for each colum in our data, start at the last point, work backwards and delete the specified point
def deletearow(i):
	for x in range (len(metalist)-1, -1, -1):
		del metalist[x][i]

#Adds the i-th deleted observation to the 'rejected.csv' file with y as a reason it was removed
#removed_metalist is the same as metalist except for the rejected points during sorting
def addToremoved(i,y):
	for x in range (len(metalist)-1, -1, -1):
		removed_metalist[x].append(metalist[x][i])
	removed_metalist[24].append("REMOVED POINT")
	removed_metalist[25].append(list_of_reasons_removed[y])

#Removes a point from the dataset if there is more than one observation from the same observer at the same night
#read as: for total number of points in a column (i.e., total number of observations), work backwards and delete necessary points
#If observer, year, month, and day are the same then delete a point based on the focal lengths metalist[11] or observation method metalist[7]
def deleteForDuplicatedDates(numdelforduplicatedate):
	for k in range (len(metalist[2])-1, -1, -1):
		if k == 0:
			print("number deleted for duplicated observation dates by same person "+str(numdelforduplicatedate))
			break
		if (metalist[23][k]==metalist[23][k-1]) and (metalist[3][k] == metalist[3][k-1]) and (metalist[4][k] == metalist[4][k-1]) and (math.floor(float(metalist[5][k])) == math.floor(float(metalist[5][k-1]))):
			if (float(metalist[11][k])<float(metalist[11][k-1])):
				numdelforduplicatedate = numdelforduplicatedate +1
				addToremoved(k-1,reasonForDelete)
				deletearow(k-1)
				continue
			elif (metalist[7][k] == "S"):
				numdelforduplicatedate = numdelforduplicatedate +1
				addToremoved(k-1,reasonForDelete)
				deletearow(k-1)
				continue
			elif (metalist[7][k-1] == "S"):
				numdelforduplicatedate = numdelforduplicatedate +1
				addToremoved(k,reasonForDelete)
				deletearow(k)
				continue
			elif (metalist[7][k] == "M"):
				numdelforduplicatedate = numdelforduplicatedate +1
				addToremoved(k-1,reasonForDelete)
				deletearow(k-1)
				continue
			elif (metalist[7][k-1] == "M"):
				numdelforduplicatedate = numdelforduplicatedate +1
				addToremoved(k,reasonForDelete)
				deletearow(k)
				continue
			elif ((metalist[7][k] == "B") or (metalist[7][k] == "I")):
				numdelforduplicatedate = numdelforduplicatedate +1
				addToremoved(k-1,reasonForDelete)
				deletearow(k)
				continue
			elif ((metalist[7][k-1] == "B") or (metalist[7][k-1] == "I")):
				numdelforduplicatedate = numdelforduplicatedate +1
				addToremoved(k-1,reasonForDelete)
				deletearow(k-1)
				continue
			else:
				numdelforduplicatedate = numdelforduplicatedate +1
				addToremoved(k,reasonForDelete)
				deletearow(k-1)
				continue

#Removes points if no magnitude (metalist[8]) is reported
def nomagreported(numdeltefornomagreported):
	for k in range (len(metalist[2])-1,-1,-1):
		if (metalist[8][k] == ""):
			numdeltefornomagreported = numdeltefornomagreported +1
			addToremoved(k-1,reasonForDelete)
			deletearow(k)
	print("number deleted for no magnitudereported "+ str(numdeltefornomagreported))

#Removes points if the reverse binocular method is used in the magnitude methods notes columns of ICQ Data metalist[6] or metalist[22]
def reversebinocular(numdeletedforreversebinocmethod):
	for k in range (len(metalist[2])-1,-1,-1):
		if ((metalist[6][k] == 'r' ) or (metalist[22][k] == 'r')):
			numdeletedforreversebinocmethod = numdeletedforreversebinocmethod +1
			addToremoved(k-1,reasonForDelete)
			deletearow(k)
	print("number deleted for using reverse binocular method "+ str(numdeletedforreversebinocmethod))

#Removes points if a bad extinction correction is used, denoted in magnitude methods notes columns of ICQ Data metalist[6] or metalist[22]
def badExtinctionCorrection(numdeletedforbadextinctioncorrection):
	for k in range (len(metalist[2])-1,-1,-1):
		if(metalist[6][k] == '&') or (metalist[22][k] == '&'):
			numdeletedforbadextinctioncorrection = numdeletedforbadextinctioncorrection +1
			addToremoved(k-1, reasonForDelete)
			deletearow(k)
	print("number deleted for a poor extinction correction " + str(numdeletedforbadextinctioncorrection))
	
#Removes points if poor weather is reported
def poorweather(numdeletedforpoorweather):
	for k in range (len(metalist[2])-1,-1,-1):
		if (metalist[9][k] == ':'):
			numdeletedforpoorweather = numdeletedforpoorweather +1
			addToremoved(k-1,reasonForDelete)
			deletearow(k)
	print("number deleted for poor weather conditions " + str(numdeletedforpoorweather))
	
#Removes points if a telescope is used for mags brighter than 5.4 or binoculars brighter than 1.4
#Please see ICQ webpage for meaning each letter in the Instrument column of ICQ data (metalist[12])
def checkTelescopesandBinocMethods(numberremovedfortelescopeunder5_4, numberremovedforbinocularunder1_4):
	telescopemethod = ['C','R','D','I','J','L','M','q','Q','r','S','T','U','W','Y']
	binocularmethod = ['A','B','N','O']
	for k in range (len(metalist[2])-1,-1,-1):
		if (metalist[8][k] == ""):
			continue
		if ((metalist[12][k] in telescopemethod) and (float(metalist[8][k])<5.4)):
			reasonForDelete = 5
			addToremoved(k-1,reasonForDelete)
			deletearow(k)
			numberremovedfortelescopeunder5_4 = numberremovedfortelescopeunder5_4 +1
		if ((metalist[12][k] in binocularmethod) and (float(metalist[8][k])<1.4)):
			reasonForDelete = 6
			addToremoved(k-1, reasonForDelete)
			deletearow(k)
			numberremovedforbinocularunder1_4 = numberremovedforbinocularunder1_4 +1
	print("number deleted for using a telescope under m = 5.4 " + str(numberremovedfortelescopeunder5_4))
	print("number deleted for using a binocular under m = 1.4 " + str(numberremovedforbinocularunder1_4))

#Removes points for magnitude acquisition methods that are not ideal.
def removeForMagMethod(numberremovedfornotspecifiedmagmethod):
	allowedmagnitudemethod = ['S','B','M','I','E']
	#ccdlist = ['C','c','g','H','k','r','u','Y','l']
	for k in range (len(metalist[2])-1,-1,-1):
		if (metalist[7][k] not in allowedmagnitudemethod):
			addToremoved(k-1,reasonForDelete)
			deletearow(k)
			numberremovedfornotspecifiedmagmethod = numberremovedfornotspecifiedmagmethod +1
	print("number deleted for using a method not specified by green (i.e. column 27 not being S, B, M, I, or E), prioritizing S then M " + str(numberremovedfornotspecifiedmagmethod))

#Some catalogs reported in ICQ's Recommended or Condemned start magnitude catalogs have requirements on when they can be used.
#For instance, the SC catalog which appears in the test data should only be used when the comet is brighter than 8.1
#TODO: add all checks for all catalogs with such requirements.
def checkSCcatalog(numberremovedforSCcatalog):
	for k in range(len(metalist[2])-1,-1,-1):
		if (metalist[10][k] == "SC"):
			if(float(metalist[8][k]) > 8.1):
				addToremoved(k-1,reasonForDelete)
				deletearow(k)
				numberremovedforSCcatalog = numberremovedforSCcatalog +1
	print("number deleted for SC catalog being used on object dimmer than 8.1 ", numberremovedforSCcatalog)
	
#This functions will take a decimal date as reported in ICQ and convert it to YYYY:MM:DD HH:MM:SS format.
#That is, for each date in the data (metalist[3], metalist[4], and metalist[5]) it will convert it to the above format.
#For example if metalist[3][0] == 1996, metalist[4][0] == 04, and metalist[5][0] == 30.50 
#then it will append to the output list 1996:04:30 12:00:00
#These dates are then used to compare to the result of the JPL HORIZONS Ephemerides query to find the nearest time in the query
def decimaldate2hhmmss():
	global date_compare_to_JPL
	date_compare_to_JPL = []
	for i in range(0, len(metalist[3])):
		tmp_date = ""
		tmp_decimal = ""
		past_decimal = 0
		for j in range(0, len(metalist[5][i])):
			if (metalist[5][i][j] != ".") and (past_decimal == 0):
				tmp_date = tmp_date + metalist[5][i][j]
			if metalist[5][i][j] == ".":
				past_decimal = 1
			if (past_decimal ==1):
				tmp_decimal = tmp_decimal + metalist[5][i][j]
		tmp_hours = str(float(tmp_decimal) * 24)
		past_decimal = 0
		hours = ""
		tmp_decimal2 = ""
		for j in range(0, len(tmp_hours)):
			if (tmp_hours[j] != ".") and (past_decimal == 0):
				hours = hours + tmp_hours[j]
			if tmp_hours[j] == ".":
				past_decimal = 1
			if (past_decimal ==1):
				tmp_decimal2 = tmp_decimal2 + tmp_hours[j]
		tmp_minutes = str(float(tmp_decimal2) * 60)
		past_decimal = 0
		minutes = ""
		tmp_decimal3 = ""
		for j in range(0, len(tmp_minutes)):
			if (tmp_minutes[j] != ".") and (past_decimal == 0):
				minutes = minutes + tmp_minutes[j]
			if tmp_minutes[j] == ".":
				past_decimal = 1
			if (past_decimal == 1):
				tmp_decimal3 = tmp_decimal3 + tmp_minutes[j]
		tmp_seconds = str(float(tmp_decimal3) * 60)
		seconds = ""
		for j in range(0, len(tmp_seconds)):
			if (tmp_seconds[j] != "."):
				seconds = seconds + tmp_seconds[j]
			if (tmp_seconds[j] == "."):
				break
		date_compare_to_JPL.append(metalist[3][i] + "-"+metalist[4][i] + "-" + tmp_date + " " + str("%02d"%(float(hours),)) + ":"+str("%02d"%(float(minutes),))+":"+ str("%02d"%(float(seconds),)))

#This will query JPL HORIZONS and pull the ephemerides of the object inputted above.
#The epoch range will be from the first time in your 'kept' observation (i.e., points remaining after the previous sorting) to the last time
#at the increment range also inputted (default is every 30 minutes). 
#That is if your first date is 1996:01:19 00:00, final date is 1996:01:19 01:00 and your increment size is every 30 minutes then it will
#query JPL Horizons for the ephemerides of your object at 1996:01:19 00:00, 1996:01:19 00:30, and 1996:01:19 01:00
#It then stores the required information from this query in global lists to be used later
#
#
#There were two choices for this block, query JPL HORIZONS at each point in the data or query once over the entire date/time range.
#I went with the latter as each individual query to JPL HORIZONS takes quite a bit of time, although this way
#does require much more sorting later on. Additionally, for date ranges over long periods of time you can only
#query JPL at 30 minute increments before running into their 100,000 epoch query limit. However, the change
#in delta and phase angle during a 30 minute (or even 1 hour) period is insignificant to the fact that amateurs report
#these magnitudes to one decimal place.
def queryJPL():
	global OBJDelta
	global OBJDates
	global OBJPhase
	global OBJr
	initial_date = str(metalist[3][0]) + "/" + str(metalist[4][0]) + "/" + str(math.floor(float(metalist[5][0])))
	final_date = str(metalist[3][0]) + "/" + str(metalist[4][0]) + "/" + str(math.floor(float(metalist[5][0])))
	place_in_list_initial = 0
	place_in_list_final = 0
	for j in range (0, len(metalist[3])):
		check_date = str(metalist[3][j]) + "/" + str(metalist[4][j]) + "/" + str(math.floor(float(metalist[5][j])))
		newdate1 = time.strptime(initial_date, "%Y/%m/%d")
		newdate2 = time.strptime(check_date, "%Y/%m/%d")
		newdate3 = time.strptime(final_date, "%Y/%m/%d")
		if newdate2 < newdate1:
			initial_date = check_date
		elif newdate2 > newdate3:
			final_date = check_date
	initial_date = initial_date.replace("/","-") + " 00:00"
	final_date = final_date.replace("/","-") + " 23:00"
	#Queries JPL HORIZONS
	small_body = callhorizons.query(small_body_designation)
	#By default time_string = '30m' that is from initial_date to final_date at 30 minute increments
	#Change increment at top of code for increments from every 1 to 60 minutes.
	#Code currently does not allow for increments in time periods greater than 60 minutes
	#TODO fix the above, will require changing the line below as well as the if statements in the optional argument blocks that compare the times.
	time_string = str(JPL_Time_Increment) + 'm'
	small_body.set_epochrange(initial_date, final_date , time_string)
	print (small_body.get_ephemerides(500), 'epochs queried')
	OBJDelta = small_body['delta']
	OBJDates = small_body['datetime']
	OBJPhase = small_body['alpha']
	OBJr = small_body['r']
	tmp_dates =[]
	hours = ""
	minutes = ""
	seconds = ""
	date_compare_to_JPL = decimaldate2hhmmss()
	for k in range(0, len(OBJDates)):
		OBJDates[k] = OBJDates[k].replace("Jan","01")
		OBJDates[k] = OBJDates[k].replace("Feb","02")
		OBJDates[k] = OBJDates[k].replace("Mar","03")
		OBJDates[k] = OBJDates[k].replace("Apr","04")
		OBJDates[k] = OBJDates[k].replace("May","05")
		OBJDates[k] = OBJDates[k].replace("Jun","06")
		OBJDates[k] = OBJDates[k].replace("Jul","07")
		OBJDates[k] = OBJDates[k].replace("Aug","08")
		OBJDates[k] = OBJDates[k].replace("Sep","09")
		OBJDates[k] = OBJDates[k].replace("Oct","10")
		OBJDates[k] = OBJDates[k].replace("Nov","11")
		OBJDates[k] = OBJDates[k].replace("Dec","12")

#Add headers to the columns in the output csv files.
def add_headers(x):
	x[0].insert(0,'col 1-3 : short period comet designation')
	x[1].insert(0, 'col 4-9 : Standard comet designation')
	x[2].insert(0, 'col 10 : multiple nuclei present?')
	x[3].insert(0,'col 12-15 : year observed')
	x[4].insert(0,'col 17-18 : month observed')
	x[5].insert(0, 'col 20-24')
	x[6].insert(0,'col 26 : special note / extinction note')
	x[7].insert(0, 'col 27 : Magnitude collection method')
	x[8].insert(0, 'col 28-32 : visual magnitude estimate')
	x[9].insert(0, 'col 33 : poor conditions?')
	x[10].insert(0, 'col 34 - 35 : reference catalog')
	x[11].insert(0, 'col 36-40 : instrument aperture in centimeters')
	x[12].insert(0, 'col 41 : instrument type')
	x[13].insert(0, 'col 42 - 43 : focal ratio')
	x[14].insert(0, 'col 44-47 : magnification used')
	x[15].insert(0, 'col 49 : error estimate for coma diameter')
	x[16].insert(0, 'col 50 - 54 : coma diameter in arcminutes')
	x[17].insert(0, 'col 55 : special note on central condensation of comet')
	x[18].insert(0, 'col : 56 -57 : degree of condensation (note / means estimate)')
	x[19].insert(0, 'col 59 - 64 : error of tail approximation and tail approximation')
	x[20].insert(0, 'col 65 - 67 : direction tail is pointed')
	x[21].insert(0, 'col 69-74 : ICQ reference publication')
	x[22].insert(0, 'col 75 : second special note / extinction note ')
	x[23].insert(0, 'col 76-80 : observer name')
				
#Sorts r from largest to smallest while keeping track of all of ICQ metadata (listPreorPost) for each observation.
#Used twice, first in stats_shifts function (firstpass = 0) where keeping track of metadata information is important
#It is also used in plotting, (firstpass = 1) where we do not need the meta information so we skip this step if that is the case
def sortbyr(listPreorPost,r,mags, firstpass):
	sorted_metalist = []
	r_sorted = []
	mag_sorted = []
	r = np.array(r)
	r_sorted = -np.sort(-r)
	for i in range (0,len(r_sorted)):	#for each object in the sorted list
		for j in range(0,len(r)):		#find its corresponding entry in unsorted list, this position is the position of all corresponding metadata
			if firstpass == 0:
				if (r_sorted[i] == r[j]) and (listPreorPost[j] not in sorted_metalist):
					sorted_metalist.append(listPreorPost[j])
					mag_sorted.append(mags[j])
					break
			if (r_sorted[i] == r[j]) and (firstpass == 1):
				mag_sorted.append(mags[j])
				break
	return sorted_metalist, mag_sorted, r_sorted
	
def getcolumn(matrix, i):
    return [row[i] for row in matrix]
		
#Performs the procedures to iterate a polynomial to convergence of tolerance 0.0001 in given data (see file 'Statistics_method_appendix.txt' in the GitHub repository)
#Inputs: preorpost - String stating whether this is pre-perihelion or post-perihelion data (determined later in the code)
#listoflists - ICQ metadata
#dateThours - date in YYYY-MM-DDTHH:MM:SS format used for submitting to NASA PDS and checking whether pre or post perihelion 
#deltas - list of geocentric distances (from JPL)
#phases - list of phase angles (from JPL)
#helio_distances - list of r- values
#condemned_list - List of observers who have failed the stationary test, not to be used on future convergence tests
#other_mag Last calculated magnitude (either mhelio or mphase depending on which combination of the two the user used)
#first_pass - If 1 then this is the first time the data are having a polynomial fit to them (so that r-values do not have their natural log taken twice upon being read in). 
#Primary Return is mshift - the magnitudes shifted by the mean of an observer's residuals between a global polynomial fit and their data (iterated to convergence)
def stats_shifts(preorpost, listoflists, corrected_mag, dateThours, deltas, phases, helio_distances, condemned_list, other_mag, first_pass):
	mshift = []
	obs_list = []
	sorted_stats = []
	stats = []
	mags = []
	r = []
	new_poly_fit = []
	original_poly_fit = []
	mags_sorted_stat = []
	r_sorted_stat = []
	stdev_per_observer = []
	resid_per_obs = []
	stdev_resid_per_observer =[]
	count_per_observer = []
	residuals = []
	sum_resid_per_observer = []
	mean_resid_per_observer = []
	tolerance = 0.0001
	
	#Checks data are in pre-perihelion range, sorts the data by r (large to small), and rearranges all metadata accordingly
	if preorpost == 'pre':
		#Next if loop necessary in case user only supplied post-perihelion data. Read: if no data are found for pre-perihelion then
		#fill it with blank spaces that the program will know to skip later.
		if (len(other_mag) == 1) or (len(other_mag) == 0) or ('' in other_mag):
			for j in range(0, len(listoflists[0])):
				other_mag.append('')
		#for each observer, check if it is pre-perihelion, if so take log(r and sort)
		print(len(listoflists[0]))
		print(len(corrected_mag))
		for j in range (0, len(listoflists[0])):
			tmprow = []
			newdate4 = time.strptime(perihelion, "%Y/%m/%d")
			datetocheck = str(listoflists[3][j]) + "/" + str(listoflists[4][j]) + "/" + str(math.floor(float(listoflists[5][j])))
			newdate5 = time.strptime(datetocheck, "%Y/%m/%d")

			if (newdate5 <= newdate4) and (listoflists[23][j] not in condemned_list):
				for k in range (0,len(listoflists)):
					tmprow.append(listoflists[k][j])
				tmprow.append(dateThours[j])
				tmprow.append(corrected_mag[j])
				tmprow.append(deltas[j])
				tmprow.append(phases[j])
				mags.append(float(corrected_mag[j]))
				if (len(other_mag) != 1) and (len(other_mag) != 0) and ('' not in other_mag):
					tmprow.append(other_mag[j])
				else:
					other_mag.append('')
					tmprow.append(other_mag[j])
				if first_pass ==1:
					r.append(math.log10(float(helio_distances[j])))
				elif first_pass == 0:
					r.append(float(helio_distances[j]))
				stats.append(tmprow)
			
	#Reads in only pre-perihelion data from data in input function			
	if preorpost =='post':
		#Next if loop necessary in case user only supplied pre-perihelion data. Read: if no data are found for post-perihelion then
		#fill it with blank spaces that the program will know to skip later.
		if (len(other_mag) == 1) or (len(other_mag) == 0) or ('' in other_mag):
			for j in range(0, len(listoflists[0])):
				other_mag.append('')
		#for each observer, check if it is pre-perihelion, if so take log(r and sort)
		for j in range (0, len(listoflists[0])):	
			tmprow = []
			newdate4 = time.strptime(perihelion, "%Y/%m/%d")
			datetocheck = str(listoflists[3][j]) + "/" + str(listoflists[4][j]) + "/" + str(math.floor(float(listoflists[5][j])))
			newdate5 = time.strptime(datetocheck, "%Y/%m/%d")
			if (newdate5 > newdate4) and (listoflists[23][j] not in condemned_list):
				for k in range (0,len(listoflists)):
					tmprow.append(listoflists[k][j])
				tmprow.append(dateThours[j])
				tmprow.append(corrected_mag[j])
				tmprow.append(deltas[j])
				tmprow.append(phases[j])
				mags.append(float(corrected_mag[j]))
				if (len(other_mag) != 1) and (len(other_mag) != 0) and ('' not in other_mag):
					tmprow.append(other_mag[j])
				else:
					other_mag.append('')
					tmprow.append(other_mag[j])
				if first_pass == 1:
					r.append(math.log10(float(helio_distances[j])))
				elif first_pass == 0:
					r.append(float(helio_distances[j]))
				stats.append(tmprow)
			
	#stats is the "metalist" containing all of the information in the function's input arguments.
	if len(stats) != 0:
	
		sorted_stats, mags_sorted_stat, r_sorted_stat = sortbyr(stats,r,mags,0)
				
		for i in range (0, len(sorted_stats)):
			if sorted_stats[i][23].strip() not in obs_list:
				obs_list.append(sorted_stats[i][23].strip())
				
		#Beginning iterating polynomial fits to convergance
		for k in range (0, 21):
		
			#if first_pass == 1:
			#	file_name = preorpost +'.'+str(k)+'.'+'out'
			#	file_writer = csv.writer(open(file_name, 'w'), delimiter =',')
			
			#if first_pass != 1:
			#	file_name = preorpost +'.'+str(k)+'.'+'With_Dropped_Obs'
			#	file_writer = csv.writer(open(file_name, 'w'), delimiter =',')
			
			#for the first polynomial fit, all weights are set to 1.0 as we do not have standard deviations yet
			if k == 0:
				#initializing lists needed through routine
				for i in range(0,len(obs_list)):
					sum_resid_per_observer.append(0)
					stdev_per_observer.append(0)
					stdev_resid_per_observer.append(0)
					count_per_observer.append(0)
					mean_resid_per_observer.append(0)
					resid_per_obs.append(0)
				
				A = np.zeros((len(mags_sorted_stat), 6))
				b = np.zeros((len(mags_sorted_stat), 1))
				
				#Inputing A_matrix and b_vector (see Numerical Recipes)
				for i in range (0, len(mags_sorted_stat)):
					for j in range(0,6):
						A[i,j] = (r_sorted_stat[i]**j) / 1.0
					
					b[i,0] = (mags_sorted_stat[i] / 1.0)
								
				U, S, Vh = np.linalg.svd(A, full_matrices = False)
				
				S_matrix = np.zeros((6,6))
				for i in range(0, len(S)):
					S_matrix[i,i] = S[i]
									
				tmp = np.zeros((6,1))
				tmp = np.matmul(np.matmul(np.matrix.transpose(Vh), linalg.inv(S_matrix)) , np.matmul( np.matrix.transpose(U), b))

				#new_poly_fit is a vector of the coeffeficients of fifth order fit, p is correpsonding function
				tmp = np.matrix.transpose(tmp)[0]	#current polynomial fit
				new_poly_fit = []
				for i in range(len(tmp)-1, -1, -1):
					new_poly_fit.append(tmp[i])
				p = np.poly1d(new_poly_fit)
				
				print(preorpost, k, new_poly_fit)

				#calculating residuals between polynomial fit and data point
				for i in range(0,len(mags_sorted_stat)):
					residuals.append(p(r_sorted_stat[i]) - mags_sorted_stat[i])
				
				#For each observer's data, calculate the mean and stdev of their residuals
				#For each observer mshift = mags_sorted_stat + mean_resid_per_observer
				#For instance, on first iteration: mshift = mph + mean_resid_per_observer
				for o in range (0, len(obs_list)):
					tmp_list = []
					for h in range(0,len(mags_sorted_stat)):
						if sorted_stats[h][23] == obs_list[o]:
							count_per_observer[o] = count_per_observer[o] + 1
							tmp_list.append(residuals[h])
							resid_per_obs.append(residuals[h])
					mean_resid_per_observer[o] = statistics.mean(tmp_list)
				for h in range(0, len(mags_sorted_stat)):
					for o in range (0, len(obs_list)):
						if sorted_stats[h][23] == obs_list[o]:
							mshift.append(mags_sorted_stat[h] + mean_resid_per_observer[o])
							break
				print('here?')
				
				#if first_pass ==1:
				#	for h in range(0, len(mags_sorted_stat)):
				#		file_writer.writerow([str(10**(r_sorted_stat[h])), str(sorted_stats[h][23]), str(mags_sorted_stat[h]),str(mags_sorted_stat[h]), str(0), str(0), str(1)])	

				#if first_pass !=1:
				#	for h in range(0, len(mags_sorted_stat)):
				#		file_writer.writerow([str(10**(r_sorted_stat[h])), str(sorted_stats[h][23]), str(mags_sorted_stat[h]),str(mags_sorted_stat[h]), str(0), str(0), str(1)])	
							
			#For each successive iteration we repeat the same things
			#Except we use stdev_resid_per_observer accordingly in calculating A and b
			if k!=0:
				for o in range (0, len(obs_list)):
					tmp_list = []
					stdev_resid_per_observer[o] = 0
					mean_resid_per_observer[o] = 0
					sum_resid_per_observer[o] = 0
					resid_per_obs[o] = 0
					for h in range(0, len(mags_sorted_stat)):
						if sorted_stats[h][23] == obs_list[o]:
							tmp_list.append(residuals[h])
							resid_per_obs.append(residuals[h])
					stdev_resid_per_observer[o] = statistics.stdev(tmp_list)
					mean_resid_per_observer[o] = statistics.mean(tmp_list)
					
				A = np.zeros((len(mshift), 6))
				b = np.zeros((len(mshift),1))
				
				for i in range (0, len(mshift)):
					for j in range(0,6):
						for o in range (0, len(stdev_resid_per_observer)):
							if sorted_stats[i][23] == obs_list[o]:
								A[i,j] = (r_sorted_stat[i]**j) / stdev_resid_per_observer[o]
					
					for o in range(0,len(stdev_resid_per_observer)):
						if sorted_stats[i][23] == obs_list[o]:
							b[i,0] = (mshift[i] / stdev_resid_per_observer[o])
								
				U, S, Vh = np.linalg.svd(A, full_matrices = False)
				
				S_matrix = np.zeros((6,6))
				for i in range(0, len(S)):
					S_matrix[i,i] = S[i]
				
				tmp = np.zeros((6,1))
				tmp = np.matmul(np.matmul(np.matrix.transpose(Vh), linalg.inv(S_matrix)) , np.matmul( np.matrix.transpose(U), b))
						
				old_poly_fit = new_poly_fit					#previous polynomial fit
				tmp = np.matrix.transpose(tmp)[0]			#current polynomial fit
				new_poly_fit = []
				for i in range(len(tmp)-1, -1, -1):
					new_poly_fit.append(tmp[i])
				p = np.poly1d(new_poly_fit)

				#Compares coefficients between the k_th and k_th - 1 polynomial fits
				converge_test = []
				for i in range(0,len(old_poly_fit)):
					converge_test.append(old_poly_fit[i] - new_poly_fit[i])
				
				for i in range(0,len(mshift)):
					residuals[i] = p(r_sorted_stat[i]) - mshift[i]
									
				for o in range (0, len(obs_list)):
					tmp_list = []
					for h in range(0,len(mshift)):
						if sorted_stats[h][23] == obs_list[o]:
							tmp_list.append(residuals[h])
							resid_per_obs.append(residuals[h])
							#tmp_list.append(mshift[h])
					mean_resid_per_observer[o] = statistics.mean(tmp_list)
				for h in range(0, len(mshift)):
					for o in range (0, len(obs_list)):
						if sorted_stats[h][23] == obs_list[o]:
						#	if first_pass !=1:
						#		file_writer.writerow([str(10**(r_sorted_stat[h])), str(sorted_stats[h][23]), str(mags_sorted_stat[h]),str(mshift[h]), str(p(r_sorted_stat[h])), str(residuals[h]), str(stdev_resid_per_observer[o])])	
						#	if first_pass ==1:
						#		file_writer.writerow([str(10**(r_sorted_stat[h])), str(sorted_stats[h][23]), str(mags_sorted_stat[h]),str(mshift[h]), str(p(r_sorted_stat[h])), str(residuals[h]), str(stdev_resid_per_observer[o])])	
							mshift[h] = mshift[h] + mean_resid_per_observer[o]
							break
							
				print(preorpost, k, new_poly_fit)
				
				#for h in range(0, len(mags_sorted_stat)):
				#	file_writer.writerow([str(r_sorted_stat[h]), str(sorted_stats[h][23]), str(mags_sorted_stat[h]),str(mshift[h]), str(p(r_sorted_stat[h]), str(mshift[h]-p(r_sorted_stat[h]), str(1)])	

						
				#if each of the coefficients are within 0.0001 then we say the polynomial has converged and are done calculating mshift
				if (abs(converge_test[0]) < tolerance) and (abs(converge_test[1])< tolerance) and (abs(converge_test[2])< tolerance) and (abs(converge_test[3])< tolerance) and (abs(converge_test[4])< tolerance) and (abs(converge_test[5])< tolerance):
					print(preorpost,': The polynomail fit converged to within tolerance of ', tolerance, ' after ', k, ' iterations')
					print('The final poly_fit is ', new_poly_fit)
					break
					
					
	if len(mshift) == 0 :
		for i in range (0,30):
			sorted_stats.append(0)
	return mshift, obs_list, sorted_stats, r_sorted_stat, new_poly_fit, original_poly_fit, stdev_resid_per_observer, mean_resid_per_observer, count_per_observer, residuals,  sorted_stats[28], mags_sorted_stat, resid_per_obs

#Assigns headers for output stats files
def add_headers_stats(inputlist, inputpreorpost,other):
	inputlist.insert(0, ['col 1-3 : short period comet designation'])
	inputlist[0].append('col 4-9 : Standard comet designation')
	inputlist[0].append('col 10 : multiple nuclei present?')
	inputlist[0].append('col 12-15 : year observed')
	inputlist[0].append('col 17-18 : month observed')
	inputlist[0].append('col 20-24')
	inputlist[0].append('col 26 : special note / extinction note')
	inputlist[0].append('col 27 : Magnitude collection method')
	inputlist[0].append('col 28-32 : visual magnitude estimate')
	inputlist[0].append('col 33 : poor conditions?')
	inputlist[0].append('col 34 - 35 : reference catalog')
	inputlist[0].append('col 36-40 : instrument aperture in centimeters')
	inputlist[0].append('col 41 : instrument type')
	inputlist[0].append('col 42 - 43 : focal ratio')
	inputlist[0].append('col 44-47 : magnification used')
	inputlist[0].append('col 49 : error estimate for coma diameter')
	inputlist[0].append('col 50 - 54 : coma diameter in arcminutes')
	inputlist[0].append('col 55 : special note on central condensation of comet')
	inputlist[0].append('col : 56 -57 : degree of condensation (note / means estimate)')
	inputlist[0].append('col 59 - 64 : error of tail approximation and tail approximation')
	inputlist[0].append('col 65 - 67 : direction tail is pointed')
	inputlist[0].append('col 69-74 : ICQ reference publication')
	inputlist[0].append('col 75 : second special note / extinction note ')
	inputlist[0].append('col 76-80 : observer name')
	inputlist[0].append('Date YYYY-MM-DDTHH:MM:SS')
	inputlist[0].append(inputpreorpost)
	inputlist[0].append('Delta (au)')
	inputlist[0].append('Phase Angle')
	inputlist[0].append(other)
	

def main():
	global metalist
	global list_of_reasons_removed
	global removed_metalist
	global reasonForDelete
	global to_report_r
	global heliocentric_corrected_magnitudes
	global phase_corrected_magnitudes
	global dates_pds_format
	global to_report_delta
	global to_report_phase
	global dates_pds_format
	global last_mag_calculated_sans_condemned_pre
	global last_mag_calculated_sans_condemned_post
	last_mag_calculated_sans_condemned_pre = []
	last_mag_calculated_sans_condemned_post =[]
	dates_pds_format = []
	#Instantiates each element of metalist (i.e., each column in the input data).
	metalist = []
	shortperapparition = []		#metalist[0]
	designation = []			#metalist[1]
	splitnuc = []				#metalist[2]
	yearobs = []				#metalist[3]
	monthobs = []				#metalist[4]
	dayobs = []					#metalist[5]
	speicalnotes = []			#metalist[6]
	magmethod = []				#metalist[7]
	mag = []					#metalist[8]
	poorconditions = []			#metalist[9]
	referencecat = []			#metalist[10]
	instaperture = []			#metalist[11]
	insttype = []				#metalist[12]
	focalratio = []				#metalist[13]
	magnification = []			#metalist[14]
	comadiamestimate = []		#metalist[15]
	comadiameter = []			#metalist[16]
	centralcondensation = []	#metalist[17]
	degreeofcondensation = []	#metalist[18]
	taillength = []				#metalist[19]
	positionangleoftail = []	#metalist[20]
	ICQPublication = []			#metalist[21]
	specialnotestwo = []		#metalist[22]
	obs = []					#metalist[23]
	removed = []
	removed_reason = []
	reasonForDelete = 0

	list_of_reasons_removed = ["Two entries on the same date by same observer", "No magnitude reported", "Used reverse binocular observing method", "Poor Weather Reported", "Used a tier 3 or 4 Source Catalog", "Used a telescope under 5.5 magnitude", "Used binoculars under 3.3 magnitude", "Did not use a magnitude method reported by Green (i.e. column 27 not being S, B, M, I, or E), prioritizing S then M", "Bad Extinction Correction used", "Observer used SC Catalog for object dimmer than 8.1"]

	#Reads in the 80 column format from ICQ or COBS data
	for line in open(input_file, encoding='utf8').readlines():
		shortperapparition.append(line[0:3].strip(' '))
		designation.append(line[3:9].strip(' '))
		splitnuc.append(line[9].strip(' '))
		yearobs.append(line[11:15].strip(' '))
		monthobs.append(line[16:18].strip(' '))
		dayobs.append(line[19:24].strip(' '))
		speicalnotes.append(line[25].strip(' '))
		magmethod.append(line[26].strip(' '))
		mag.append(line[28:32].strip(' '))
		poorconditions.append(line[32].strip(' '))
		referencecat.append(line[33:35].strip(' '))
		instaperture.append(line[35:40].strip(' '))
		insttype.append(line[40].strip(' '))
		focalratio.append(line[41:43].strip(' '))
		magnification.append(line[43:47].strip(' '))
		comadiamestimate.append(line[48].strip(' '))
		comadiameter.append(line[49:54].strip(' '))
		centralcondensation.append(line[54].strip(' '))
		degreeofcondensation.append(line[55:57].strip(' '))
		taillength.append(line[58:63].strip(' '))
		positionangleoftail.append(line[64:67].strip(' '))
		ICQPublication.append(line[68:74].strip(' '))
		specialnotestwo.append(line[74].strip(' '))
		obs.append(line[75:80].strip(' '))

	#Places each of the lists into one list for organization
	metalist = [shortperapparition,designation,splitnuc,yearobs,monthobs,dayobs,speicalnotes,
		   magmethod,mag,poorconditions,referencecat,instaperture,insttype,focalratio,magnification,
		   comadiamestimate,comadiameter,centralcondensation,degreeofcondensation,
		   taillength,positionangleoftail,ICQPublication,specialnotestwo,obs]

	removed_metalist = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]

	#Total number of initial datapoints
	print('initial number of points ', len(metalist[0]))
		  
	#Removes data from metalist based on specific criteria
	numdeltefornomagreported =0
	reasonForDelete = 1
	nomagreported(numdeltefornomagreported)

	numdeletedforreversebinocmethod =0
	reasonForDelete = 2
	reversebinocular(numdeletedforreversebinocmethod)

	numdeletedforpoorweather = 0
	reasonForDelete = 3
	poorweather(numdeletedforpoorweather)
	
	numdeletedforbadextinctioncorrection = 0
	reasonForDelete = 8
	badExtinctionCorrection(numdeletedforbadextinctioncorrection)

	numberremovedfortelescopeunder5_4 = 0
	numberremovedforbinocularunder1_4 = 0
	checkTelescopesandBinocMethods(numberremovedfortelescopeunder5_4, numberremovedforbinocularunder1_4)

	numberremovedfornotspecifiedmagmethod = 0
	if CCD_Bool == 1:
		reasonForDelete = 7
		removeForMagMethod(numberremovedfornotspecifiedmagmethod)

	numberremovedforSCcatalog = 0
	if CCD_Bool == 1:
		reasonForDelete = 9
		checkSCcatalog(numberremovedforSCcatalog)

	numdelforduplicatedate = 0
	reasonForDelete = 0
	deleteForDuplicatedDates(numdelforduplicatedate)

	#How many points are left in our data after sorting out 'rejected' points
	print("final remaining points " + str(len(metalist[2])))		

	#If you are not doing any further corrections to data then output "kept" points as is
	if "--heliocentric" not in sys.argv and '--phase' not in sys.argv:
		add_headers(metalist)
		file_writer = csv.writer(open(ouput_file_kept_points, 'w'), delimiter =',')
		for k in range (0,len(metalist[2])):
			file_writer.writerow([metalist[0][k],metalist[1][k],metalist[2][k],metalist[3][k],metalist[4][k],metalist[5][k],metalist[6][k],metalist[7][k],metalist[8][k],metalist[9][k],metalist[10][k],metalist[11][k],metalist[12][k],metalist[13][k],metalist[14][k],metalist[15][k],metalist[16][k],metalist[17][k],metalist[18][k],metalist[19][k],metalist[20][k],metalist[21][k],metalist[22][k],metalist[23][k]])

	#Outputs removed data points in separate csv along with reason it was deleted.
	add_headers(removed_metalist)
	removed_metalist[24].insert(0, 'Point removed')
	removed_metalist[25].insert(0, 'Reason Point was Removed')
		
	file_writer = csv.writer(open(output_file_rejected_points, 'w'), delimiter =',')
	for k in range (0,len(removed_metalist[2])):
		file_writer.writerow([removed_metalist[0][k],removed_metalist[1][k],removed_metalist[2][k],removed_metalist[3][k],removed_metalist[4][k],removed_metalist[5][k],removed_metalist[6][k],removed_metalist[7][k],removed_metalist[8][k],removed_metalist[9][k],removed_metalist[10][k],removed_metalist[11][k],removed_metalist[12][k],removed_metalist[13][k],removed_metalist[14][k],removed_metalist[15][k],removed_metalist[16][k],removed_metalist[17][k],removed_metalist[18][k],removed_metalist[19][k],removed_metalist[20][k],removed_metalist[21][k],removed_metalist[22][k],removed_metalist[23][k], removed_metalist[24][k], removed_metalist[25][k]])	
		
	#Optional command line argument --heliocentric to perform just heliocentric corrections to 'kept' data
	if "--heliocentric" in sys.argv and '--phase' not in sys.argv:
		print('Performing Heliocentric Corrections to the Data')
		queryJPL()
		heliocentric_corrected_magnitudes = []
		to_report_delta = []
		to_report_r = []
		to_report_phase = []
		#queryJPL will report an r, delta, and phase angle at every 30 minute increment in the ephemerides
		#It will also take each date/time of an observation in the dataset and convert it to YYYY:MM:DD HH:MM:SS format.
		#The next few lines will compare the date/time of each point in the observation and find the nearest 30 minute increment in the ephemerides
		#It will then take the delta at that nearest increment and use it to apply a heliocentric correction to the magnitude,
		#and repeat this for each point in the data.
		for i in range (0, len(date_compare_to_JPL)):
			for j in range(0, len(OBJDates)):
				if (date_compare_to_JPL[i][0:4] == OBJDates[j][0:4]) and (date_compare_to_JPL[i][5:7] == OBJDates[j][5:7]) and (date_compare_to_JPL[i][8:10] == OBJDates[j][8:10]) and (date_compare_to_JPL[i][11:13] == OBJDates[j][11:13]) and ((JPL_Time_Increment*round(float(date_compare_to_JPL[i][14:16])/JPL_Time_Increment))%60 == float(OBJDates[j][14:16])):
					heliocentric_corrected_magnitudes.append(str(float(metalist[8][i]) - 5 * float(math.log10(OBJDelta[j]))))
					to_report_r.append(OBJr[j])
					to_report_delta.append(OBJDelta[j])
					to_report_phase.append(OBJPhase[j])
					continue
		for k in range(0,len(date_compare_to_JPL)):
			dates_pds_format.append(date_compare_to_JPL[k].replace(" ","T"))

		#writes out final heliocentric corrected data.
		add_headers(metalist)
		to_report_r.insert(0, 'Heliocentric Distance (au)')
		dates_pds_format.insert(0, 'Dates YYYY:MM:DDTHH:MM:SS')
		to_report_delta.insert(0, 'Delta (au)')
		to_report_phase.insert(0, 'Phase angle')
		heliocentric_corrected_magnitudes.insert(0, 'magnitdues with only heliocentric correction (mhelio)')
		file_writer = csv.writer(open(ouput_file_kept_points, 'w'), delimiter =',')
		for k in range (0,len(metalist[2])):
			file_writer.writerow([metalist[0][k],metalist[1][k],metalist[2][k],metalist[3][k],metalist[4][k],metalist[5][k],metalist[6][k],metalist[7][k],metalist[8][k],metalist[9][k],metalist[10][k],metalist[11][k],metalist[12][k],metalist[13][k],metalist[14][k],metalist[15][k],metalist[16][k],metalist[17][k],metalist[18][k],metalist[19][k],metalist[20][k],metalist[21][k],metalist[22][k],metalist[23][k], to_report_r[k], heliocentric_corrected_magnitudes[k]])

	#Optional command line argument --phase to perform just phase corrections to 'kept' data
	if '--phase' in sys.argv and '--heliocentric' not in sys.argv:
		print('Performing Phase Angle Corrections to the Data')
		queryJPL()
		phase_corrected_magnitudes = []
		to_report_delta = []
		to_report_r = []
		to_report_phase = []
		#queryJPL will report an r, delta, and phase angle at every 30 minute increment in the ephemerides
		#It will also take each date/time of an observation in the dataset and convert it to YYYY:MM:DD HH:MM:SS format.
		#The next few lines will compare the date/time of each point in the observation and find the nearest 30 minute increment in the ephemerides
		#It will also read in Schleicher's composite phse function which has inputs at every whole number degree.
		#It will then take the phase at that nearest increment, round it to the nearest whole number, compare that to Shcleicher's data to get the
		#composite phase function normalized to 0 degrees, and use that to calculate the phase corrected magnitudes.
		#This is repeated for each point in the data.
		with open('Schleicher_Composite_Phase_Function.txt') as f:
			lines1 = f.readlines()
			phase_angles = [line.split()[0] for line in lines1]
			deg_0_normalized = [float(line.split()[1]) for line in lines1]
		for i in range(0, len(date_compare_to_JPL)):
			for j in range(0, len(OBJPhase)):
				if (date_compare_to_JPL[i][0:4] == OBJDates[j][0:4]) and (date_compare_to_JPL[i][5:7] == OBJDates[j][5:7]) and (date_compare_to_JPL[i][8:10] == OBJDates[j][8:10]) and (date_compare_to_JPL[i][11:13] == OBJDates[j][11:13]) and ((JPL_Time_Increment*round(float(date_compare_to_JPL[i][14:16])/JPL_Time_Increment))%60 == float(OBJDates[j][14:16])):
					to_report_r.append(OBJr[j])
					to_report_delta.append(OBJDelta[j])
					to_report_phase.append(OBJPhase[j])
					for l in range (0, len(phase_angles)):
						if (round(float(OBJPhase[j])) == float(phase_angles[l])):
							phase_corrected_magnitudes.append(str(float(metalist[8][i]) + 2.5 * float(math.log10(deg_0_normalized[l]))))
							continue
				
		for k in range(0,len(date_compare_to_JPL)):
			dates_pds_format.append(date_compare_to_JPL[k].replace(" ","T"))
				
		#writes out final phase corrected data
		add_headers(metalist)
		to_report_r.insert(0, 'Heliocentric Distance (au)')
		dates_pds_format.insert(0, 'Dates YYYY:MM:DDTHH:MM:SS')
		to_report_delta.insert(0, 'Delta (au)')
		to_report_phase.insert(0, 'Phase angle')
		phase_corrected_magnitudes.insert(0, 'magnitudes with only phase correction (mph*)')
		file_writer = csv.writer(open(ouput_file_kept_points, 'w'), delimiter =',')
		for k in range (0,len(metalist[2])):
			file_writer.writerow([metalist[0][k],metalist[1][k],metalist[2][k],metalist[3][k],metalist[4][k],metalist[5][k],metalist[6][k],metalist[7][k],metalist[8][k],metalist[9][k],metalist[10][k],metalist[11][k],metalist[12][k],metalist[13][k],metalist[14][k],metalist[15][k],metalist[16][k],metalist[17][k],metalist[18][k],metalist[19][k],metalist[20][k],metalist[21][k],metalist[22][k],metalist[23][k], to_report_r[k], phase_corrected_magnitudes[k]])
	
	#Performs a heliocentric correction to the raw data and then a phase correction to the heliocentric corrected data
	#See the above two blocks to understand how the heliocentric and phase corrections work
	if '--heliocentric' in sys.argv and '--phase' in sys.argv:
		print('Performing heliocentric and Phase Angle Corrections to the Data')
		queryJPL()
		heliocentric_corrected_magnitudes = []
		phase_corrected_magnitudes = []
		to_report_r = []
		to_report_delta = []
		to_report_phase = []
		with open('Schleicher_Composite_Phase_Function.txt') as f:
			lines1 = f.readlines()
			phase_angles = [line.split()[0] for line in lines1]
			deg_0_normalized = [float(line.split()[1]) for line in lines1]
		for i in range (0, len(date_compare_to_JPL)):
			for j in range(0, len(OBJDates)):
				if (date_compare_to_JPL[i][0:4] == OBJDates[j][0:4]) and (date_compare_to_JPL[i][5:7] == OBJDates[j][5:7]) and (date_compare_to_JPL[i][8:10] == OBJDates[j][8:10]) and (date_compare_to_JPL[i][11:13] == OBJDates[j][11:13]) and ((JPL_Time_Increment*round(float(date_compare_to_JPL[i][14:16])/JPL_Time_Increment))%60 == float(OBJDates[j][14:16])):
					to_report_r.append(OBJr[j])
					to_report_delta.append(OBJDelta[j])
					to_report_phase.append(OBJPhase[j])
					heliocentric_corrected_magnitudes.append(str(float(metalist[8][i]) - 5. * float(math.log10(OBJDelta[j]))))
					for l in range (0, len(phase_angles)):
						if (round(float(OBJPhase[j])) == float(phase_angles[l])):
							phase_corrected_magnitudes.append(str(float(heliocentric_corrected_magnitudes[i]) + 2.5 * float(math.log10(deg_0_normalized[l]))))

		for k in range(0,len(date_compare_to_JPL)):
			dates_pds_format.append(date_compare_to_JPL[k].replace(" ","T"))
		#writes out the heliocentric and phase corrected magnitudes
		add_headers(metalist)
		to_report_r.insert(0, 'Heliocentric Distance (au)')
		heliocentric_corrected_magnitudes.insert(0, 'heliocentric corrected magnitudes (mhelio)')
		phase_corrected_magnitudes.insert(0, 'magnitudes with heliocentric and phase corrections applied (mph)')
		dates_pds_format.insert(0, 'Dates YYYY:MM:DDTHH:MM:SS')
		to_report_delta.insert(0, 'Delta (au)')
		to_report_phase.insert(0, 'Phase angle')
		file_writer = csv.writer(open(ouput_file_kept_points, 'w'), delimiter =',')
		for k in range (0,len(metalist[2])):
			file_writer.writerow([metalist[0][k],metalist[1][k],metalist[2][k],metalist[3][k],metalist[4][k],metalist[5][k],metalist[6][k],metalist[7][k],metalist[8][k],metalist[9][k],metalist[10][k],metalist[11][k],metalist[12][k],metalist[13][k],metalist[14][k],metalist[15][k],metalist[16][k],metalist[17][k],metalist[18][k],metalist[19][k],metalist[20][k],metalist[21][k],metalist[22][k],metalist[23][k], to_report_r[k], heliocentric_corrected_magnitudes[k], phase_corrected_magnitudes[k], dates_pds_format[k], to_report_delta[k], to_report_phase[k]])
		
	#Performs all statistical corrections outlined in 'Statistics_method_appendix.txt' in GitHub repository
	#All TRY-EXCEPT blocks are case scenarios depending on whether the user calculated mph, mehlio, or both.
	if '--stats' in sys.argv:
		last_mag_calculated = []
		pre_condemned_obs = []
		post_condemned_obs = []
		other_mag = []
		magsfound = 0 	#checks what was last magnitude calculated (either mph or mhelio depending on if user calculated one, neither, or both)
		first_pass = 1
		other = ''
		try:
			if magsfound == 0:
				last_mag_calculated = phase_corrected_magnitudes
				print('Will perform statistical corrections on Phase Corrected Magnitudes')
				magsfound = 'mph'
				try:
					other = 'mhelio'
					other_mag = heliocentric_corrected_magnitudes
					try:
						float(other_mag[0])
					except:
						del other_mag[0]
				except:
					pass
		except:
			print('Searching for heliocentric corrected magnitudes...')
			
		try:
			if magsfound == 0:
				last_mag_calculated = heliocentric_corrected_magnitudes
				print('Will perform statistical corrections on Heliocentric Corrected Magnitudes')
				magsfound = 'mhelio'
				try:
					other = 'mph'
					other_mag = phase_corrected_magnitudes
					try:
						float(other_mag[0])
					except:
						del other_mag[0]
				except:
					pass
		except:
			print('Searching for raw magnitudes...')
			
		if magsfound == 0:
			print('Please run either --heliocentric or --phase or both to perform statistical corrections')
			
		#deleting headers now so we dont have to worry about them in the stats function
		del last_mag_calculated[0]
		del to_report_r[0]
		del to_report_delta[0]
		del to_report_phase[0]
		del dates_pds_format[0]
		deletearow(0)
		
		print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')	
		
		#Initializes and defines pre-perihelion statistical outputs
		pre_mshift = []  
		pre_r = []		
		pre_meta = []   
		pre_final_polyfit = [] 
		pre_final_stdevs = []	
		pre_last_mag_correction = []	
		pre_original_polyfit = []	
		pre_obs_list = []   
		pre_final_mean_resid = []	
		pre_count_per_obs = []
		pre_other_mag = []
		pre_last_mag_calculated = []
		pre_mshift, pre_obs_list, pre_meta, pre_r,  pre_final_polyfit, pre_original_polyfit, pre_final_stdevs, pre_final_mean_resid, pre_count_per_obs, pre_last_mag_correction, pre_other_mag, pre_last_mag_calculated, pre_resid_per_obs= stats_shifts('pre', metalist, last_mag_calculated, dates_pds_format, OBJDelta, OBJPhase, to_report_r,pre_condemned_obs, other_mag, first_pass)
				
		print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
		
		#Initializes and defines post-perihelion statistical outputs
		post_mshift = []
		post_r = []
		post_meta = []
		post_final_polyfit = []
		post_final_stdevs = []
		post_last_mag_correction = []
		post_original_poly_fit = []
		post_final_mean_resid = []
		post_obs_list = []
		post_count_per_obs = []
		post_other_mag = []
		post_last_mag_calculated =[]
		post_mshift, post_obs_list, post_meta, post_r,  post_final_polyfit, post_original_polyfit, post_final_stdevs, post_final_mean_resid, post_count_per_obs, post_last_mag_correction, post_other_mag, post_last_mag_calculated, post_resid_per_obs= stats_shifts('post',  metalist, last_mag_calculated, dates_pds_format, OBJDelta, OBJPhase, to_report_r,post_condemned_obs, other_mag,first_pass)	
		
		#sets this to 0 so that if we need to run stats_shifts again then we wont accidentally log the r-values again
		first_pass = 0
		
		#Block begins to perform stationary t-test and two tail probability test. This block is for pre-perihelion t-test
		#terminate iterations = 0 means at least one observer failed the t-test, so we must reconverge a polynomial fit (i.e., run stats_shifts) with their data removed
		#number_t_pre keeps track of how many times we have at least one observer fail a t-test on a given iteration
		terminate_iterations = 0
		number_t_pre = 0
		while terminate_iterations == 0:
			drop_observers = 0
			thismagsfound = magsfound
			p_func = np.poly1d(pre_final_polyfit)
			for o in range (0, len(pre_obs_list)):
				N1 = []
				N2 = []
				mean_1 = 0
				mean_2 = 0
				sum_1 = 0
				sum_2 = 0
				for i in range (0, len(pre_r)):
								
					if pre_meta[i][23] == pre_obs_list[o]:
						if len(N1) < math.floor(pre_count_per_obs[o] /2):
							N1.append(p_func(pre_r[i]) - pre_mshift[i])
						else:
							N2.append(p_func(pre_r[i]) - pre_mshift[i])
				
				#returns t-statistic and corresponding p-statistics
				t2, p2 = stats.ttest_ind(N1, N2, equal_var = False)
				print(pre_obs_list[o], ' t = ',t2,' p = ', p2)
				#Adds observers who failed p-test to 'condemned list' to be avoided on future polynomial fit convergeances.
				if p2 < 0.05:
					drop_observers = 1
					pre_condemned_obs.append(pre_obs_list[o])
			#If we did deleted one observer, prepare their data to be in write format for stats_shifts and then rerun stats_shifts
			#i.e., one observer failed the stationary test so we reconverge the polynomial fit / mshift values without the bias from their data present
			if drop_observers ==1:
				tmp = []
				last_mag_calculated_sans_condemned_pre = []
				for z in range(0,len(pre_last_mag_calculated)):
					last_mag_calculated_sans_condemned_pre.append(pre_last_mag_calculated[z])
				print(len(last_mag_calculated_sans_condemned_pre))
				for h in range (len(pre_last_mag_calculated)-1, -1, -1):
					if pre_meta[h][23] in pre_condemned_obs:
						del pre_last_mag_calculated[h]
				for i in range(0,len(pre_meta[1])):
					for j in range(0,len(pre_meta)):
						if j == 0:
							tmp.append([])
						tmp[i].append(pre_meta[j][i])
				tmp_dates = tmp[24]
				tmp_deltas = tmp[26]
				tmp_phase = tmp[27]
				tmp_other_mags = tmp[28]
				del tmp[28]
				del tmp[27]
				del tmp[26]
				del tmp[25]
				del tmp[24]
				number_t_pre = number_t_pre +1
				print('###########################################################################################')
				print('PRE: THE FOLLOWING OBSERVERS WERE REJECTED BY T-TEST: ')
				print(pre_condemned_obs)
				print('DROPPING THESE OBSERVERS AND REPEATING THE ITERATIVE PROCESS')
				print('###########################################################################################')
				pre_mshift, pre_obs_list, pre_meta, pre_r,  pre_final_polyfit, pre_original_polyfit, pre_final_stdevs, pre_final_mean_resid, pre_count_per_obs, pre_last_mag_correction,pre_other_mag, tmp, pre_resid_per_obs = stats_shifts('pre', tmp, last_mag_calculated_sans_condemned_pre, tmp_dates, tmp_deltas, tmp_phase, pre_r, pre_condemned_obs, tmp_other_mags, first_pass)
			else:
				terminate_iterations = 1
		
		#Repeats the above t- and p-tests for Post-perihelion data
		terminate_iterations = 0
		number_t_post = 0
		print('###########################################################################################')
		while terminate_iterations == 0:
			drop_observers = 0
			thismagsfound = magsfound
			p_func = np.poly1d(post_final_polyfit)
			for o in range (0, len(post_obs_list)):
				N1 = []
				N2 = []
				mean_1 = 0
				mean_2 = 0
				sum_1 = 0
				sum_2 = 0
				for i in range (0, len(post_r)):
					if post_meta[i][23] == post_obs_list[o]:
						if len(N1) < math.floor(post_count_per_obs[o] /2):
							N1.append(p_func(post_r[i]) - post_mshift[i])
						else:
							N2.append(p_func(post_r[i]) - post_mshift[i])
			
				t2, p2 = stats.ttest_ind(N1, N2, equal_var = False)
				print(post_obs_list[o], ' t = ', t2,' p = ', p2) #t1, 2*p1)
				if p2 < 0.05:
					drop_observers = 1
					post_condemned_obs.append(post_obs_list[o])
			if drop_observers ==1:
				last_mag_calculated_sans_condemned_post = []
				for z in range(0,len(post_last_mag_calculated)):
					last_mag_calculated_sans_condemned_post.append(post_last_mag_calculated[z])
				print(len(last_mag_calculated_sans_condemned_post))
				for h in range (len(post_last_mag_calculated)-1, -1, -1):
					if post_meta[h][23] in post_condemned_obs:
						del post_last_mag_calculated[h]
				tmp = []
				for i in range(0,len(post_meta[1])):
					for j in range(0,len(post_meta)):
						if j == 0:
							tmp.append([])
						tmp[i].append(post_meta[j][i])
				tmp_dates = tmp[24]
				tmp_deltas = tmp[26]
				tmp_phase = tmp[27]
				tmp_other_mags = tmp[28]
				del tmp[28]
				del tmp[27]
				del tmp[26]
				del tmp[25]
				del tmp[24]
				last_mag_calculated_sans_condemned = []
				print('#################################################################################')
				print('POST: THE FOLLOWING OBSERVERS WERE REJECTED BY T-TEST: ')
				print(post_condemned_obs)
				print('DROPPING THESE OBSERVERS AND REPEATING THE ITERATIVE PROCESS FOR POSTPERIHELION')
				print('#################################################################################')
				post_mshift, post_obs_list, post_meta, post_r,  post_final_polyfit, post_original_polyfit, post_final_stdevs, post_final_mean_resid, post_count_per_obs, post_last_mag_correction,post_other_mag, tmp, post_resid_per_obs = stats_shifts('post', tmp, last_mag_calculated_sans_condemned_post, tmp_dates, tmp_deltas, tmp_phase, post_r, post_condemned_obs, tmp_other_mags, first_pass)
			else:
				terminate_iterations = 1
				
				
		print(number_t_pre, 'pre t tests', number_t_post, 'post t tests')
		print('PRE: observers who failed t-test and were removed: ', pre_condemned_obs)
		print('POST: observers who failed t-test and were removed: ', post_condemned_obs)				
	
		#Adds headers and writes out pre-perihelion data to file 'pre-stats.csv'
		if (len(pre_meta) != 0) and (pre_meta[0] != 0):
			magsfound = 'mshift (no dropped observers)'
			add_headers_stats(pre_meta, magsfound,other)
			if type(pre_r) != list:
				pre_r = pre_r.tolist()
			pre_r.insert(0, 'r (au)')
			pre_last_mag_calculated.insert(0, thismagsfound)
			pre_last_mag_correction.insert(0, 'residual of mshift from polyfit')
			pre_meta[0][28] = other + ' (BLANK if you did not ask to calculate this value)'
			pre_mshift.insert(0, 'mshift with dropped observers')
			file_writer = csv.writer(open('pre-stats.csv', 'w'), delimiter =',')
			for m in range (1, len(pre_r)):
				pre_r[m] = str((-1.0)*10**(float(pre_r[m])))
			for k in range (0,len(pre_meta)):
				file_writer.writerow([pre_meta[k][0],pre_meta[k][1],pre_meta[k][2],pre_meta[k][3],pre_meta[k][4],pre_meta[k][5],pre_meta[k][6],pre_meta[k][7],pre_meta[k][8],pre_meta[k][9],pre_meta[k][10],pre_meta[k][11],pre_meta[k][12],pre_meta[k][13],pre_meta[k][14],pre_meta[k][15],pre_meta[k][16],pre_meta[k][17],pre_meta[k][18],pre_meta[k][19],pre_meta[k][20],pre_meta[k][21],pre_meta[k][22],pre_meta[k][23], pre_meta[k][24],pre_r[k],pre_meta[k][28], pre_last_mag_calculated[k], pre_meta[k][25], pre_mshift[k], pre_meta[k][26], pre_meta[k][27], pre_last_mag_correction[k]])	
		else:
			print('No preperihelion data to perform statistics on')
			print('##########################################################################################')

		#Adds headers and writes out post-perihelion data to file 'post-stats.csv'
		if (len(post_meta) != 0) and (post_meta[0] != 0):
			magsfound = 'mshift (no dropped observers)'
			add_headers_stats(post_meta, magsfound,other)
			if type(post_r) != list:
				post_r = post_r.tolist()
			post_r.insert(0, 'r (au)')
			post_last_mag_calculated.insert(0, thismagsfound)
			post_last_mag_correction.insert(0, 'residual of mshift from polyfit')
			post_meta[0][28] = other + ' (BLANK if you did not ask to calculate this value)'
			post_mshift.insert(0, 'mshift with dropped observers')
			file_writer = csv.writer(open('post-stats.csv', 'w'), delimiter =',')
			for m in range (1, len(post_r)):
				post_r[m] = str(10**(float(post_r[m])))
			for k in range (0,len(post_meta)):
				file_writer.writerow([post_meta[k][0],post_meta[k][1],post_meta[k][2],post_meta[k][3],post_meta[k][4],post_meta[k][5],post_meta[k][6],post_meta[k][7],post_meta[k][8],post_meta[k][9],post_meta[k][10],post_meta[k][11],post_meta[k][12],post_meta[k][13],post_meta[k][14],post_meta[k][15],post_meta[k][16],post_meta[k][17],post_meta[k][18],post_meta[k][19],post_meta[k][20],post_meta[k][21],post_meta[k][22],post_meta[k][23], post_meta[k][24],post_r[k],post_meta[k][28],post_last_mag_calculated[k], post_meta[k][25], post_mshift[k], post_meta[k][26], post_meta[k][27], post_last_mag_correction[k]])	
		else:
			print('No postperihelion data to perform statistics on')
			print('###########################################################################################')

	#Plots all available data. That is, any combination of mraw, mhelio, mph, and mshifts depending on 
	#which combination of those the user has calculated (i.e., running --heliocentric --shifts will only plot mraw, mhelio, and mshift).
	#All if statements and TRY - EXCEPT blocks are checks
	#for each possible scenario.
	if '--plot' in sys.argv:

		if ('--heliocentric' not in sys.argv) and ('--phase' not in sys.argv):
			print('Please perform --heliocentric, --phase, or both before attempting to plot')
			sys.exit()
			
		mags_to_plot_meta = []
		max_x_values = []
		min_x_values = []
		max_y_values = []
		min_y_values = []
		to_check_min_max_x = []
		to_check_min_max_y = []
		titles = []
		axis = []
		count = 0
				
		if '--stats' not in sys.argv:
			deletearow(0)
			del to_report_r[0]
		for j in range(0, len(to_report_r)):
			newdate4 = time.strptime(perihelion, "%Y/%m/%d")
			datetocheck = str(metalist[3][j]) + "/" + str(metalist[4][j]) + "/" + str(math.floor(float(metalist[5][j])))
			newdate5 = time.strptime(datetocheck, "%Y/%m/%d")
			if (newdate5 <= newdate4):
				to_report_r[j] = float(-1. * to_report_r[j])
			
		try:
			tmpmeta, tmp_mags, tmp_r = sortbyr(metalist,to_report_r,metalist[8],1)
			for i in range (0, len(tmp_mags)):
				tmp_mags[i] = float(tmp_mags[i])
				tmp_r[i] = float(tmp_r[i])
			mags_to_plot_meta.append(tmp_mags)
			mags_to_plot_meta.append(tmp_r)
			count = count +1
			print('mraw found, adding to plot')
			titles.append('Reported Visual Magnitude Vs. Heliocentric Distance')
			axis.append('mraw')
		except:
			print('Please perform --heliocentric, --phase, or both before plotting')
			
		try:
			try:
				float(heliocentric_corrected_magnitudes[0])
			except:
				del heliocentric_corrected_magnitudes[0]
			tmpmeta, tmp_mags, tmp_r = sortbyr(metalist,to_report_r,heliocentric_corrected_magnitudes,1)
			for i in range (0, len(tmp_mags)):
				tmp_mags[i] = float(tmp_mags[i])
				tmp_r[i] = float(tmp_r[i])
			mags_to_plot_meta.append(tmp_mags)
			mags_to_plot_meta.append(tmp_r)
			print('mhelio found, adding to plot')
			titles.append('Heliocentric Corrected Magnitudes')
			axis.append('mhelio')
			count = count + 1
		except:
			print('mhelio not found, looking for other magnitudes to plot...')
			
		try:
			try:
				float(phase_corrected_magnitudes[0])
			except:
				del phase_corrected_magnitudes[0]
			tmpmeta, tmp_mags, tmp_r = sortbyr(metalist,to_report_r,phase_corrected_magnitudes,1)
			for i in range (0, len(tmp_mags)):
				tmp_mags[i] = float(tmp_mags[i])
				tmp_r[i] = float(tmp_r[i])
			mags_to_plot_meta.append(tmp_mags)
			mags_to_plot_meta.append(tmp_r)
			count = count + 1
			print('mph found, adding to plot')
			if '--heliocentric' in sys.argv:
				titles.append('Heliocentric and Phase Corrected Magnitudes')
			else:
				titles.append('Phase Corrected Magnitudes')
			axis.append('mph')
		except:
			print('mph not found, looking for other magnitudes to plot...')
			
		try:
			try:
				del pre_mshift[0]
			except:
				pass
			try:
				del post_mshift[0]
			except:
				pass
			try:
				del pre_r[0]
			except:
				pass
			try:
				del post_r[0]
			except:
				pass
			pre_and_post_shift_mags = pre_mshift + post_mshift
			pre_and_post_shift_r = pre_r + post_r
			mags_to_plot_meta.append(pre_and_post_shift_mags)
			mags_to_plot_meta.append(pre_and_post_shift_r)
			for i in range(0,len(pre_and_post_shift_mags)):
				pre_and_post_shift_mags[i] = float(pre_and_post_shift_mags[i])
				pre_and_post_shift_r[i] = float(pre_and_post_shift_r[i])
			count = count + 1
			print('Statistically corrected magnitudes found, adding to plot')
			axis.append('mshift')
			if ('--heliocentric' in sys.argv )and ('--phase' in sys.argv):
				titles.append('Heliocentric, Phase, and Statistical Corrected Magnitudes')
			elif '--heliocentric' in sys.argv:
				titles.append('Heliocentric and Statistically Corrected Magnitudes')
			elif '--phase' in sys.argv:
				titles.append('Phase and Statistically Corrected Magnitudes')
		except:
			print('No statistical corrections found, now plotting...')
			
		for i in range(0, len(mags_to_plot_meta)):
			if len(mags_to_plot_meta[i]) < len(mags_to_plot_meta[0]):
				continue
			if i%2 == 0:
				for j in range(0, len(mags_to_plot_meta[0])):
					to_check_min_max_y.append(float(mags_to_plot_meta[i][j]))
			if i%2 == 1:
				for j in range(0, len(mags_to_plot_meta[1])):
					to_check_min_max_x.append(float(mags_to_plot_meta[i][j]))
				
		true_max_x = math.ceil(max(to_check_min_max_x))
		true_min_x = math.floor(min(to_check_min_max_x))
		true_max_y = math.ceil(max(to_check_min_max_y))
		true_min_y = math.floor(min(to_check_min_max_y))
		
		#print(len(mags_to_plot_meta))
		for i in range(0,len(mags_to_plot_meta)):
			if i%2 == 1:
				continue
			
			fig = plt.figure(figsize=(1406./96., 897./96.))
			m1 = MultipleLocator(0.25)
			m2 = MultipleLocator(0.25)
			ax = fig.add_subplot(111)
			ax.yaxis.set_minor_locator(m2)
			ax.xaxis.set_minor_locator(m1)
			ax.set_xticks(np.arange(true_min_x -1, true_max_x+1,1))
			ax.set_yticks(np.arange(true_min_y -1, true_max_y+1,1))
			ax.tick_params(which='major', length=10, width=1, bottom=True, top=True, left=True, right=True)
			ax.tick_params(which='minor', length=5, width=1,labelbottom=True,bottom=True,top=True, left=True, right=True)
			ax.tick_params(which='both', direction='in', labelsize='15')
			ax.set_xlim(true_min_x-1, true_max_x+1)
			ax.set_ylim(true_min_y-1, true_max_y+1)
			ax.plot(mags_to_plot_meta[i+1], mags_to_plot_meta[i], '+', mew=1.1, ms=13, color='blue')
			ax.invert_yaxis()
			ax.set_title(titles[int(i/2)], fontsize='20')
			ax.set_xlabel('r (au)', fontsize='20')
			ax.set_ylabel(axis[int(i/2)], fontsize='20')
			dir_path = os.path.dirname(os.path.realpath(__file__))
			#title = dir_path+ "\" + target_nickname + "_graph_"+str(int(i/2))
			title = dir_path + '\_' + target_nickname + '_graph_'+str(int(i/2))
			canvas = FigureCanvas(fig)
			canvas.print_figure(title, dpi=96, bbox_inches='tight')
			ax.cla()
			
if __name__ == '__main__':
	main()
