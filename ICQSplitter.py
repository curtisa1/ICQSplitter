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

curtisa1 (at) mail.usf.edu, latest version: v1.0, 2018-19-19


*	v1.0: Sorts problematic entries from data, performs heliocentric and phase corrections.



"""

###############################
####### Input Arguments #######
###############################

input_file = 'input_data.txt'			#Name of your input file
small_body_designation = '902008;'			#Name of your small body ex) 'ceres' or 'eris'
JPL_Time_Increment = 30 					#How much to increment JPL queries in minutes up to 60.
ouput_file_kept_points = 'keepers.csv'		#Name of output file for points that meet all sorting criterion
output_file_rejected_points = 'removed.csv'	#Name of output file for points that were removed from the data

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
			elif ((metalist[7][k-1] == "B") or (metalist[7][k-1] == "I")):
				numdelforduplicatedate = numdelforduplicatedate +1
				addToremoved(k-1,reasonForDelete)
				deletearow(k-1)
				continue
			elif ((metalist[7][k] == "B") or (metalist[7][k] == "I")):
				numdelforduplicatedate = numdelforduplicatedate +1
				addToremoved(k-1,reasonForDelete)
				deletearow(k)
				continue
			else:
				numdelforduplicatedate = numdelforduplicatedate +1
				addToremoved(k-1,reasonForDelete)
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
	return date_compare_to_JPL

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
	
def main():
	global metalist
	global list_of_reasons_removed
	global removed_metalist
	global reasonForDelete
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
	for line in open(input_file).readlines():
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
	reasonForDelete = 7
	removeForMagMethod(numberremovedfornotspecifiedmagmethod)

	numberremovedforSCcatalog = 0
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
		to_report_r = []
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
					continue
		
		#writes out final heliocentric corrected data.
		add_headers(metalist)
		to_report_r.insert(0, 'Heliocentric Distance (au)')
		heliocentric_corrected_magnitudes.insert(0, 'magnitdues with only heliocentric correction (mhelio)')
		file_writer = csv.writer(open(ouput_file_kept_points, 'w'), delimiter =',')
		for k in range (0,len(metalist[2])):
			file_writer.writerow([metalist[0][k],metalist[1][k],metalist[2][k],metalist[3][k],metalist[4][k],metalist[5][k],metalist[6][k],metalist[7][k],metalist[8][k],metalist[9][k],metalist[10][k],metalist[11][k],metalist[12][k],metalist[13][k],metalist[14][k],metalist[15][k],metalist[16][k],metalist[17][k],metalist[18][k],metalist[19][k],metalist[20][k],metalist[21][k],metalist[22][k],metalist[23][k], to_report_r[k], heliocentric_corrected_magnitudes[k]])

	#Optional command line argument --phase to perform just phase corrections to 'kept' data
	if '--phase' in sys.argv and '--heliocentric' not in sys.argv:
		print('Performing Phase Angle Corrections to the Data')
		queryJPL()
		phase_corrected_magnitudes = []
		to_report_r = []
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
					for l in range (0, len(phase_angles)):
						if (round(float(OBJPhase[j])) == float(phase_angles[l])):
							phase_corrected_magnitudes.append(str(float(metalist[8][i]) + 2.5 * float(math.log10(deg_0_normalized[l]))))
							continue
							
		#writes out final phase corrected data
		add_headers(metalist)
		to_report_r.insert(0, 'Heliocentric Distance (au)')
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
		with open('Schleicher_Composite_Phase_Function.txt') as f:
			lines1 = f.readlines()
			phase_angles = [line.split()[0] for line in lines1]
			deg_0_normalized = [float(line.split()[1]) for line in lines1]
		for i in range (0, len(date_compare_to_JPL)):
			for j in range(0, len(OBJDates)):
				if (date_compare_to_JPL[i][0:4] == OBJDates[j][0:4]) and (date_compare_to_JPL[i][5:7] == OBJDates[j][5:7]) and (date_compare_to_JPL[i][8:10] == OBJDates[j][8:10]) and (date_compare_to_JPL[i][11:13] == OBJDates[j][11:13]) and ((JPL_Time_Increment*round(float(date_compare_to_JPL[i][14:16])/JPL_Time_Increment))%60 == float(OBJDates[j][14:16])):
					to_report_r.append(OBJr[j])
					heliocentric_corrected_magnitudes.append(str(float(metalist[8][i]) - 5 * float(math.log10(OBJDelta[j]))))
					for l in range (0, len(phase_angles)):
						if (round(float(OBJPhase[j])) == float(phase_angles[l])):
							phase_corrected_magnitudes.append(str(float(heliocentric_corrected_magnitudes[i]) + 2.5 * float(math.log10(deg_0_normalized[l]))))

		#writes out the heliocentric and phase corrected magnitudes
		add_headers(metalist)
		to_report_r.insert(0, 'Heliocentric Distance (au)')
		heliocentric_corrected_magnitudes.insert(0, 'heliocentric corrected magnitudes (mhelio)')
		phase_corrected_magnitudes.insert(0, 'magnitudes with heliocentric and phase corrections applied (mph)')
		file_writer = csv.writer(open(ouput_file_kept_points, 'w'), delimiter =',')
		for k in range (0,len(metalist[2])):
			file_writer.writerow([metalist[0][k],metalist[1][k],metalist[2][k],metalist[3][k],metalist[4][k],metalist[5][k],metalist[6][k],metalist[7][k],metalist[8][k],metalist[9][k],metalist[10][k],metalist[11][k],metalist[12][k],metalist[13][k],metalist[14][k],metalist[15][k],metalist[16][k],metalist[17][k],metalist[18][k],metalist[19][k],metalist[20][k],metalist[21][k],metalist[22][k],metalist[23][k], to_report_r[k], heliocentric_corrected_magnitudes[k], phase_corrected_magnitudes[k]])
		
if __name__ == '__main__':
	main()
