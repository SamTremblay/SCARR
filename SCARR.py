## SCARR v1.0
## Copyright 2018, Samuel Tremblay-Belzile

import numpy
import sys
import getopt
from operator import itemgetter

#Default input data file
#blastfile = 'blast.txt'

#Default output results file
summaryfile = 'AnalysisSummary.txt'
statfile = 'JunctionStats.txt'
singjuncfile = 'SingleJunctions.txt'
doubjuncfile = 'DoubleJunctions.txt'
unlabelled = 'UnlabelledReads.txt'

#Default maximum number of alignments to keep per read to keep the
#data array smaller
alignlimit = 200

#Default number of characters to remove from beginning of Read ID
#to make the data array smaller
nbchar = 0

#Maximum distance in bases to consider U-turns
uturndist = 49

#Default scoring options
mismatchweight = 1.25
gapweight = 1.5
aligngapweight = 1.0
samechrbonus = 0.0
threshold2 = 80.0
threshold3 = 80.0

helpfile=('Input a BLAST+ output file in tabular format to identify\n'+
          ' any rearrangement junctions present in the aligned reads.'+
          '\nAnalyzeBLAST will create 5 output files:\n- A single '+
          'junctions file containing BLAST+ lines,\n  a score, and '+
          'an analysis of rearrangements.\n- A double junctions file '+
          'containing BLAST+ lines,\n  a score, and an analysis of '+
          'rearrangements.\n- A file for unlabelled reads containing '+
          'their\n  identifier and the best score obtained for each '+
          'category.\n- A summary file giving only the result and '+
          'score for each read.\n- A statistics file containing the '+
          'number of each type\n  of rearrangement in each case.\n\n'+
          'Options:\n\n-h\t      Displays the help file.\n-i '+
          '<filename> Provide the name of the BLAST+ file to analyze.'+
          '\n-n <integer>  Remove a number of non-unique characters '+
          'from the beginning\n\t      of read identifiers. Names '+
          'remaining over 40 characters\n\t      will be cut from '+
          'the end, which may create problems.\n\t      '+
          'Default = %d' % nbchar + '\n-L <integer>  Maximum number '+
          'of BLAST+ alignments to consider for each read.\n          '+
          '    Default = %d' % alignlimit + '\n-m <float>    Penalty '+
          ' weight of each mismatch to determine the score.\n\t      '+
          'Default = %.2f' % mismatchweight + '\n-g <float>    '+
          'Penalty weight of each gap to determine the score.\n\t'+
          '      Default = %.2f' % gapweight + '\n-c <float>    Bonus '+
          'given to two reads for being on the same chromosome.\n\t'+
          '      Default = %.2f' % samechrbonus + '\n-t <float>    '+
          'Minimum score required for a read to be assigned.\n\t'+
          '      Single junctions default '+
          '= %.2f' % threshold2 + '\n\t      Double junctions '+
          'default = %.2f' % threshold3)

try:
    opts,args = getopt.getopt(sys.argv[1:],'hi:m:g:c:t:L:n:')
except getopt.GetoptError:
    print 'Invalid arguments\nUse -h for help'
    sys.exit()
for opt,arg in opts:
    if opt == '-h':
        print helpfile
        sys.exit()
    elif opt in ("-i"):
        blastfile = arg
    elif opt in ("-m"):
        mismatchweight = float(arg)
    elif opt in ("-g"):
        gapweight = float(arg)
    elif opt in ("-c"):
        samechrbonus = float(arg)
    elif opt in ("-t"):
        threshold2 = float(arg)
        threshold3 = float(arg)
    elif opt in ("-L"):
        alignlimit = int(arg)
    elif opt in ("-n"):
        nbchar = int(arg)

#--------------------------------------------------

#Count the number of lines in the file.

f1 = open(blastfile,'r')
for a,b in enumerate(f1):
    pass
lines = a+1
f1.close()

#Count the number of unique read identifiers in the file,
#find the highest number of alignments for a single read
#and find the maximum read length.

f2 = open(blastfile,'r')
prevline = f2.readline()
prevline = prevline.strip().split()
reads = 1
highestnumber = 1
readlen = int(prevline[7])
runningtally = 1
n=int((lines-1)/1000000) #Split the loop to avoid memory errors
rest=(lines-1)-(1000000*n)

for h in range(n):
    for i in range(1000000):
        newline = f2.readline()
        newline = newline.strip().split()
        if newline[0] != prevline[0]:
            reads += 1
            if runningtally > highestnumber:
                highestnumber = runningtally
            runningtally = 1
        else:
            runningtally += 1
            if runningtally > highestnumber:
                highestnumber = runningtally
        if int(newline[7]) > readlen:
            readlen = int(newline[7])
        prevline = newline

for i in range(rest):
    newline = f2.readline()
    newline = newline.strip().split()
    if newline[0] != prevline[0]:
        reads += 1
        if runningtally > highestnumber:
            highestnumber = runningtally
        runningtally = 1
    else:
        runningtally += 1
        if runningtally > highestnumber:
            highestnumber = runningtally
    if int(newline[7]) > readlen:
        readlen = int(newline[7])
    prevline = newline
f2.close()

#Limit the size of the array by removing alignments
if highestnumber >= alignlimit and alignlimit != 0:
    highestnumber = alignlimit

#Open results output file and write headers
g = open(summaryfile,'w')
g.write('Read Identifier\tResult\tScore\tSingle Read Score\t'+
        'Two Read Score\tThree Read Score\n')
h1 = open(singjuncfile,'w')
h1.write('Read Identifier\tChr\tId%\tLength\tMismatch\tGap\tReadStart'+
        '\tReadEnd\tChrStart\tChrEnd\tEvalue\tBlast Score\tScore' +
        '\tMicrohomology\tNo Microhomology\tU-Turn\n')
h2 = open(doubjuncfile,'w')
h2.write('Read Identifier\tChr\tId%\tLength\tMismatch\tGap\tReadStart'+
        '\tReadEnd\tChrStart\tChrEnd\tEvalue\tBlast Score\tScore' +
        '\tMicrohomology\tNo Microhomology\tU-Turn\n')
h3 = open(unlabelled,'w')
h3.write('Read Identifier\tSingle Read Score\tTwo Read Score'+
         '\tThree Read Score\n')

#Initialize rearrangement counts
uturn2 = 0
uturn3 = 0
micro2 = 0
micro3 = 0
nomicro2 = 0
nomicro3 = 0

#Open the data file to load alignments for each read into an array
f3 = open(blastfile,'r')
prevline = f3.readline()
prevline = prevline.strip().split()
prevline[0] = prevline[0][nbchar:]
alignid = 0
linecount = 0

for j in range(reads):
    data = numpy.zeros((highestnumber,12),dtype = 'a40')
    while linecount < lines:
        linecount += 1
        if alignid < highestnumber:
            data[alignid] = prevline
        newline = f3.readline()
        newline = newline.strip().split()
        if len(newline) > 0:
            newline[0] = newline[0][nbchar:]
            if newline[0] != prevline[0]:
                alignid = 0
                prevline = newline
                break
            elif alignid < highestnumber:
                alignid += 1
                prevline = newline
    #Initialize start position, end position and score variables
    posmin = readlen
    posmax = 1
    highscore1 = 0
    highscore2 = 0
    highscore3 = 0
    #Initialize read result flags
    singlereadalign = False
    tworeadalign = False
    threereadalign = False
    #Determine the overall start and end positions for the read
    for k in range(highestnumber):
        if data[k,6] != '':
            if int(data[k,6]) < posmin:
                posmin = int(data[k,6])
            if int(data[k,7]) > posmax:
                posmax = int(data[k,7])
            readlength = posmax-posmin+1.0
    #Determine whether the read is aligned in a single position
    for k in range(highestnumber):
        if data[k,6] != '':
            alignstart = int(data[k,6])
            alignend = int(data[k,7])
            alignlength = alignend-alignstart+1.0
            penalty = ((mismatchweight*int(data[k,4]))+
                       (gapweight*int(data[k,5])))
            onereadscore = ((100*alignlength/readlength)-penalty)
            if onereadscore > highscore1:
                highscore1 = onereadscore
            if (posmax-posmin)-(alignend-alignstart) <= 5:
                singlereadalign = True
                g.write(str(data[k,0])+'\tSingle Read Match\t'+str(highscore1)
                        +'\t'+str(highscore1)+'\n')
                break
    #Determine whether the read is aligned in two segments
    if singlereadalign == False:
        tworeadscore = numpy.zeros((highestnumber,highestnumber))
        for k in range(highestnumber-1):
            if data[k+1,0] == '':
                break
            start1 = int(data[k,6])
            end1 = int(data[k,7])
            for L in range(k+1,highestnumber,1):
                if data[L,0] == '':
                    break
                start2 = int(data[L,6])
                end2 = int(data[L,7])
                #Determine the number of bases found in the 
                #combination of alignments and number of gaps
                #Determine the overlap between pairs of reads
                coveredbases = 0.0
                overlap = 0
                aligngap = 0
                for n in range(posmin,posmax+1,1):
		    if ((n >= start1 and n <= end1) and
                            (n >= start2 and n <= end2)):
			overlap += 1
			coveredbases += 1
                    elif ((n >= start1 and n <= end1) or
                            (n >= start2 and n <= end2)):
                        coveredbases += 1
                    elif coveredbases >= 1:
                        aligngap += 1
                #Subtract any missing bases at the end of the 
                #alignments from aligngap penalty and apply the
                #appropriate weight to the penalty.
                aligngap = aligngap-(posmax-max(end1,end2))
                aligngap = aligngapweight*aligngap
                #Determine if the alignments are on the same 
                #chromosome to apply the given bonus
                if (data[k,1] == data[L,1]):
                    samechr = samechrbonus
                else:
                    samechr = 0
                #Determine if the alignments are very similar to 
                #apply a large penalty (If the overlap between reads
                #is similar to the length of one of the reads.)
                if overlap >= (end1-start1)-4 or overlap >= (end2-start2)-4:
                    samealignpenalty = 100
                else:
                    samealignpenalty = 0
                #Calculate the score
                scorelen = 100*coveredbases/readlength
                mismatchpenalty = mismatchweight*(int(data[k,4])+
                                                  int(data[L,4]))
                gappenalty = gapweight*(int(data[k,5])+int(data[L,5]))
                tworeadscore[k,L] = (scorelen+samechr-
                                     (mismatchpenalty+gappenalty+
                                      aligngap+samealignpenalty))
        #Check which pair has the highest score
        align1 = 0
        align2 = 0
        for k in range(highestnumber-1):
            for L in range(k+1,highestnumber,1):
                if tworeadscore[k,L] > highscore2:
                    highscore2 = tworeadscore[k,L]
                    #Sort the alignments in the order in which 
                    #they appear on the read
                    if int(data[k,6]) < int(data[L,6]):
                        align1 = k
                        align2 = L
                    else:
                        align1 = L
                        align2 = k
        #Write read pairs that pass the given threshold to the 
        #result and rearrangement files. Score must be higher than single
        #read score to be considered.
        if highscore2 >= threshold2:
            #Identify the type of rearrangement
            if (int(data[align1,7])-int(data[align2,6])) >= 4:
                micro2 += 1
                ismicro = True
            else:
                nomicro2 += 1
                ismicro = False
            if ((int(data[align1,9])-int(data[align1,8]))*
                (int(data[align2,9])-int(data[align2,8])) < 0
                and (abs(int(data[align2,8])-int(data[align1,9]))
                - (int(data[align1,7])-int(data[align2,6]))) <= uturndist
                and data[align1,1]==data[align2,1]):
                        uturn2 += 1
                        isuturn = True
            else:
                        isuturn = False
            tworeadalign = True
            g.write(str(data[align1,0])+'\tTwo Read Match\t'+str(highscore2)+
                    '\t'+str(highscore1)+'\t'+str(highscore2)+'\n')
            for i in range(0,12,1):
                h1.write(str(data[align1,i])+ '\t')
            h1.write(str(highscore2)+'\t'+str(ismicro)+'\t'+str(not(ismicro))
                     +'\t'+str(isuturn)+'\n')
            for i in range(0,12,1):
                h1.write(str(data[align2,i])+ '\t')
            h1.write(str(highscore2)+'\n')
    if singlereadalign == False and tworeadalign == False:
        threereadscore = numpy.zeros((highestnumber,highestnumber,
                                      highestnumber))
        for k in range(highestnumber-2):
            if data[k+2,0] == '':
                break
            start1 = int(data[k,6])
            end1 = int(data[k,7])
            for L in range(k+1,highestnumber-1,1):
                if data[L+1,0] == '':
                    break
                start2 = int(data[L,6])
                end2 = int(data[L,7])
                for m in range(L+1,highestnumber,1):
                    if data[m,0] == '':
                        break
                    start3 = int(data[m,6])
                    end3 = int(data[m,7])
                    #Determine the number of bases found in the 
                    #combination of alignments and number of gaps
                    #Determine the overlap between each pairs of reads
                    coveredbases = 0.0
                    overlap12 = 0
                    overlap23 = 0
                    overlap13 = 0
                    aligngap = 0
                    for n in range(posmin,posmax+1,1):
                        if ((n >= start1 and n <= end1) or
                                (n >= start2 and n <= end2) or
                                (n >= start3 and n <= end3)):
                            coveredbases += 1
                        elif coveredbases >= 1:
                            aligngap += 1
                        if ((n >= start1 and n <= end1) and
                                (n >= start2 and n <= end2)):
		            overlap12 += 1
		        if ((n >= start3 and n <= end3) and
                                (n >= start2 and n <= end2)):
		            overlap23 += 1
		        if ((n >= start1 and n <= end1) and
                                (n >= start3 and n <= end3)):
		            overlap13 += 1
                    #Subtract any missing bases at the end of the 
                    #alignments from aligngap penalty and apply the
                    #appropriate weight to the penalty.
                    aligngap = aligngap-(posmax-max(end1,end2,end3))
                    aligngap = aligngapweight*aligngap
                    #Determine if the alignments are on the same
                    #chromosome to apply the given bonus
                    if data[k,1] == data[L,1] == data[m,1]:
                        samechr = samechrbonus
                    else:
                        samechr = 0
                    #Determine if the alignments are very similar to 
                    #apply a large penalty (If the overlap between reads
                    #is similar to the length of one of the reads.)
                    if overlap12 >= (end1-start1)-4 or overlap12 >= (end2-start2)-4:
                        samealignpenalty = 100
                    elif overlap23 >= (end3-start3)-4 or overlap23 >= (end2-start2)-4:
                        samealignpenalty = 100
                    elif overlap13 >= (end1-start1)-4 or overlap13 >= (end3-start3)-4:
                        samealignpenalty = 100
                    else:
                        samealignpenalty = 0
                    scorelen = 100*coveredbases/readlength
                    mismatchpenalty = mismatchweight*(int(data[k,4])+
                                                      int(data[L,4])+
                                                      int(data[m,4]))
                    gappenalty = gapweight*(int(data[k,5])+int(data[L,5])+
                                            int(data[m,5]))
                    penalty = mismatchpenalty+gappenalty+aligngap+samealignpenalty
                    threereadscore[k,L,m] = scorelen+samechr-penalty
        #Check which pair has the highest score
        align1 = 0
        align2 = 0
        align3 = 0
        for k in range(highestnumber-2):
            for L in range(k+1,highestnumber-1,1):
                for m in range(L+1,highestnumber,1):
                    if threereadscore[k,L,m] > highscore3:
                        highscore3 = threereadscore[k,L,m]
                        #Sort the alignments in the order in which
                        #they appear on the read
                        beforesort = [(k,int(data[k,6])),(L,int(data[L,6])),
                                      (m,int(data[m,6]))]
                        aftersort = sorted(beforesort,key=itemgetter(1))
                        align1 = aftersort[0][0]
                        align2 = aftersort[1][0]
                        align3 = aftersort[2][0]
        #Write read pairs that pass the given threshold to the result
        #and rearrangement files. Score must be higher than single
        #read score to be considered.
        if highscore3 >= threshold3:
            threereadalign = True
            #Identify the type of rearrangement for the first junction
            if (int(data[align1,7])-int(data[align2,6])) >= 4:
                micro3 += 1
                ismicro = True
            else:
                nomicro3 += 1
                ismicro = False
            if ((int(data[align1,9])-int(data[align1,8]))*
                (int(data[align2,9])-int(data[align2,8])) < 0
                and (abs(int(data[align2,8])-int(data[align1,9]))
                - (int(data[align1,7])-int(data[align2,6]))) <= uturndist
                and data[align1,1]==data[align2,1]):
                        uturn3 += 1
                        isuturn = True
            else:
                        isuturn = False
            g.write(str(data[align1,0])+'\tThree Read Match\t'+str(highscore3)
                    +'\t'+str(highscore1)+'\t'+str(highscore2)+'\t'+
                    str(highscore3)+'\n')
            for i in range(0,12,1):
                h2.write(str(data[align1,i])+ '\t')
            h2.write(str(highscore3)+'\t'+str(ismicro)+'\t'+str(not(ismicro))+
                     '\t'+str(isuturn)+'\n')
            #Identify the type of rearrangement for the second junction
            if (int(data[align2,7])-int(data[align3,6])) >= 4:
                micro3 += 1
                ismicro = True
            else:
                nomicro3 += 1
                ismicro = False
            if ((int(data[align2,9])-int(data[align2,8]))*
                (int(data[align3,9])-int(data[align3,8])) < 0
                and (abs(int(data[align3,8])-int(data[align2,9])) -
                (int(data[align2,7])-int(data[align3,6]))) <= uturndist
                and data[align2,1]==data[align3,1]):
                        uturn3 += 1
                        isuturn = True
            else:
                        isuturn = False
            for i in range(0,12,1):
                h2.write(str(data[align2,i])+ '\t')
            h2.write(str(highscore3)+'\t'+str(ismicro)+'\t'+str(not(ismicro))
                     +'\t'+str(isuturn)+'\n')
            for i in range(0,12,1):
                h2.write(str(data[align3,i])+ '\t')
            h2.write(str(highscore3)+'\n')   
    if singlereadalign == tworeadalign == threereadalign == False:
        g.write(str(data[align1,0])+ '\tNo Good Match\n')
        h3.write(str(data[align1,0])+'\t'+str(highscore1)
                 +'\t'+str(highscore2)+'\t'+str(highscore3)+'\n')

f3.close()
g.close()
h1.close()
h2.close()
h3.close()

#Calculate the final stats and write them to a file
micro = micro2+micro3
nomicro = nomicro2+nomicro3
uturn = uturn2+uturn3
g2 = open(statfile, 'w')
g2.write('\tSingle Junctions\tDouble Junctions\tTotal\nMicrohomology\t')
g2.write(str(micro2)+'\t'+str(micro3)+'\t'+str(micro)+'\n')
g2.write('No Microhomology\t'+str(nomicro2)+'\t'+str(nomicro3)+'\t')
g2.write(str(nomicro)+'\nU-Turn\t'+str(uturn2)+'\t'+str(uturn3)+'\t')
g2.write(str(uturn))
g2.close()

