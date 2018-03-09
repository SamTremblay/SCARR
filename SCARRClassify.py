## SCARRClassify (v1.0)
## Copyright 2018, Samuel Tremblay-Belzile

import numpy
import sys
import os
import getopt
from operator import itemgetter

#Optional run name for input file format "SingleJunctions[runname].txt"
runname=''

try:
    opts,args = getopt.getopt(sys.argv[1:],'n:')
except getopt.GetoptError:
    print 'Invalid arguments'
    sys.exit()
for opt,arg in opts:
    if opt == '-n':
        runname = arg

singjuncfile = 'SingleJunctions%s.txt' % (runname)
singoutput = 'SingleJunctions%sDetail.txt' % (runname)

#Count the number of lines in the files.

f1 = open(singjuncfile,'r')
for a,b in enumerate(f1):
    pass
lines1 = a+1
singjuncs=(lines1-1)/2
f1.close()

#Initialize the arrays for the rearrangement data
data1=numpy.zeros(7,dtype='int64')
data2=numpy.zeros(7,dtype='int64')

f1 = open(singjuncfile, 'r')
g1 = open(singoutput, 'w')

#Create temporary files for deletions, duplications and inversions
h1 = open("DeletionsTemp.txt", 'w')
h2 = open("DuplicationsTemp.txt", 'w')
h3 = open("InversionsTemp.txt", 'w')
#Initialize rearrangement count
nbdup = 0
nbdel = 0
nbinv = 0

headers=f1.readline().strip('\n') #Copy headers
g1.write(headers + '\tSlippage\tOverlap\tInsertion\tInversion\tDeletion'
         '\tDuplication\tTranslocation\tDistance\tBreakpoint 1\tBreakpoint 2\n')

for i in range(singjuncs):
    line1=f1.readline().strip('\n')
    line2=f1.readline().strip('\n')
    align1=line1.split('\t')
    align2=line2.split('\t')
    for z in range(7):
        data1[z]=int(align1[z+3])
        data2[z]=int(align2[z+3])
    if (data1[4]-data2[3]) >= 0:
        overlap = data1[4]-data2[3]+1
        insertion = 0
    else:
        overlap = 0
        insertion = abs(data1[4]-data2[3])-1
    if align1[1] == align2[1]:
        samechr=True
        if ((data1[6]-data1[5])*(data2[6]-data2[5]))>0:
            deldup=True
            inversion=False
        else:
            inversion=True
            deldup=False
    else:
        samechr=False
        inversion=False
        deldup=False
    if inversion:
        nbinv += 1
        if (data1[5]) < (data1[6]):
            breakpoint1 = data1[6]
            breakpoint2 = data2[5]
        else:
            breakpoint1 = data2[5]
            breakpoint2 = data1[6]
        distance = abs(breakpoint2-breakpoint1)-overlap
        if distance < 0:
            distance = 0
        if data1[6]-data1[5] > 0:
            h3.write(align1[0]+'\t'+align1[1]+'\t'+str(breakpoint1)+'\t'+
                     str(breakpoint2)+'\t'+str(distance)+
                     '\t'+str(overlap)+'\t'+str(insertion)+'\n')
        else:
            h3.write(align1[0]+'\t'+align1[1]+'\t'+str(breakpoint2)+'\t'+
                     str(breakpoint1)+'\t'+str(distance)+
                     '\t'+str(overlap)+'\t'+str(insertion)+'\n')
    if deldup:
        if (data1[6]-data1[5])>0:
            forward=True
            distance=(data2[5]-data1[6])+overlap-insertion-1
        else:
            forward=False
            distance=(data1[6]-data2[5])+overlap-insertion-1
        if distance > insertion:
            deletion=True
            duplication=False
        else:
            deletion=False
            duplication=True
        if duplication:
            nbdup += 1
            if forward:
                if overlap == 0:
                    breakpoint2 = data1[6]
                    breakpoint1 = data2[5]
                else:
                    if (data1[1]+data1[2]) < (data2[1]+data2[2]):
                        breakpoint2 = data1[6] - overlap
                        breakpoint1 = data2[5]
                    else:
                        breakpoint2 = data1[6]
                        breakpoint1 = data2[5] + overlap
            else:    
                if overlap == 0:
                    breakpoint1 = data1[6]
                    breakpoint2 = data2[5]
                else:
                    if (data1[1]+data1[2]) <= (data2[1]+data2[2]):
                        breakpoint1 = data1[6]
                        breakpoint2 = data2[5] - overlap
                    else:
                        breakpoint1 = data1[6] + overlap
                        breakpoint2 = data2[5]
            h2.write(align1[0]+'\t'+align1[1]+'\t'+str(breakpoint1)+'\t'+
                     str(breakpoint2)+'\t'+str(distance)+
                     '\t'+str(overlap)+'\t'+str(insertion)+'\n')
        if deletion:
            nbdel += 1
            if forward:
                if overlap == 0:
                    breakpoint1 = data1[6] + 1
                    breakpoint2 = data2[5] - 1
                else:
                    if (data1[1]+data1[2]) < (data2[1]+data2[2]):
                        breakpoint1 = data1[6] + 1
                        breakpoint2 = data2[5] + overlap - 1
                    else:
                        breakpoint1 = data1[6] - overlap + 1
                        breakpoint2 = data2[5] - 1
            else:    
                if overlap == 0:
                    breakpoint1 = data2[5] + 1
                    breakpoint2 = data1[6] - 1
                else:
                    if (data1[1]+data1[2]) <= (data2[1]+data2[2]):
                        breakpoint1 = data2[5] - overlap + 1
                        breakpoint2 = data1[6] - 1
                    else:
                        breakpoint1 = data2[5] + 1
                        breakpoint2 = data1[6] + overlap - 1
            h1.write(align1[0]+'\t'+align1[1]+'\t'+str(breakpoint1)+'\t'+
                     str(breakpoint2)+'\t'+str(distance)+
                     '\t'+str(overlap)+'\t'+str(insertion)+'\n')
        if abs(distance)<50:
            slip50=True
        else:
            slip50=False
    else:
        deletion=False
        duplication=False
        slip50=False
    g1.write(line1 + '\t' + str(slip50) + '\t' + str(overlap) + '\t' +
             str(insertion) + '\t' + str(inversion) + '\t' + str(deletion)
             + '\t' + str(duplication) + '\t' + str(not(samechr)) + '\t')
    if deldup or inversion:
        g1.write(str(distance) + '\t' + str(breakpoint1) + '\t')
        g1.write(str(breakpoint2) + '\n')
    else:
        g1.write('\t\t\n')
    g1.write(line2 + '\n')

#Close rearrangement and output files
f1.close()
g1.close()

#Close temporary files
h1.close()
h2.close()
h3.close()

#Open temporary files to read, group and sort
h1 = open("DeletionsTemp.txt", 'r')
h2 = open("DuplicationsTemp.txt", 'r')
h3 = open("InversionsTemp.txt", 'r')
#Initialize data arrays for rearrangements
dtype1=[('id', 'a40'),('chr', 'a40'),('bp1', 'int64'),('bp2', 'int64'),
       ('dist', 'int64'),('over', 'int32'),('ins', 'int32')]
datadel=numpy.zeros(nbdel,dtype=dtype1)
datadup=numpy.zeros(nbdup,dtype=dtype1)
datainv=numpy.zeros(nbinv,dtype=dtype1)
#Load the data into the arrays
for i in range(nbdel):
    words=h1.readline().strip('\n').split('\t')
    datadel[i]=(words[0],words[1],int(words[2]),int(words[3]),int(words[4]),
                int(words[5]),int(words[6]))
for i in range(nbdup):
    words=h2.readline().strip('\n').split('\t')
    datadup[i]=(words[0],words[1],int(words[2]),int(words[3]),int(words[4]),
                int(words[5]),int(words[6]))
for i in range(nbinv):
    words=h3.readline().strip('\n').split('\t')
    datainv[i]=(words[0],words[1],int(words[2]),int(words[3]),int(words[4]),
                int(words[5]),int(words[6]))
#Close temporary files
h1.close()
h2.close()
h3.close()
    
#Sort the rearrangements by chromosome and position of first breakpoint
datadel=numpy.sort(datadel,order=['chr','bp1','bp2','over','ins'])
datadup=numpy.sort(datadup,order=['chr','bp1','bp2','over','ins'])
datainv=numpy.sort(datainv,order=['chr','bp1','bp2','over','ins'])

#Create rearrangement files and write column headers
g1=open("Deletions.txt",'w')
g2=open("Duplications.txt",'w')
g3=open("Inversions.txt",'w')
g1.write("Identifier\tChromosome\tBreakpoint 1\tBreakpoint 2\tDistance\t"+
         "Microhomology Length\tInserted Bases\n")
g2.write("Identifier\tChromosome\tBreakpoint 1\tBreakpoint 2\tDistance\t"+
         "Microhomology Length\tInserted Bases\n")
g3.write("Identifier\tChromosome\tBreakpoint 1\tBreakpoint 2\tDistance\t"+
         "Microhomology Length\tInserted Bases\n")
#Write the sorted Deletions, Duplications and Inversions to files
for i in range(nbdel):
    g1.write(datadel[i][0]+'\t'+datadel[i][1]+'\t'+str(datadel[i][2])+
             '\t'+str(datadel[i][3])+'\t'+str(datadel[i][4])+
             '\t'+str(datadel[i][5])+'\t'+str(datadel[i][6])+'\n')
for i in range(nbdup):
    g2.write(datadup[i][0]+'\t'+datadup[i][1]+'\t'+str(datadup[i][2])+
             '\t'+str(datadup[i][3])+'\t'+str(datadup[i][4])+
             '\t'+str(datadup[i][5])+'\t'+str(datadup[i][6])+'\n')
for i in range(nbinv):
    g3.write(datainv[i][0]+'\t'+datainv[i][1]+'\t'+str(datainv[i][2])+
             '\t'+str(datainv[i][3])+'\t'+str(datainv[i][4])+
             '\t'+str(datainv[i][5])+'\t'+str(datainv[i][6])+'\n')

g1.close()
g2.close()
g3.close()

#Create arrays for rearrangement groups
dtype2=[('chr', 'a40'),('bp1', 'int64'),('bp2', 'int64'),
        ('dist', 'int64'),('nb', 'int32'),('mh','int32'),
        ('ins','int32'),('id', 'a500')]
sorteddel=numpy.zeros(nbdel,dtype=dtype2)
sorteddup=numpy.zeros(nbdup,dtype=dtype2)
sortedinv=numpy.zeros(nbinv,dtype=dtype2)

#Group identical deletions
prevchr=datadel[0][1]
prevbp1=datadel[0][2]
prevbp2=datadel[0][3]
prevdist=abs(prevbp2-prevbp1)+1
prevmh=datadel[0][5]
previns=datadel[0][6]
nb = 1
identifiers = datadel[0][0]
groupsdel=1
for i in range(1, nbdel):
    if (datadel[i][1]==prevchr and datadel[i][2]==prevbp1 and
        datadel[i][3]==prevbp2 and datadel[i][5]==prevmh and
        datadel[i][6]==previns):
        nb += 1
        identifiers += ","
        identifiers += datadel[i][0]
    else:
        sorteddel[i-1]=(prevchr, prevbp1, prevbp2, prevdist, nb, prevmh,
                        previns, identifiers)
        prevchr=datadel[i][1]
        prevbp1=datadel[i][2]
        prevbp2=datadel[i][3]
        prevdist=abs(prevbp2-prevbp1)+1
        prevmh=datadel[i][5]
        previns=datadel[i][6]
        nb = 1
        identifiers = datadel[i][0]
        groupsdel += 1
sorteddel[nbdel-1]=(prevchr, prevbp1, prevbp2, prevdist, nb, prevmh,
                    previns, identifiers)
    
#Group identical duplications
prevchr=datadup[0][1]
prevbp1=datadup[0][2]
prevbp2=datadup[0][3]
prevdist=abs(prevbp2-prevbp1)+1
prevmh=datadup[0][5]
previns=datadup[0][6]
nb = 1
identifiers = datadup[0][0]
groupsdup=1
for i in range(1, nbdup):
    if (datadup[i][1]==prevchr and datadup[i][2]==prevbp1 and
        datadup[i][3]==prevbp2 and datadup[i][5]==prevmh and
        datadup[i][6]==previns):
        nb += 1
        identifiers += ","
        identifiers += datadup[i][0]
    else:
        sorteddup[i-1]=(prevchr, prevbp1, prevbp2, prevdist, nb, prevmh,
                        previns, identifiers)
        prevchr=datadup[i][1]
        prevbp1=datadup[i][2]
        prevbp2=datadup[i][3]
        prevdist=abs(prevbp2-prevbp1)+1
        prevmh=datadup[i][5]
        previns=datadup[i][6]
        nb = 1
        identifiers = datadup[i][0]
        groupsdup += 1
sorteddup[nbdup-1]=(prevchr, prevbp1, prevbp2, prevdist, nb, prevmh,
                    previns, identifiers)

#Group identical inversions
prevchr=datainv[0][1]
prevbp1=datainv[0][2]
prevbp2=datainv[0][3]
prevmh=datainv[0][5]
prevdist=abs(prevbp2-prevbp1)-prevmh
if prevdist < 0:
    prevdist = 0
previns=datainv[0][6]
nb = 1
identifiers = datainv[0][0]
groupsinv=1
for i in range(1, nbinv):
    if (datainv[i][1]==prevchr and datainv[i][2]==prevbp1 and
        datainv[i][3]==prevbp2 and datainv[i][5]==prevmh and
        datainv[i][6]==previns):
        nb += 1
        identifiers += ","
        identifiers += datainv[i][0]
    else:
        sortedinv[i-1]=(prevchr, prevbp1, prevbp2, prevdist, nb, prevmh,
                        previns, identifiers)
        prevchr=datainv[i][1]
        prevbp1=datainv[i][2]
        prevbp2=datainv[i][3]
        prevmh=datainv[i][5]
        prevdist=abs(prevbp2-prevbp1)-prevmh
        if prevdist < 0:
            prevdist = 0
        previns=datainv[i][6]
        nb = 1
        identifiers = datainv[i][0]
        groupsinv += 1
sortedinv[nbinv-1]=(prevchr, prevbp1, prevbp2, prevdist, nb, prevmh,
                    previns, identifiers)

#Sort grouped rearrangements by their frequency and reverse the array to
#make the order descending
sorteddel=numpy.sort(sorteddel,order='nb')
desdel=sorteddel[::-1]
sorteddup=numpy.sort(sorteddup,order='nb')
desdup=sorteddup[::-1]
sortedinv=numpy.sort(sortedinv,order='nb')
desinv=sortedinv[::-1]

#Create grouped rearrangement files and write column headers
g1=open("DeletionsGrouped.txt",'w')
g2=open("DuplicationsGrouped.txt",'w')
g3=open("InversionsGrouped.txt",'w')
g1.write("Chromosome\tBreakpoint 1\tBreakpoint 2\tLength\tFrequency\t"+
         "Microhomology Length\tInserted Bases\tIdentifiers\n")
g2.write("Chromosome\tBreakpoint 1\tBreakpoint 2\tLength\tFrequency\t"+
         "Microhomology Length\tInserted Bases\tIdentifiers\n")
g3.write("Chromosome\tBreakpoint 1\tBreakpoint 2\tLength\tFrequency\t"+
         "Microhomology Length\tInserted Bases\tIdentifiers\n")
#Write the grouped Deletions, Duplications and Inversions to files
for i in range(groupsdel):
    g1.write(desdel[i][0]+'\t'+str(desdel[i][1])+'\t'+str(desdel[i][2])+
             '\t'+str(desdel[i][3])+'\t'+str(desdel[i][4])+
             '\t'+str(desdel[i][5])+'\t'+str(desdel[i][6])+
             '\t'+desdel[i][7]+'\n')
for i in range(groupsdup):
    g2.write(desdup[i][0]+'\t'+str(desdup[i][1])+'\t'+str(desdup[i][2])+
             '\t'+str(desdup[i][3])+'\t'+str(desdup[i][4])+
             '\t'+str(desdup[i][5])+'\t'+str(desdup[i][6])+
             '\t'+desdup[i][7]+'\n')
for i in range(groupsinv):
    g3.write(desinv[i][0]+'\t'+str(desinv[i][1])+'\t'+str(desinv[i][2])+
             '\t'+str(desinv[i][3])+'\t'+str(desinv[i][4])+
             '\t'+str(desinv[i][5])+'\t'+str(desinv[i][6])+
             '\t'+desinv[i][7]+'\n')


#Remove temporary files for deletions, duplications and inversions
os.remove("DeletionsTemp.txt")
os.remove("DuplicationsTemp.txt")
os.remove("InversionsTemp.txt")
