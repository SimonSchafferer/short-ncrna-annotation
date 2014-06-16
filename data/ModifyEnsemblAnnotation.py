'''
Created on 13.06.2014

@author: schaffrr
'''

import csv
import sys
import getopt

class EnsemblAnnotation():
        
    def __init__(self, inputFile):
        self._inputFile = inputFile
    
    def readAndModifyEnsembleAnnotation(self, ouputFile):
        r = csv.reader(open(self._inputFile), delimiter="\t") # Here your csv file
        writer = csv.writer(open(ouputFile, "wb"), delimiter = '\t')
        
        firstRow = r.next() #header line
#         print firstRow
        transcriptID_idx = firstRow.index("transcript_id")
        start_idx = firstRow.index("start")
        end_idx = firstRow.index("end")
        strand_idx = firstRow.index("strand")
        width_idx = firstRow.index("width")
        #optional argument therefore try to find it otherwise add it for comfort
        exon_number_present = True
        try:
            exon_number_idx = firstRow.index("exon_number")
        except ValueError:
            exon_number_present = False
            print "adding exon_number column"
            firstRow.append("exon_number")
            exon_number_idx = firstRow.index("exon_number")
            
        type_idx = firstRow.index("type")
        previous = ""
        
        writer.writerow(firstRow)
        insertCount = 0
        exonCount = 1
        for row in r:
            if row:
                if not exon_number_present:
                    row.append(1)
                if row[type_idx] == "exon":
                    current = row
                    if previous == "":
                        previous = current
                
                    if current[transcriptID_idx] == previous[transcriptID_idx]:
                        exonCount = exonCount + 1
                        if not exon_number_present:
                            row[exon_number_idx] = exonCount
                        #checking strandness and adjusting insertStart and insertEnd
                        if previous[strand_idx] == "+":
                            insertStart = int(previous[end_idx]) + 1
                            insertEnd = int(current[start_idx])-1
                        elif previous[strand_idx] == "-":
                            insertEnd = int(previous[start_idx]) - 1
                            insertStart = int(current[end_idx]) + 1
                        
                        if insertStart <= insertEnd:  
                            insertCount = insertCount + 1
                            insert = list(current)
                            toChangeIdx = [ start_idx, end_idx, width_idx, type_idx, exon_number_idx]
                            toChangeVect = [insertStart, insertEnd, insertEnd-insertStart+1, "intron",insertCount]
                            for i,j in zip(toChangeIdx,range(0,len(toChangeVect))):
                                insert[i] = toChangeVect[j]
                            writer.writerow(insert)
                    else:
                        insertCount = 0
                        exonCount = 0
                    previous = current
                writer.writerow(row)


def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["help","ifile=", "ofile="])
    except getopt.GetoptError:
        print 'ModifyEnsemblAnnotation.py -i <inputfile> -o <outputfile>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'help'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
    if inputfile == "":
        print "input filename must not be empty! see -h for parameters"
        sys.exit()
    if outputfile == "":
        print "no output filename, choosing <input>.mod"
        outputfile = inputfile + ".mod"

    modifyEnsemblAnnotation = EnsemblAnnotation(inputfile)
    print "calling readAndModify"
    modifyEnsemblAnnotation.readAndModifyEnsembleAnnotation(outputfile)
    
#     print "finished correctly " + str(finValid)
    
if __name__ == "__main__":
    main(sys.argv[1:])

