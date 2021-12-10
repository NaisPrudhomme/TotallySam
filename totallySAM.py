#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#module importation
import re, sys, os
#re : re.search function
#sys : parameters list
#os : for "checking file" part

#***********************************************************************************************
#***********************************************************************************************

#Main
def main(inputSam) :

    print(" File Name : \n",inputSam) #print file name
    
    #print(extraction(inputSam)) #Uncomment for dictionnary printing    
    
    desc(extraction(inputSam)) #Count the number of unmapped reads, mapped reads, partially mapped reads and paires of mapped/unmapped reads and paires of mapped/partially mapped reads 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
#Read file and extract reads name, FLAG, CIGAR and TAG MD:Z from SAM file and import it in a dictionnary (Key = Flag) of dictionnary (Key = read names)
def extraction(inputSam) :
    
    with open(inputSam, "r") as fSam : #Open SAM file (reading)

        print("\n Extraction and file description : Start \n")#Print the file description

        headerN = 0

        
        dicoExt = {} #create a dictionnary for extraction that will contains reads data
        readNumber = 0
        corruptedRead = 0

        #Description of @HD header's line
        HDpresent = 0
        
        #Description of @SQ header's line
        SQpresent = 0

        #Description of @PG header's line
        PGpresent = 0

        #Description of @RG header's line
        RGpresent = 0

        #Extract headers lines and reads lines
        for reads in fSam :
            descReads = re.search("^@.*", reads) #research line with @ at the beginning

            #If the line begin with "@" its a description line
            if descReads :
                headerN += 1 #Count header's line number 

                #Description of header line
                resHD = re.search("^@HD",reads) 
                if resHD :
                    HDpresent = 1
                    reshVN = re.search("(VN\:.*\t|VN\:.*\n$)", reads) #Get Version format
                    if reshVN :
                        hVN = reshVN.group(0)[3:-1]
                        hVN_split = hVN.split("\t")
                        hVN = hVN_split[0]
                    
                #Description of refence sequence line
                resSQ = re.search("^@SQ",reads)
                if resSQ :
                    SQpresent = 1
                    resSN = re.search("(SN\:.*\t|SN\:.*\n$)", reads) #Get Reference sequence name
                    if resSN :
                        sSN = resSN.group(0)[3:-1]
                        sSN_split = sSN.split("\t")
                        sSN = sSN_split[0]
                    resLN = re.search("(LN\:.*\t|LN\:.*\n$)",reads) #Get reference sequence length
                    if resLN :
                        sLN = resLN.group(0)[3:-1]
                        sLN_split = sLN.split("\t")
                        sLN = sLN_split[0]

                #Description of program line
                resPG = re.search("^@PG", reads) 
                if resPG :
                    PGpresent = 1
                    respID = re.search("(ID\:.*\t|ID\:.*\n$)",reads) #Get program identifier
                    if respID :
                        pID = respID.group(0)[3:-1]
                        pID_split = pID.split("\t")
                        pID = pID_split[0]

                #Description of read group
                resRG = re.search("^@RG",reads)
                if resRG :
                    RGpresent = 1
                    resrID = re.searche("(ID\:.*\t|ID\:.*\n$)",reads) #Read group identifier
                    if resrID :
                        rID = respID.group(0)[3:-1]
                        rID_split = rID.split("\t")
                        rID = rID_split[0]

            #If the line is a read, the program proceed to extraction 
            else :
                readNumber += 1 #Count the number of reads
                #print(reads)
            
                col = reads.split("\t") #split the read after each tabulation for
                #print(col)
                #print("column number",len(col))

                #Does my file contain the right number of column (9 tabulations and 10 columns)
                if len(col) < 10 : 
                    print("WARNING : the column number for one read is not in accord with sam standard file")
                    corruptedRead += 1
                    #if my read have less than 10 columns, the program print an WARNING message 

                #Get read informations
                flag = int(col[1])  #int() convert flag to int, get flag for each read (2nd column)
                nomReads = col[0] #get read name for each read (1st column)
                cigar = col[5] #get cigar dor each read (6th column)

                #Get TagMDZ 
                resTagMDZ = re.search("((MD\:Z\:.*?)\t)", reads) #MD:Z is a optionnal information
                #if the read contains MD:Z, get MD:Z for the read
                if resTagMDZ : 
                    resTagMDZ = resTagMDZ.group(1)[:-1]# to get "MD:Z:.*" without the tabulation
                    #tagMDZ.append(resTagMDZ)
                    tagMDZ = resTagMDZ

                    completeFill(dicoExt,flag,nomReads,cigar,tagMDZ)
                    #Call the "completeFill" function to fill the dictionnary (with MD:Z in value)

                else :
                    incompleteFill(dicoExt,flag,nomReads,cigar)
                    #Call the "incompleteFill" function to fill the dictionnary (without MD:Z)


        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print("File Description :")

        if HDpresent == 1 :
            print(" - Format version :", hVN)
        
        if SQpresent == 1 :
            print(" - Reference sequence name :",sSN)
            print(" - Reference sequence lenght :",sLN)

        if PGpresent == 1 :
            print(" - Program identifier :", pID)

        if RGpresent == 1 :
            print(" - Read group identifier :",rID)

        print("\n")
        print("File informations :")
        print(" - Number of header lines :", headerN)
        print(" - Total number of reads : ",readNumber)
        print(" - Total of corrupted reads (ignored in description and analysis) : ", corruptedRead)
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")

        print(" Extraction : Finished \n") #Print for the end of the extraction

        return dicoExt #return the dictionnary that contains flag (key for the first dictionnary), the read name (key for the second dictionnary), cigar and MD:Z for each read


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#fill function that create a dictionnary dictionnary
def completeFill(dicoExt, flag, nomReads, cigar, tagMDZ) :

    #Create a dictionnary with flag for key and the second dictionnary for value
    #The second dictionnary has name read as keay and cigar and MD:Z (tagMDZ) as values

    if flag in dicoExt : 
        dicoExt[flag][nomReads] = [cigar, tagMDZ]
        #if the flag already exist, the program creates a new key for the first dictionnary
    else :
        dicoExt[flag] = {nomReads : [cigar, tagMDZ]}
        #if the flag doesn't exist in the first dictionnary, the program creates a new key for the first dictionnary

#fill function where no MDZ is add (when the read description doesn't give MDZ information
def incompleteFill(dicoExt, flag, nomReads, cigar) :

    #Create a dictionnary with flag for key and the second dictionnary for value
    #The second dictionnary has name read as keay and cigar

    if flag in dicoExt : 
        dicoExt[flag][nomReads] = [cigar]
        #if the flag already exist, the program creates a new key for the first dictionnary
    else :
        dicoExt[flag] = {nomReads : [cigar]}
        #if the flag doesn't exist in the first dictionnary, the program creates a new key for the first dictionnary

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Description function that count the reads number for each classes
def desc(extraction) :

    print(" Reads Analysis : Start \n") # Print for the begining of the reads analysis

    dicoExt = extraction #call the return of extraction" function
    #print(dicoExt)

    #flag for each class : read description
    mapped = [67, 73, 83, 89, 97, 99, 115, 121, 131, 137, 145, 147, 153, 163, 179, 185]
    unmapped = [63, 69, 77, 101, 117, 133, 141, 165, 181]

#####
    #flag for each class : pairs description
    pairedMapped1 = [99, 83, 67, 115, 81, 97, 65, 113, 147, 163, 131, 179, 161, 145, 129, 177] #1st read puis 2nd read
    pairedMapped2 = [147, 163, 131, 179, 161, 145, 129, 177] #2nd read
    pairedMapUnmap = [73, 89, 121, 137, 153, 185] #just the unmapped read
    pairedUnmapped = [77]
                    
    #read counters
    rm = 0 #mapped perfectly
    ru = 0 #unmapped
    rp = 0 #partially mapped
    uncommon = 0 #read with uncommon flag, these reads are ignored in pair analysis

    #pair counters
    pm = 0 #pair perfectly mapped
    pmu = 0 #one mapped, one unmapped
    pmp = 0 #one mapped, one partially
    pmm = 0 #pair mapped mapped
    puu = 0 #pair unmapped unmapped
    ppu = 0 #pair partial unmapped
    
    #List of reads name for each class (for pairs counting)
    nameMappedPerfect = [] #stock perfectly mapped name
    nameMappedPartially = [] #stock partially mapped name
    nameUnmapped = [] #stock unmapped name

    #Counting loop
    for flag in dicoExt : #for each flag
        for nomReads in dicoExt[flag]: #for each read
            
            if flag in mapped :
            #if flag matches mapped score

                #Get mapped reads
                if (dicoExt[flag][nomReads][1] == 'MD:Z:100' or dicoExt[flag][nomReads][0] == '100M'):
                    rm += 1
                    nameMappedPerfect.append(nomReads)
                    #if the MD:Z score and cigar show perfect mapping :
                       #(1) count one more perfectly mapped read
                       #(2) add perfect mapped read names to the list

                #Get partially mapped reads
                else :
                    rp += 1 
                    nameMappedPartially.append([nomReads,dicoExt[flag][nomReads][1]])
                    #if alignment score if not 100 (read length) :
                       #(1) count one more partially mapped read ;
                       #(2) add read name and MD:Z score to the list (MD:Z score for variant analysis option)


                    #Get perfectly mapped/partially mapped pairs
                    if flag in pairedMapped1 and nomReads in nameMappedPerfect : 
                        pmp += 1
                        #if the flag of one of the read is in mapped read and the read name in the list which contains perfect mapped read, count one more pair "perfectly mapped/partially"

            #Get unmapped reads            
            if flag in unmapped :
                ru +=1
                nameUnmapped.append(nomReads)
                #if flag matches unmapped score : count one more unmapped read
                
            #Get read with uncommon flag (ignored for the analysis)
            if not flag in mapped and not flag in unmapped :
                uncommon += 1
                #if the flag doesn't match both mapped and unmapped score, count one more read with uncommon flag (ignored within analysis)

            #Get Mapped perfectly/unmapped pair
            if flag in pairedMapUnmap and nomReads in nameMappedPartially:
                ppu += 1 
                #if flag in map/unmapped list, count one more 

            #Get mapped perfectly/unmapped pair
            if flag in pairedMapUnmap and nomReads in nameMappedPerfect:
                pmu += 1
                #if flag in map/unmapped, count one more

            #Get both unmapped read pairs
            if flag in pairedUnmapped :
                puu +=1
                #if flag in unmapped list, count one more
    
    #This bloc print results
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("Reads description")
    print(" - Reads perfectly mapped : \t",rm)
    print(" - Reads partially mapped : \t", rp)
    print(" - Reads unmapped : \t \t",ru)
    print(" - Reads with uncommon flag : \t",uncommon)

    print("\n")
    print("Pairs description")
    print(" - Pairs partially mapped/unmapped : \t", ppu)
    print(" - Pairs mapped/unmapped : \t \t", pmu)
    print(" - pairs mapped/partially mapped : \t",pmp)# pmu + pmp
    print(" - pairs unmapped/unmapped : \t \t", puu)
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")

    print("Reads Analysis : Finished") # Print for the end of the reads analysis

    #If the user want to know if some partially mapped reads have substitutions, it call variantANalysis (the user specifies an other parameter
    if len(sys.argv) == 3 :
        if (sys.argv[2]) == "-v" :
            VariantAnalysis(nameMappedPartially)
    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def VariantAnalysis(nameMappedPartially) :

    print("\n")
    print("====================================================================")
    print("====================================================================")
    print("\n")
    print("Variant analysis")
    print("\n")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

    #If a nucleotide doesn't match, the MD:Z contains the new nucleotide
    possibilities = ("A", "T", "G", "C")

    for variantRead in nameMappedPartially : #for each read in partially mapped list

        #Get read name (vRead) and MD:Z score (vMDZ)
        vRead = variantRead[0] 
        vMDZ = variantRead[1][5:] #Get MD:Z score without "MD:Z"

        #itialise SNP at 0
        SNP = 0

        #for each character in MD:Z score, if it contains substitution (nucleotide), count one more mutation
        for i in vMDZ :
            if i in possibilities :
                SNP += 1

        #if the read contains substitution, print number of substitution per read
        if SNP != 0 :
            print("SNP for",vRead, ":", SNP)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
#Cheking file : Error handling

print(" Checking file : Start")

#Does the user specified file name ?
if len(sys.argv) >= 2 :  #(check if there is an argument after program name)

    #Does my file is empty ?
    if os.stat(sys.argv[1]).st_size > 1 : #if file size is higher than 1, continue (an empty sam file can have a size about 1)

        #Does my file is a file ? Does my file is in Sam format ?
        if os.path.isfile(sys.argv[1]) and sys.argv[1].endswith('.sam') : #if the file is a file in Sam format, continue

            #Does my Sam file has header ?
            file = open(sys.argv[1], "r") #"Open the file for checking"
            ligne = file.readline() #Read the first line of the file
            res = re.search("^@.*", ligne) #Check if there is a header
            
            if res : #If the header is in the file, print a message for processing, run main program
                print(" Checking file : End (Sam file, not empty and uncorrupted)  \n TotallySam Processing ... \n")
                file.close()
                main(sys.argv[1]) #Run main program
                

            else : #If there is no header, print WARNING message and exit program
                print(" WARNING : Corrupted file (no Sam header)")
                sys.exit()
            
        else : #If my file is not a file or not in SAM format, print error message
            print("\n WARNING : This is not a file or the file is not a SAM extension")
        

    else : #If my file is empty, print error message
            print("\n WARNING : The file is empty (file size < 1)")
        
else : #If the user don't specify file name, print error message
        print("\n WARNING : The file name should be specified on first argument (after program name)")
