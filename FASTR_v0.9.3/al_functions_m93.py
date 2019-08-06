#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 12:26:59 2019

@author: zzhu3
"""

#!python
#cython: language_level=3

from time import time
import subprocess
import os
import pickle
from collections import Counter
import fcntl
import pandas as pd
import numpy as np
from statistics import median
import argparse
import gzip





#%%
def alignFastq(fastqPath, samPath, starGenomeDir, nthreads = 48, pairedSeqPath = ''):
    '''
    use STAR to align reads in a fastq
    '''
    
    t0 = time()
    print(f'aligning {fastqPath} using {nthreads} processes...')
    if pairedSeqPath:
        bcmd = f'STAR --runThreadN {nthreads} --outSAMmultNmax 1 --outSAMunmapped Within --outSAMorder PairedKeepInputOrder --genomeDir {starGenomeDir} --readFilesIn {fastqPath} {pairedSeqPath} --outFileNamePrefix {samPath}'
    else:
        bcmd = f'STAR --runThreadN {nthreads} --outSAMmultNmax 1 --outSAMunmapped Within --outSAMorder PairedKeepInputOrder --genomeDir {starGenomeDir} --readFilesIn {fastqPath} --outFileNamePrefix {samPath}'
    subprocess.run(bcmd, shell=True)
    print(f'finished aligning {fastqPath}, {time() - t0} s \n')


#%%
def readBarcode(lane, chunk, barcodeLen, umiLen, barOutputPath, chunkOutputPath = '', isV1chem=False, tsoLen=13, ignoreSampleInd=False):
    '''
    read barcodes and merger to seqs to be aligned by STAR
    '''
    t = time()
    print(f'reading lane {lane} chunk {chunk} barcodes...')
    if isV1chem:
        print("splitting r1 for 5' reads")
    r1file = f'{lane}_R1_chunk{str(chunk).zfill(4)}'
    r2file = f'{lane}_R2_chunk{str(chunk).zfill(4)}'
    fR1 = open(r1file)
    fR2 = open(r2file)
    if not ignoreSampleInd:
        I1file = f'{lane}_I1_chunk{str(chunk).zfill(4)}'
        fI1 = open(I1file)
    try:
        counter2 = 0
        counter3 = 0
        counter = 0
        pathFreqDict = {}
        pathset = set()
        sampleFreqDict = {}
        sampleset = set()
        barFreqDict = {}
        barset = set()
        notEOF = True
        while notEOF:
            rawR1 = []
            rawI1 = []
            rawR2 = []
            iBlock = 0
            while iBlock < 20000:           #5,000 reads 
                r1Read = fR1.readline()
                r2Read = fR2.readline()
                if not r1Read:
                    notEOF = False
                    break
                rawR1.append(r1Read)
                rawR2.append(r2Read)
                if not ignoreSampleInd:
                    i1Read = fI1.readline()
                    rawI1.append(i1Read)
                counter2 += 1
                iBlock += 1
            counter3 += 1
            #in case last block before EOF is empty
            if not rawR1:
                break
            #barcode and umi from R1    
            seqR1 = rawR1[1::4]
            cellBar = [x[:barcodeLen] for x in seqR1]
            if ignoreSampleInd:
                sampleInd = ['sampleID' for x in seqR1]             #placeholder sample index
            else:
                seqI1 = rawI1[1::4]
                sampleInd = [x[:-1] for x in seqI1]
            cellPath = [f'sample_{x[0]}/cell_{x[1]}' for x in zip(sampleInd, cellBar)]
            #check barcodes
            n = len(cellPath)
            counter += n
            i = 0
            while i < n:
                if cellPath[i] in pathset:
                    pathFreqDict[cellPath[i]] += 1
                else:
                    pathset.add(cellPath[i])
                    pathFreqDict[cellPath[i]] = 1
                if sampleInd[i] in sampleset:
                    sampleFreqDict[sampleInd[i]] += 1
                else:
                    sampleset.add(sampleInd[i])
                    sampleFreqDict[sampleInd[i]] = 1
                if cellBar[i] in barset:
                    barFreqDict[cellBar[i]] += 1
                else:
                    barset.add(cellBar[i])
                    barFreqDict[cellBar[i]] = 1
                i += 1
            #v1 chemistry has 5' reads after barcodes: split out to new file 'R5'
            if isV1chem:
                rawR1[1::4] = [x[(barcodeLen + umiLen + tsoLen):] for x in seqR1]           #seq
                rawR1[3::4] = [x[(barcodeLen + umiLen + tsoLen):] for x in rawR1[3::4]]     #quality score
                writeStr = ''
                writeStr = writeStr.join(rawR1)
                with open(f'{chunkOutputPath}/{lane}_R5_chunk{chunk}', 'a') as f:
                    f.write(writeStr)
    finally:
        with open(f'{barOutputPath}/paths{lane}_{chunk}','wb') as f:
            pickle.dump(pathFreqDict, f)
        with open(f'{barOutputPath}/sampleInd{lane}_{chunk}','wb') as f:
            pickle.dump(sampleFreqDict, f)
        with open(f'{barOutputPath}/cellBarcode{lane}_{chunk}','wb') as f:
            pickle.dump(barFreqDict, f)
        fR1.close()
        if not ignoreSampleInd:
            fI1.close()
        print(f'finished processing lane {lane} chunk {chunk}, {counter} reads, {time() - t} s')    
    
    
    

#%%
def selectBarcode(laneNumber, inputDir, outputDir, barcodeLen, sampleIndexLen, cellReadsCutoff, sampleReadsCutoff):
    '''
    iterate through all cell barcodes from reads from a sequencing lane and
    returns sample_index/barcode combinations that meet the specified cutoffs 
    in a pickle
    
    inputs are pickles of dicts of cell barcode and sample index frequencies
    from readMerge(). 
    '''
    t0 = time()
    print(f'processing lane {laneNumber} barcodes...')
    #load barcode counts from file
    pathsFreqDict = {}
    pathsFreqDict = Counter(pathsFreqDict)
    pathsset = set()
    sampleFreqDict = {}
    sampleFreqDict = Counter(sampleFreqDict)
    sampleset = set()
    barFreqDict = {}
    barFreqDict = Counter(barFreqDict)
    barset = set()
    i = 0
    n = str(i).zfill(4)
    while os.path.isfile(f'{inputDir}/paths{laneNumber}_{n}'):
        with open(f'{inputDir}/paths{laneNumber}_{n}','rb') as f:
            pd = pickle.load(f)
            pathsFreqDict += pd
        with open(f'{inputDir}/sampleInd{laneNumber}_{n}','rb') as f:
            sd = pickle.load(f)
            sampleFreqDict += sd
        with open(f'{inputDir}/cellBarcode{laneNumber}_{n}','rb') as f:
            bd = pickle.load(f)
            barFreqDict += bd        
        i += 1
        n = str(i).zfill(4)
    #select barcodes based on cutoff
    nreadsSelected = 0
    nreadsTotal = 0
    for barcode, freq in pathsFreqDict.items():
        if freq > cellReadsCutoff:
            pathsset.add(barcode)
            nreadsSelected += freq
        nreadsTotal += freq
    for sampleIndex, freq in sampleFreqDict.items():
        if freq > sampleReadsCutoff:
            sampleset.add(sampleIndex)
    for barcode, freq in barFreqDict.items():
        if freq > cellReadsCutoff:
            barset.add(barcode)
    #dump results in pickles for later use
    with open(f'{outputDir}/paths_set{laneNumber}','ab') as f1:
        pickle.dump(pathsset, f1)
    with open(f'{outputDir}/sampleIndex_set{laneNumber}','ab') as f2:
        pickle.dump(sampleset, f2)
    with open(f'{outputDir}/cellBarcode_set{laneNumber}','ab') as f3:
        pickle.dump(barset, f3)
    #record summary information
    t1 = time()    
    nUniqueBarcodesComb = len(pathsFreqDict)
    nSelectedBarcodes = len(pathsset)
        #nUniqueBars = len(barFreqDict)    
    outputStr = f'found {nUniqueBarcodesComb} unique barcodes in lane {laneNumber},\n\t{nSelectedBarcodes} selected  with over {cellReadsCutoff} reads,\n\t{nreadsSelected} out of {nreadsTotal} reads to be used, {t1 - t0} s \n'
    print(outputStr)
    with open('summary.txt','a') as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        f.write(outputStr)
        fcntl.flock(f, fcntl.LOCK_UN)
    return nreadsTotal

     
#%%
def lxSplit(inputname, outputname, inputDir, outputDir, lines, skip=0):
    '''
    wrapper for linux command "split" to split large fastq for multiprocessing
    '''
    if skip != 0:
        splitcmd = f'tail {inputDir}/{inputname} -n +{skip + 1} | split -a 4 -d -l {lines} - {outputDir}/{outputname}_chunk'
    else:
        splitcmd = f'split -a 4 -d -l {lines} {inputDir}/{inputname} {outputDir}/{outputname}_chunk'
    subprocess.run(splitcmd, shell=True)
    
#%%
def lxDecompress(filename, inputDir, outputDir):
    '''
    wrapper for linux gzip to extract fastq inputs
    '''
    nFileName = filename[:-3]
    extractcmd = f'gzip -d -c {inputDir}/{filename}> {outputDir}/{nFileName}'
    subprocess.run(extractcmd, shell=True)    
    
#%%
def lxDelDir(dirpath):
    '''
    wrapper for linux delete command to recursively delete directory
    '''
    
    delcmd = f'rm -rf {dirpath}'
    subprocess.run(delcmd, shell=True)


#%%
def getFilepaths(directory):
    """
    adapted from "Johnny" (https://stackoverflow.com/users/1779256/johnny)
    
    This function will generate the file names in a directory 
    tree by walking the tree either top-down or bottom-up. 
    """
    file_paths = []
    # Walk the tree.
    for root, directories, files in os.walk(directory):
        for filename in files:
            # Join the two strings in order to form the full filepath.
            filepath = os.path.join(root, filename)
            file_paths.append(filepath)  # Add it to the list.
    return file_paths  

#%%
def sortSam(inputDir, outputDir, chunkDir, sam_name, r1chunk, i1chunk, set_name, sampleIndLen, barLen, umiLen, ignoreSampleInd, paired = False):
    '''
    sort reads in SAM files by sample index and cell barcode and store as 
    seperate SAM files
    '''
    path0 = os.getcwd()
    os.chdir(inputDir)
    fs = open(sam_name)
    fR1 = open(f'{chunkDir}/{r1chunk}')
    if not ignoreSampleInd:
        fI1 = open(f'{chunkDir}/{i1chunk}')
    try:
        print(f'sorting {sam_name}...')
        t = time()    
        #load inputs
        os.chdir(inputDir)
        if ignoreSampleInd:
            with open(set_name, 'rb') as f1:    #bar_set
                bar_set = pickle.load(f1)
        else:
            with open(set_name, 'rb') as f1:    #path_set
                path_set = pickle.load(f1)
        os.chdir(outputDir)
        #seperate reads by barcodes and write in blocks
        #   sample line: 'D000684:779:H53GNBCXY:1:1110:3316:42789_sample_GGGCAAAT_cell_TTCTACAGTGACTCAT_umi_AGGCGGTTGG\t16\tchr5\t82276115\t1\...'
        writeDict = {}
        count = 0
        n = 0
        writeCount = 0
        notEOF = True
        while notEOF:
            line = fs.readline()
            if paired:
                line2 = fs.readline()
            if line:
                count += 1
                rdump = fR1.readline()
                r1seq = fR1.readline()
                rdump = fR1.readline()
                rdump = fR1.readline()
                if not ignoreSampleInd:
                    rdump = fI1.readline()
                    i1seq = fI1.readline()
                    rdump = fI1.readline()
                    rdump = fI1.readline()
                    sampleName = i1seq[:sampleIndLen]
                else:
                    sampleName = 'sampleID'
                cellName = r1seq[:barLen]
                umi = r1seq[barLen:(barLen + umiLen)]
                line = umi + line
                if paired:
                    line2 = umi + line2
                cellPath = f'sample_{sampleName}/cell_{cellName}'
                n += 1
                if ignoreSampleInd:
                    if cellName in bar_set:
                        if cellName in writeDict:
                            writeDict[cellName].append(line)
                            if paired:
                                writeDict[cellName].append(line2)
                                writeCount += 1
                        else:
                            writeDict[cellName] = [line]
                            if paired:
                                writeDict[cellName].append(line2)
                                writeCount += 1
                        writeCount += 1
                else:
                    if cellPath in path_set:
                        if cellPath in writeDict:
                            writeDict[cellPath].append(line)
                            if paired:
                                writeDict[cellPath].append(line2)
                                writeCount += 1
                        else:
                            writeDict[cellPath] = [line]
                            if paired:
                                writeDict[cellPath].append(line2)
                                writeCount += 1
                        writeCount += 1
                if n == 500000:
                    n = 0
                    for cell in writeDict:          
                        w = ''
                        w = w.join(writeDict[cell])
                        with open(f'{cell}.sam', 'a') as f:
                            fcntl.flock(f, fcntl.LOCK_EX)
                            f.write(w)
                            fcntl.flock(f, fcntl.LOCK_UN)
                    writeDict = {}
            else:
                notEOF = False
                if writeDict:
                    for cell in writeDict:
                        w = ''
                        w = w.join(writeDict[cell])
                        with open(f'{cell}.sam', 'a') as f:
                            fcntl.flock(f, fcntl.LOCK_EX)
                            f.write(w)
                            fcntl.flock(f, fcntl.LOCK_UN)    
    finally:
        fs.close()
        fR1.close()
        if not ignoreSampleInd:
            fI1.close()
        print(f'{count} reads total from {sam_name}, {writeCount} written, {time() - t} s')
        os.chdir(path0)

#%%
def readSelection(inputSamPath, outputDir, umiLen, cellBarcodeLen, sampleIndLen, maxReadsOffset, selectOne=False, compress=True, paired=False, strict=True, maxPairOffset=100000000, checker=False):
    '''
    use alignment results to select valid reads from reads with the same UMI
    
    valids reads must be successfully mapped to the referance genome, 
    as well as be mapped to the most popular chromosome within reads that 
    share the same UMI, and have mapped locations on the chromosome within 
    maxReadsOffset of the median mapped location
    '''    
    path0 = os.getcwd()
    colNames = ['qname','flag','rname','pos','mapq','CIGAR','rnext','pnext','templ_len','seq','qual','AS','XN','XM','XO','XG','NM','MD','YS','YT']
    rml=['flag','mapq','CIGAR','rnext','pnext','templ_len','AS','XN','XM','XO','XG','NM','MD','YS','YT']
    sam = pd.read_csv(inputSamPath, sep='\t', header=None, names=colNames, index_col=False, lineterminator='\n')
    for x in rml:
        del sam[f'{x}']
    headers = sam.qname
    umi = [line[:umiLen] for line in headers] 
    headers = [('@' + x[umiLen:]) for x in headers]
    pos = sam.pos
    rname = sam.rname
    seqs = sam.seq
    quals = sam.qual
    del sam
    seqLens = [len(x) for x in seqs]
    saml = list(zip(umi, headers, rname, seqLens, pos, seqs, quals))
    nReads0 = len(saml)
    #remove reads that failed to align   
    saml = [x for x in saml if x[2] != '*']
    saml.sort()
    rname = list(zip(*saml))[2]
    umi = list(zip(*saml))[0]
    umiN = Counter(umi)
    lcell = len(saml)
    #determine which reads to select
    selectl = np.zeros(lcell, dtype = bool)
    i = 0
    while i < lcell:
        if paired:
            nreads = umiN[saml[i][0]]
            chrs = rname[i:(i + nreads)]
            chrscounter = Counter(chrs)
            chrsmode = chrscounter.most_common()[0][0]
            j = 0            
            poslistOdd = []
            TposlistOdd = []
            poslistEven = []
            TposlistEven = []
            while j < nreads:
                if (j % 2) == 1:
                    TposlistOdd.append(saml[i + j][4])
                    if chrs[j] == chrsmode:
                        poslistOdd.append(saml[i + j][4])
                    j += 1
                else:
                    TposlistEven.append(saml[i + j][4])
                    if chrs[j] == chrsmode:
                        poslistEven.append(saml[i + j][4])
                    j += 1
            npairs = len(poslistOdd)
            medposOdd = median(poslistOdd)
            medposEven = median(poslistEven)
            if selectOne:
                k = 0
                prevreadOffset = 1000000
                while k < npairs:
                    totalOffset = abs(TposlistOdd[k] - medposOdd) + abs(TposlistEven[k] - medposEven)
                    if (rname[i + (2 * k)] == chrsmode) and (totalOffset < prevreadOffset):
                        medPairInd = k
                        prevreadOffset = totalOffset
                    k += 1
                selectl[i + (2 * medPairInd)] = True
                selectl[i + (2 * medPairInd) + 1] = True
            else:
                k = 0
                while k < npairs:
                    oddOffset = abs(TposlistOdd[k] - medposOdd)
                    evenOffset = abs(TposlistEven[k] - medposEven)
                    pairOffset = abs(TposlistOdd[k] - poslistEven[k])
                    #strict qc: both reads in a pair must align to within cutoff
                    if strict:
                        if (rname[i + (2 * k)] == chrsmode) and (oddOffset < maxReadsOffset) and (evenOffset < maxReadsOffset) and (pairOffset < maxPairOffset):
                            selectl[i + (2 * k)] = True
                            selectl[i + (2 * k) + 1] = True
                    #lenient qc: only one read in a pair must align to within cutoff
                    else:
                        if (rname[i + (2 * k)] == chrsmode) and ((oddOffset < maxReadsOffset) or (evenOffset < maxReadsOffset)):
                            if pairOffset < maxPairOffset:
                                selectl[i + (2 * k)] = True
                                selectl[i + (2 * k) + 1] = True
                    k += 1                
        else:    
            nreads = umiN[saml[i][0]]
            chrs = rname[i:(i + nreads)]
            chrscounter = Counter(chrs)
            chrsmode = chrscounter.most_common()[0][0]
            j = 0
            poslist = []
            while j < nreads:
                if chrs[j] == chrsmode:
                    poslist.append(saml[i + j][4])
                j += 1
            medpos = median(poslist)
            if selectOne:
                k = 0
                prevreadOffset = 1000000
                while k < nreads:
                    offset = abs(saml[i + k][4] - medpos)
                    if (rname[i + k] == chrsmode) and (offset < prevreadOffset):
                        medInd = k
                        prevreadOffset = offset
                    k += 1
                selectl[i + medInd] = True
            else:
                k = 0
                while k < nreads:
                    if (rname[i + k] == chrsmode) and ((saml[i + k][4] - medpos) < maxReadsOffset):
                        selectl[i + k] = True
                    k += 1
        i += nreads

    nReadsSelected = 0
    if paired:
        writeBufR1 = []
        writeBufR2 = []
        i = 0
        while i < lcell:
            if selectl[i]:
                if i % 2 == 0:
                    currentRead = saml[i][1] + '\n' + saml[i][5] + '\n' + '+\n' + saml[i][6] + '\n'
                    writeBufR1.append(currentRead)
                    nReadsSelected += 1
                else:
                    currentRead = saml[i][1] + '\n' + saml[i][5] + '\n' + '+\n' + saml[i][6] + '\n'
                    writeBufR2.append(currentRead)
                    nReadsSelected += 1
            i += 1
        writeBufR1 = ''.join(writeBufR1)
        writeBufR2 = ''.join(writeBufR2)
    else:
        writeBuf = []
        i = 0
        for i in range(lcell):
            if selectl[i]:
                currentRead = saml[i][1] + '\n' + saml[i][5] + '\n' + '+\n' + saml[i][6] + '\n'
                writeBuf.append(currentRead)
                nReadsSelected += 1
        writeBuf = ''.join(writeBuf)

    os.chdir(outputDir)
    
    #testing
    if checker:
        pick = list(zip(selectl, saml))
        with open(f'{inputSamPath}_selection', 'wb') as f:
            pickle.dump(pick, f)
  
    
    cellName = inputSamPath[-(cellBarcodeLen + 4):-4]
    if compress:
        if paired:
            bwriteBufR1 = writeBufR1.encode()
            with gzip.open(f'cell_{cellName}_R1.fastq.gz', 'ab') as f:
                f.write(bwriteBufR1)
            bwriteBufR2 = writeBufR2.encode()
            with gzip.open(f'cell_{cellName}_R2.fastq.gz', 'ab') as f:
                f.write(bwriteBufR2)
        else:
            bwriteBuf = writeBuf.encode()
            with gzip.open(f'cell_{cellName}.fastq.gz', 'ab') as f:
                f.write(bwriteBuf)
    else:
        if paired:
            with open(f'cell_{cellName}_R1.fastq', 'a') as f:
                f.write(writeBufR1)
            with open(f'cell_{cellName}_R2.fastq', 'a') as f:
                f.write(writeBufR2)
        else:
            with open(f'cell_{cellName}.fastq', 'a') as f:
                f.write(writeBuf)
    os.chdir(path0)
    summary = [cellName, nReads0, lcell, nReadsSelected]
    with open('readSelection', 'ab') as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        pickle.dump(summary, f)
        fcntl.flock(f, fcntl.LOCK_UN)
    
#%%
def str2bool(v):
    '''
    allows input of boolean values from command line
    '''
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')
        
  