# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 11:54:09 2019

@author: g06t
"""

import al_functions_m93 as af
from time import time
import os
from multiprocessing import Pool
import pickle
import psutil
import argparse
import sys
#import fnmatch
from glob import glob
import csv


#%%
def main():    
    t0 = time()   
    #parse user inputs
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, 
                             epilog='example: \n \t python FASTR.py --f input_names.txt --i /home/seqs/example --r /home/results/example_run --ig True --u False')
    parser.add_argument('--f', type=str, default='inputNames.txt',
                        help='file of filenames')
    parser.add_argument('--i', type=str, default='',
                        help='path of the directory of the input fastq files')
    parser.add_argument('--r', type=str, default='',
                        help='path of the directory for output')
    parser.add_argument('--ind', type=str, default='',
                        help='path of the STAR reference index')
    parser.add_argument('--chver', type=int, default=2,
                        help='10x chemistry version; options: 1 (v1 chemistry paired), 2 (v2 chemistry); default: 2')
    parser.add_argument('--qc', type=str, default='lenient',
                        help='read QC mode for v1 chemistry paired reads; options:"lenient," "strict"; default: "lenient"')
    parser.add_argument('--cb', type=int, default=16,
                        help='length of cell barcode; default = 16 bp')
    parser.add_argument('--sl', type=int, default=8,
                        help='length of sample index; default = 8 bp')
    parser.add_argument('--ul', type=int, default=10,
                        help='length of UMI; default = 10 bp')
    parser.add_argument('--tsol', type=int, default=13,
                        help='length of TSO, for older 10x Genomics v1 chemistry paired end sequencing; default = 13 bp')
    parser.add_argument('--cc', type=int, default=3000,
                        help='cell reads cutoff: minimum number of reads required for a cell barcode to be valid; default = 3000')
    parser.add_argument('--sc', type=int, default=1000000,
                        help='sample reads cutoff: minimum number of reads required for a sample index to be valid; default = 1,000,000')
    parser.add_argument('--ro', type=int, default=500,
                        help='read offset: maximum deviation from the median alignment position of a group of reads with identical UMIs for a read to be valid; default = 500 bp')
    parser.add_argument('--po', type=int, default=50000000,
                        help='paired read offset: maximum difference in alignment position between two mates for a paired set of reads to be valid; default = 50,000,000 bp')
    parser.add_argument('--so', type=af.str2bool, nargs='?', const=True, default=False,
                        help='select one: choose whether to select one read per UMI (True) or keep all valid reads (False)')
    parser.add_argument('--ig', type=af.str2bool, nargs='?', const=True, default=False,
                        help='force ignore sample index: if set to "True," sort reads as if they have identical sample indices; default = False, automatically set to True if I1 filenames are not provided')
    parser.add_argument('--u', type=af.str2bool, nargs='?', const=True, default=False,
                        help='pause program for user to confirm parameters before starting; default = False')
    parser.add_argument('--gz', type=af.str2bool, nargs='?', const=True, default=True,
                        help='compressed output; default = True (use ''False'' to keep results as uncompressed fastq')
    parser.add_argument('--sam', type=af.str2bool, nargs='?', const=True, default=False,
                        help=f"(debugging) keep SAM: if set to ""True"", keeps SAM files from alignment for checking; default: False")
    args = parser.parse_args()
    #transfer inputs to variables
    nameFile = args.f
    inputDir0 = args.i
    outputDir0 = args.r
    refind = args.ind
    chemVer = args.chver
    if chemVer != 2 and chemVer != 1:
        print(f'Error -1: unsupported 10x chemistry version {chemVer}')
        return -1
    barcodeLen = args.cb 
    sampleIndLen = args.sl  
    umiLen = args.ul 
    tsoLen = args.tsol
    cellReadsCutoff = args.cc 
    sampleReadsCutoff = args.sc 
    maxReadOffset = args.ro 
    selectOne = args.so 
    ignoreSampleInd = args.ig 
    waitForConfirm = args.u
    gzresults = args.gz
    keepSAM = args.sam
    maxPairCutoff = args.po
    qcMode = args.qc
    if qcMode == 'strict':
        isStrict = True
    else:
        isStrict = False
    print('\n\n\n\n\n\n')
    #set input directory
    p0 = os.getcwd()
    if not inputDir0:
        inputDir0 = p0
        inputstr = 'input directory: current working directory'
        print(inputstr)
    else:
        inputstr = f'input directory: {inputDir0}'
        print(inputstr)
    #import filenames of fastq inputs
    inputNames = []
    with open(nameFile, newline='') as f:
        nameReader = csv.reader(f, delimiter='\t')
        for row in nameReader:
            inputNames.append(row)
    #check if inputs are compressed    
    if '.gz' in inputNames[0][0]:
        compressedInputs = True
    else:
        compressedInputs = False
    lanes = len(inputNames)
    n0 = len(inputNames[0])
    if n0 == 2:
        if not ignoreSampleInd:
            print('caution: one of the three read files (R1, R2, or I1) was not specified in inputNames\n\
                  \t ***the run will proceed assuming I1 is not provided and all sample indices are presumed to be identical***')
        ignoreSampleInd = True
    if not ignoreSampleInd:
        if isStrict:
            print('warning: qc mode is set to ''strict'' but is not relavent for unpaired reads; strict mode is reset to False\n')
            isStrict = False
    if (n0 > 3) or (n0 < 2):
        print(f'\n\nERROR -2:invalid number of input files: {n0}\n\n')
        return -2
    i = 0
    n1 = 0
    while i < lanes:
        n1 = len(inputNames[i])
        if n1 != n0:
            print(f'\n\nERROR -3: filename format error in row {i + 1} of {nameFile}: {n1} file names found, expected {n0}\n\n')
            return -1
        n0 = n1
        i += 1            
    if not outputDir0:
        outputDir0 = f'{inputDir0}/split_reads'
        os.mkdir(outputDir0)
        outputstr = 'output directory: new directory in input directory'
        print(outputstr)
    else:
        outputstr = f'output directory: {outputDir0}'
        print(outputstr)        
    #make output directory if it does not exist
    if not os.path.isdir(outputDir0):
        os.mkdir(outputDir0)
    os.chdir(outputDir0)
    #record input/output directory locations
    with open('summary.txt','a') as f:
        wrtstr = inputstr + '\n' + outputstr + '\n'
        f.write(wrtstr)
    #set bowtie2 reference index 
    '''
    if not refind:
        btind = f'{p0}/bowtie2_index/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index'
    else:
        btind = refind
    '''
    if not glob(refind + '/SA'):
        print(f'Error -2: no genome index found at {refind}!')
        return -2
    #display/record user inputs 
    if chemVer == 1:   
        versionStr = '10x v1 chemistry paired'
    else:
        versionStr = '10x v2 chemstry'
    l = {'fastqs':f'{inputNames[0][0]}, {inputNames[0][1]}, ...', 'STAR index': refind, 'barcode_length':barcodeLen,'sampl.ind.len': sampleIndLen,\
         'UMI len': umiLen,'cellReadsCutoff': cellReadsCutoff,'sampleReadsCutoff': sampleReadsCutoff,\
         'maxReadOffset': maxReadOffset,'select one read per UMI': selectOne,\
         'ignore sample indices': ignoreSampleInd, 'number of lanes': lanes,\
         'user confirmation of inputs': waitForConfirm, 'compressed inputs': compressedInputs, 'keep SAMs': keepSAM,\
         'version': versionStr, 'QC mode': qcMode, 'max pair offset': maxPairCutoff}
    for x, y in l.items():
        sumstr0 = f'{x} :\t {y} \n'
        print(sumstr0)
        with open('summary.txt','a') as f:
            f.write(sumstr0)
    if waitForConfirm:
        p = input('proceed? (y to continue) \n')
        if p != 'y' and p != 'Y':
            sys.exit(0)            
    #calculate chunk size based on system resources
    nprocess = os.cpu_count()
    availMem = psutil.virtual_memory()[1] / 1073741824
    '''
    #unused (7/26/19): STAR uses so much memory that memory use is independent of chunk size
    processMem = availMem / nprocess
    maxnReads = processMem / 900        #max mem consumption is ~0.7 KB per read per process
    if maxnReads > 1000000:
        chunkSize = 1000000             #samSort slows down dramatically with excessive chunk size
    else:
        chunkSize = int(maxnReads)
    '''
    chunkSize = 5000000    
    #record system resource information
    iniStr = f'{availMem} GB memory available, running with {nprocess} processes \n\n'
    '''
    STAR memory usage (from author):
        "The required RAM is ~10*GenomeSize bytes for the genome indices, 
        plus per-thread buffer memory (150MB per thread, 
        but this can be reduced with --limitIObufferSize option)."
    '''
    #check available mamory (using hg38 ref index ~31GB needed for index)
    startMemRequirement = 32 + 150/1024*nprocess
    if availMem < startMemRequirement:
        warnStr= '***Warning: possible insuficient free memory detected! Run may hang during alignment***'
        iniStr += warnStr
    print(iniStr)
    with open('summary.txt','a') as f:
        f.write(iniStr)
    #make folders
    path0 = outputDir0
    if not os.path.isdir('chunks'):
        os.mkdir('chunks')
    os.chdir('chunks')
    pathC = os.getcwd()
    os.chdir(path0)
    if not os.path.isdir('barcode_pickles'):
        os.mkdir('barcode_pickles')
    os.chdir('barcode_pickles')
    pathB = os.getcwd()
    os.chdir(path0)
    if not os.path.isdir('sams'):
        os.mkdir('sams')
    os.chdir('sams')
    pathS = os.getcwd()
    os.chdir(path0)
    if not os.path.isdir('sam_chunks'):
        os.mkdir('sam_chunks')
    os.chdir('sam_chunks')
    pathSC = os.getcwd()
    os.chdir(path0)
    if not os.path.isdir('sorted_sams'):
        os.mkdir('sorted_sams')
    os.chdir('sorted_sams')
    pathS2 = os.getcwd()
    os.chdir(path0)
    if not os.path.isdir('results'):
        os.mkdir('results')
    os.chdir('results')
    pathR = os.getcwd()
    os.chdir(path0)
    #if inputs are compressed fastqs, decompress to a new directory and set that as input directory
    if compressedInputs:
        t1 = time()
        print('decompressing inputs...')
        if not os.path.isdir('extracted_inputs'):
            os.mkdir('extracted_inputs')
        os.chdir('extracted_inputs')
        pathEx = os.getcwd()
        os.chdir(path0)    
        inputlist = []
        for lane in inputNames:
            for name in lane:
                inputlist.append((name,inputDir0, pathEx))
        with Pool() as pool:
            pool.starmap(af.lxDecompress, inputlist)    
        inputNames = [[filename[:-3] for filename in lane] for lane in inputNames]
        #set directory of extracted files as input directory
        inputDir0 = pathEx
        #log time
        t2 = time()
        summstrEx = f'finished decompressing inputs, {t2 - t1} s, {t2 - t0} s total\n\n'
        t1 = t2
        print(summstrEx)
        with open('summary.txt','a') as f:
            f.write(summstrEx)  
    #%%split inputs into chunks for multiprocessing
    t1 = time()
    print('making chunks...')
    inputlist = []
    lines = chunkSize * 4
    inputlist = []
    laneNum = 1
    if ignoreSampleInd:
        readType = ['R1','R2']
    else:
        readType = ['R1','R2','I1']
    readTypeInd = 0
    for lane in inputNames:
        if ignoreSampleInd:
            i = 0
            while i < 2:
                inputlist.append((lane[i],f'{laneNum}_{readType[i]}', inputDir0, pathC, lines))
                i += 1
        else:
            for name in lane:
                inputlist.append((name,f'{laneNum}_{readType[readTypeInd]}', inputDir0, pathC, lines))
                readTypeInd += 1
        laneNum += 1
        readTypeInd = 0
    with Pool() as pool:
        pool.starmap(af.lxSplit, inputlist)    
    #log time
    t2 = time()
    print(f'finished making chunks, {t2 - t1} s, {t2 - t0} s total')
    t1 = t2
    #delete extracted inputs to free disk space
    if compressedInputs:
        print('deleting extracted inputs...')
        af.lxDelDir(pathEx)
        print('finished deleting\n')
    os.chdir(pathC)
    #find number of chunks made
    nchunks = []
    laneNum = 0
    while laneNum < lanes:
        n = 0
        while os.path.isfile(f'{laneNum + 1}_{readType[0]}_chunk{str(n).zfill(4)}'):
            n += 1
        nchunks.append(n)
        laneNum += 1
    os.chdir(path0)    
    #%%process barcodes
    print('processing barcodes...')
    os.chdir(pathC)
    inputlist = []
    laneNum = 0
    if chemVer == 1:
        isV1chem = True
    else:
        isV1chem = False
    while laneNum < lanes:
        chunkNum = 0
        while chunkNum < nchunks[laneNum]:
            inputlist.append(((laneNum + 1), str(chunkNum).zfill(4), barcodeLen, umiLen, pathB, pathC, isV1chem, tsoLen, ignoreSampleInd))
            chunkNum += 1
        laneNum += 1
    with Pool() as pool:
        pool.starmap(af.readBarcode, inputlist)
    os.chdir(path0)
    #select barcodes
    t1 = time()
    print('selecting barcodes...')
    inputlist = []
    laneNum = 0
    while laneNum < lanes:
        inputlist.append(((laneNum + 1), pathB, pathSC, barcodeLen, sampleIndLen, cellReadsCutoff, sampleReadsCutoff))
        laneNum += 1
    with Pool() as pool:
        totalReadsList = pool.starmap(af.selectBarcode, inputlist)
    totalReads = sum(totalReadsList)
    if chemVer == 1:
        totalReads *= 2
    t2 = time()
    print(f'finished barcode selection, {t2 - t1} s, {t2 - t0} s total')
    t1 = t2
    #%%alignment
    os.chdir(pathC)
    print('starting reads alignment...')
    i = 0 
    while i < lanes:
        r2chunkNames = []
        r5chunkNames = []
        n = 0
        while n < nchunks[i]:
            r2chunkNames.append(f'{i + 1}_R2_chunk{str(n).zfill(4)}')
            r5chunkNames.append(f'{i + 1}_R5_chunk{str(n).zfill(4)}')
            n += 1
        read1str = ','.join(r5chunkNames)
        read2str = ','.join(r2chunkNames)
        if isV1chem:
            af.alignFastq(read1str, f'{pathS}/lane_{i + 1}', refind, nprocess, read2str)
        else:
            af.alignFastq(read2str, f'{pathS}/lane_{i + 1}', refind, nprocess)
        i += 1
    t2 = time()
    print(f'finished read alignment, {t2 - t1} s, {t2 - t0} s total')
    t1 = t2
    #%%split alignment result for multiprocessing
    print('splitting sam...')
    inputlist = []
    laneNum = 0
    lines = int(lines / 4)
    if isV1chem:
        #pair end alignments produce two lines per pair
        lines = int(lines * 2)
    os.chdir(pathS)
    while laneNum < lanes:
        samName = f'lane_{laneNum + 1}Aligned.out.sam'
        #sam header lines need to be skipped to match order with previous chunks
        with open(samName) as fsm:
            l = fsm.readline()
            skipLines = 0
            while l[0] == '@':
                l = fsm.readline()
                skipLines += 1
        inputlist.append((samName, f'lane_{laneNum + 1}_sam', pathS, pathSC, lines, skipLines))
        laneNum += 1
    os.chdir(path0)
    #split SAM
    with Pool() as pool:
        pool.starmap(af.lxSplit, inputlist)    
    #log time
    t2 = time()
    print(f'finished making SAM chunks, {t2 - t1} s, {t2 - t0} s total')
    t1 = t2
    #%%sorting sam files based on barcodes
    print(f'started sorting sams...')
    os.chdir(pathSC)
    inputlist = []
    sampleIndexSet = set()
    laneNum = 0
    while laneNum < lanes:
        chunkNum = 0
        while chunkNum < nchunks[laneNum]:
            samChunkName = f'lane_{laneNum + 1}_sam_chunk{str(chunkNum).zfill(4)}'
            r1chunkName = f'{laneNum + 1}_R1_chunk{str(chunkNum).zfill(4)}'
            i1chunkName = f'{laneNum + 1}_I1_chunk{str(chunkNum).zfill(4)}'
            pathSet = f'paths_set{laneNum + 1}'
            sampleIndSetName = f'sampleIndex_set{laneNum + 1}'
            barSet = f'cellBarcode_set{laneNum + 1}'
            if ignoreSampleInd:
                inputlist.append((pathSC, pathS2, pathC, samChunkName, r1chunkName, i1chunkName, barSet, sampleIndLen, barcodeLen, umiLen, ignoreSampleInd, isV1chem))
            else:
                inputlist.append((pathSC, pathS2, pathC, samChunkName, r1chunkName, i1chunkName, pathSet, sampleIndLen, barcodeLen, umiLen, ignoreSampleInd, isV1chem))
            with open(sampleIndSetName, 'rb') as f:
                sampleIndexSet.update(pickle.load(f))    
            chunkNum += 1
        laneNum += 1
    if not ignoreSampleInd:
        os.chdir(pathS2)    
        for sample in sampleIndexSet:
            if not os.path.isdir(f'sample_{sample}'):
                os.mkdir(f'sample_{sample}')
    os.chdir(path0)
    with Pool() as pool:    
        pool.starmap(af.sortSam, inputlist)
    t2 = time()
    print(f'finished sorting sams, {t2 - t1} s, {t2 - t0} s total')
    t1 = t2
    print('deleting intermediate files...')
    af.lxDelDir(pathC)
    af.lxDelDir(pathB)
    af.lxDelDir(pathS)
    af.lxDelDir(pathSC)
    #%%selecting reads
    print('selecting reads based on alignment...')
    if ignoreSampleInd:
        cellSamPaths = af.getFilepaths(pathS2)
        inputlist = []
        for path in cellSamPaths:
            inputlist.append((path, pathR, umiLen, barcodeLen, sampleIndLen, maxReadOffset, selectOne, gzresults, isV1chem, isStrict, maxPairCutoff, keepSAM))
        with Pool() as pool:
            pool.starmap(af.readSelection, inputlist)
    else:
        for sample in sampleIndexSet:
            sampleDir = f'sample_{sample}'
            outputDir = f'{pathR}/{sampleDir}'
            if not os.path.isdir(outputDir):
                os.mkdir(outputDir)
            inputDir = f'{pathS2}/{sampleDir}'
            cellSamPaths = af.getFilepaths(inputDir)
            inputlist = []
            for path in cellSamPaths:
                inputlist.append((path, outputDir, umiLen, barcodeLen, sampleIndLen, maxReadOffset, selectOne, gzresults, isV1chem, isStrict, maxPairCutoff, keepSAM))
            with Pool() as pool:
                pool.starmap(af.readSelection, inputlist)
    t2 = time()
    print(f'finished read selection, {t2 - t1} s, {t2 - t0}s total')
    if not keepSAM:
        print('deleting intermediate files...')
        af.lxDelDir(pathS2)
    #%%calculate and record summary
    f = open('readSelection', 'rb')
    summaryList = []
    while True:
        try:
            summaryList.append(pickle.load(f))
        except EOFError:
            break 
    f.close()
    nReads_correct_barcode = 0
    nReads_aligned = 0
    nReads_selected = 0
    for cell in summaryList:
        nReads_correct_barcode += cell[1]
        nReads_aligned += cell[2]
        nReads_selected += cell[3]
    percentCBar = format((nReads_correct_barcode / totalReads * 100), '.2f')
    percentAln = format((nReads_aligned / totalReads * 100), '.2f')
    percentSort = format((nReads_selected / nReads_aligned * 100), '.2f')
    percentOverall = format((nReads_selected / totalReads * 100), '.2f')
    t2 = time()
    summStr4 = f'of {totalReads} total reads, {nReads_correct_barcode} found with valid barcodes ({percentCBar} %), {nReads_aligned} aligned to genome ({percentAln} %),\n\t {nReads_selected} selected and sorted ({percentSort} % of aligned reads, {percentOverall} % of total)\n\nrun comepleted, {t2 - t0} s total'
    with open('summary.txt','a') as f:
        f.write(summStr4) 
    print(summStr4)


    
if __name__ == '__main__':
    main()
