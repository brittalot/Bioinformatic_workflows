import sys

# open files
fq1 = open(sys.argv[1])
fq2 = open(sys.argv[2])
fq1out = open(sys.argv[3],'w')
fq2out = open(sys.argv[4],'w')
singleReads = open(sys.argv[5],'w')

readcount = 0
pairscount = 0
singlereadscount = 0

while True:
    try:

        # get read1 from fq1
        header1 = fq1.readline().rstrip()
        seq1    = fq1.readline().rstrip()
        junk1   = fq1.readline().rstrip()
        qual1   = fq1.readline().rstrip()
        
        # get read2 from fq2
        header2 = fq2.readline().rstrip()
        seq2    = fq2.readline().rstrip()
        junk2   = fq2.readline().rstrip()
        qual2   = fq2.readline().rstrip()

        if header1 == '': sys.stderr.write('Header one is empty exiting.\n');break
        
        readcount += 1
        
        # check if read-sequence is zero
        if len(seq1) and len(seq2):
            # print to outfiles
            fq1out.write(header1+'\n'+seq1+'\n'+junk1+'\n'+qual1+'\n')
            fq2out.write(header2+'\n'+seq2+'\n'+junk2+'\n'+qual2+'\n')
            pairscount += 1
            
        else:
            
            if len(seq1):
                singleReads.write(header1+'\n'+seq1+'\n'+junk1+'\n'+qual1+'\n')
                singlereadscount += 1
                
            elif len(seq2):
                singleReads.write(header2+'\n'+seq2+'\n'+junk2+'\n'+qual2+'\n')
                singlereadscount += 1
                
            elif not len(seq1) and not len(seq2):
                sys.stderr.write('Skipping read '+header1+' as sequence length is zero.\n')

    except EOFError:
        break

sys.stderr.write(str(readcount)+' reads processed.\n')
sys.stderr.write(str(pairscount)+' read pairs to outfiles '+fq1out.name+' and '+ fq2out.name+'.\n')
sys.stderr.write(str(singlereadscount)+' single reads to outfile '+singleReads.name+'.\n')

fq1.close()
fq2.close()
fq1out.close()
fq2out.close()
singleReads.close()
