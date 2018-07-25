#!/usr/bin/env python2
import argparse
import subprocess

def get_picard_dup_stats(picard_dup_file, paired_status):
    '''
    Parse Picard's MarkDuplicates metrics file
    '''
    mark = 0
    dup_stats = {}
    with open(picard_dup_file) as fp:
        for line in fp:
            if '##' in line:
                if 'METRICS CLASS' in line:
                    mark = 1
                continue

            if mark == 2:
                line_elems = line.strip().split('\t')
                ndup_idx,pdup_idx= (5,7) if line_elems[5]==0 else (6,8)
                dup_stats['PERCENT_DUPLICATION'] = line_elems[pdup_idx]
                dup_stats['READ_PAIR_DUPLICATES'] = line_elems[ndup_idx]
                dup_stats['READ_PAIRS_EXAMINED'] = line_elems[2]
                if paired_status == 'Paired-ended':
                    return 2*int(line_elems[ndup_idx]), float(line_elems[pdup_idx])
                else:
                    return int(line_elems[4]), float(line_elems[7])

            if mark > 0:
                mark += 1
    return None

def parse_args():
    ## parse arg 
    parser = argparse.ArgumentParser(description='correct dup number form dup.qc')
    parser.add_argument('--lib', help='lib id')
    
    args = parser.parse_args()
    return args.lib


def main():
    
    '''
    correct "Read count after removing duplicate reads" in lib_qc file
    by re-parsing dup_qc file
    '''
    
    lib = parse_args()
    libqc_dir='/home/zhc268/data/outputs/libQCs/'
    cmd="find "+libqc_dir+lib+"/"+" -name '*_qc.txt'"
    libqc_file=subprocess.check_output(cmd,shell=True).strip('\n')
    
    cmd="find "+libqc_dir+lib+"/"+" -name '*.dup.qc'"    
    dupqc_file=subprocess.check_output(cmd,shell=True).strip('\n')

    with open(libqc_file,'r') as f: 
        libqc= f.readlines()

    ## get the two lines
    ss =['Read count after filtering for mapping quality\t',
         'Read count after removing duplicate reads\t',
         'Duplicates (after filtering)']
    
    qc_ss= [[[i,int(e.strip(s))] for i,e in enumerate(libqc) if s in e] for s in ss[0:2]]
    
    
    ## correct & write file
    if qc_ss[0][0][1] == qc_ss[1][0][1]:

        ## update ndup 
        paired_status='Paired-ended' if "_R1" in libqc_file else 'Single-end'
        ndup,pdup=get_picard_dup_stats(dupqc_file,paired_status)
        nv=str(qc_ss[1][0][1] - ndup)
        libqc[qc_ss[1][0][0]] = libqc[qc_ss[1][0][0]].replace(str(qc_ss[0][0][1]),nv)

        # update ndup fraction 
        idx= [i  for i,e in enumerate(libqc) if ss[2] in e]
        libqc[idx[0]]='Duplicates (after filtering)\t'+str(ndup)+'\t'+str(pdup)+'\n'
        with open(libqc_file,'w') as f: f.writelines(libqc)    
    return None
     


if __name__ == '__main__':
    main()


