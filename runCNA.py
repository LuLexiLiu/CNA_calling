"""
Main program to run CNA calling.
Author: LL
Inspired by origin scripts written by Yusi Fu & Xiannian Zhang.
Version 1.0
"""

import os, json, re, zipfile
from subprocess import Popen, PIPE, check_output, call
from optparse import OptionParser
import seaborn as sns, pandas as pd, matplotlib.pyplot as plt,numpy as np

def bamMapping(bowtieRef,seqfile1,seqfile2,outputpath):
    print("--------------------------------------")
    print("Bowtie2 Mapping...")
    print("--------------------------------------")
    samfile=outputpath+'/Bowtie2.sam'
    my_env = os.environ.copy()
    bowtie_cmd = ['bowtie2', '-p', '10', '-x', bowtieRef, '-1', seqfile1, '-2', seqfile2, '-S', samfile]
    align = Popen(bowtie_cmd, env=my_env, stdout=PIPE, cwd=outputpath)
    output, err = align.communicate()

    samstats = Popen(['samtools', 'flagstat', samfile], stdout=PIPE, env=my_env)
    output, err = samstats.communicate()
    mapinfo = [[x for x in re.split('\s+', line) if len(x) > 0] for line in output.decode().split("\n")]
    mapres = [int(mapinfo[0][0]), int(mapinfo[3][0]), int(mapinfo[4][0])]

    return mapres

def binCounting(outputpath,scriptdir,dynamic_bin,binflag):
    print("--------------------------------------")
    print("Bin Counting...")
    print("--------------------------------------")
    samfile=outputpath+'/Bowtie2.sam'
    call(['perl', scriptdir + '/uniq2bin.50k.pl', samfile, dynamic_bin, outputpath + '/CNV', binflag])

    return "Bin Counting Processed."

def CBScalling(haplotype,scriptdir,countbinfile,dynamic_bin,bad_bin,sample,outpath,binflag,SD='1',alpha='0.0001'):
    print("--------------------------------------")
    print("CBS calling...")
    print("--------------------------------------")
    if haplotype == "Diploid":
        call(['Rscript', scriptdir+'/cbs.copy.R','', countbinfile, dynamic_bin, bad_bin, sample, outpath, 'notcallpeak',SD,alpha,binflag])
    elif haplotype == "Non-Diploid":
        call(['Rscript', scriptdir+'/cbs.copy.R','', countbinfile, dynamic_bin, bad_bin, sample, outpath, 'callpeak',SD,alpha,binflag])
    else:
        print("Haplotype Error. Possible haplotype values: (Diploid, Non-Diploid).")

    return "CBS Calling Processed."

def CNAplot(name, path, SD='1',alpha='0.0001',binflag='50K', ylimUp=6):
    print("--------------------------------------")
    print("CNA Profile plotting...")
    print("--------------------------------------")
    plt.style.use('default')
#     sns.set(context='talk', style='white')
    t = pd.read_csv(path + "/SD_"+SD+"_alpha_"+alpha+"_Copy." + binflag + ".txt", sep="\t")
    x = chrabspos(t)
    chrlist = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, "X", "Y")
    fig,ax=plt.subplots(figsize=(20, 2), facecolor='w')
    ax.scatter(t.abspos, t['cn.ratio'], s=8, c='slategrey', alpha=1)
    ax.plot(t.abspos, t['cn.seg'], color='#ff7f0e', alpha=1, lw=3.5)
    for i in x:
        ax.axvline(x=i, c='black', lw=0.35, ls='dashed', alpha=0.5)

    ax.set_xticks(c_loc(chrlist, t))
    ax.set_xticklabels(chrlist, color='k', size=19)
    ax.set_yticks((0, 2, 4, 6))
    ax.set_yticklabels((0, 2, 4, 6), color='k',size=20)
    ax.set_xlim(0, list(t.abspos)[-1])
    ax.set_ylim(0, ylimUp)
    ax.set_ylabel("Copy Number",color='k',size=20)
    ax.set_title(str(name),color='k',size=20)
    ax.set_xlabel("Chromosome",color='k',size=20)
    fig.savefig(path + '/CNAplot.pdf', bbox_inches = 'tight')

def chrabspos(t):
    chrpos = list(t.groupby('chrom').apply(lambda x: x['abspos'].max()))
    return chrpos

def c_loc(chrlist, t):
    binpos = chrabspos(t)
    chrloc = []
    for i in chrlist:
        try:
            if i == 1:
                chrloc.append(int(binpos[i - 1] / 2))
            else:
                chrloc.append(int((binpos[i - 2] + binpos[i - 1]) / 2))
        except:
            if i == "X":
                chrloc.append(int((binpos[21] + binpos[22]) / 2))
            if i == "Y":
                chrloc.append(int((binpos[22] + binpos[23]) / 2))
    return chrloc

def findDFsegments(DF,columnname,chrom=True):
    DF['block']=(DF[columnname].shift(1) != DF[columnname]).astype(int).cumsum()
    if chrom:
        matrix=DF.reset_index().groupby([columnname,'block','chrom'])['index'].apply(list).values
    else:
        matrix=DF.reset_index().groupby([columnname,'block'])['index'].apply(list).values
    return matrix

def GetCNAposition(path,SegFile,sex,CNcriterion=0.5):
    DF=pd.read_csv(path + "/"+SegFile, sep="\t")
    Segdic={'chrom':[],
            'begin':[],
            'end':[],
            'cntype':[]}
    normw=[2]*len(DF[DF.chrom!=24])+[0]*len(DF[DF.chrom==24])
    normm=[2]*len(DF[(DF.chrom!=24)&(DF.chrom!=23)])+[1]*len(DF[(DF.chrom==24)|(DF.chrom==23)])
    if sex.lower()=='m':
        DF['DiffCN']=DF['cn.seg']-normm
    elif sex.lower()=='f':
        DF['DiffCN']=DF['cn.seg']-normw
    cnadf=DF[[abs(i)>=CNcriterion for i in DF['DiffCN']]]
    matrix=findDFsegments(cnadf,'cn.seg',chrom=True)
    for segIndex in matrix:
#         slicedf=DF[DF.index.isin(segIndex)]
        Segdic['chrom'].append(DF.chrom[segIndex[0]])
        Segdic['begin'].append(DF.chrompos[segIndex[0]])
        if DF.chrom[segIndex[-1]+1]:
            if DF.chrom[segIndex[-1]]==DF.chrom[segIndex[-1]+1]:
                Segdic['end'].append(DF.chrompos[segIndex[-1]+1])
            else:
                Segdic['end'].append(DF.chrompos[segIndex[-1]])
        else:
            Segdic['end'].append(DF.chrompos[segIndex[-1]])
        Segdic['cntype'].append("<CN_"+str(round(DF['cn.seg'][segIndex[0]],1))+">")
    return pd.DataFrame(Segdic).sort_values(by='chrom')

def annotSVcalling(path,annotSVpath,genomebuild):
    call([
        annotSVpath , '-SVinputFile', path+'/CNAposition.bed', '-SVinputInfo', '1',
        '-outputFile',path+'/CNA.annotated.tsv',
        '-svtBEDcol', '4', '-genomeBuild', genomebuild
    ]
    )

    return "AnnotSV processed."

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-c", "--config", dest="config",
                      help="configure files", metavar="config")
    parser.add_option("-i", "--sampleid", dest="sampleid",
                      help="running commands for samples", metavar="sampleid")
    parser.add_option("-m","--mode",dest="mode",
                      help="Choose the program mode. Possible values: (CNAcalling, Annotation). Default: CNAcalling.",
                      metavar="mode")
    parser.add_option("-s", "--segfile", dest="segfile",
                      help="Name of Segment File.",
                      metavar="segfile")

    (options, args) = parser.parse_args()

    with open(options.config, 'r') as infile:
        cfg = json.load(infile)
    if options.sampleid:
        outdir = cfg['OutputPath']
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        sample = options.sampleid
        path = os.path.join(outdir, sample)
        if not os.path.exists(path):
            os.mkdir(path)
        if options.mode=='Annotation':
            if options.segfile:
                SegFile=options.segfile
            else:
                SegFile="SD_"+cfg["SD"]+"_alpha_"+cfg["alpha"]+"_Copy." + cfg["BinFlag"] + ".txt"
            GetCNAposition(path, SegFile, cfg['Samples'][sample]['sex'],
                           float(cfg['CN_Diff_Criterion'])).to_csv(path+'/CNAposition.bed',sep='\t',header=None,index=None)
            annotSVcalling(path, cfg['AnnotSVpath'], cfg['GenomeBuild'])
        else:
            mapping=bamMapping(cfg['BowtieRef'], cfg['Samples'][sample]['seq1'], cfg['Samples'][sample]['seq2'], path)
            logfile=os.path.join(path+'/mapping_log.json')
            out = {
               'sampleID': sample,
               'total reads':mapping[0],
               'mappred reads':mapping[2]
            }
            with open(logfile, 'a') as logjson:
               json.dump(out,logjson)
            binCounting(path, cfg['ScriptDir'], cfg['Dynamic_Bin'], cfg['BinFlag'])
            CBScalling(cfg['Samples'][sample]['ploidy'],cfg['ScriptDir'],
                       os.path.join(path,'CNV.dyna_'+cfg['BinFlag']+'.bin'),
                       cfg['Dynamic_Bin'], cfg['Bad_Bin'],
                       sample, path, cfg['BinFlag'],
                       cfg["SD"],cfg["alpha"]
                       )
            CNAplot(sample, path, cfg["SD"],cfg["alpha"], cfg['BinFlag'])



