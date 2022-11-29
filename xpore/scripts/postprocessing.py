import os

def readAnnotation(gtf_or_gff):
    gtf=open(gtf_or_gff,"r")
    dict,is_gff={},0
    for ln in gtf:
        if not ln.startswith("#"):
            ln=ln.strip("\n").split("\t")
            if is_gff == 0:
                if ln[-1].startswith("ID") or ln[-1].startswith("Parent"):
                    is_gff = 1
                else:
                    is_gff = -1
            if is_gff < 0:
                if ln[2] == "transcript" or ln[2] == "exon":
                    chr,type,start,end=ln[0],ln[2],int(ln[3]),int(ln[4])
                    attrList=ln[-1].split(";")
                    attrDict={}
                    for k in attrList:
                        p=k.strip().split(" ")
                        if len(p) == 2:
                            attrDict[p[0]]=p[1].strip('\"')
                    ##tx_id=ln[-1].split('; transcript_id "')[1].split('";')[0]
                    ##g_id=ln[-1].split('gene_id "')[1].split('";')[0]
                    tx_id = attrDict["transcript_id"]
                    g_id = attrDict["gene_id"]
                    g_name = attrDict["gene_name"]
                    if tx_id not in dict:
                        dict[tx_id]={'chr':chr,'g_id':g_id,'g_name':g_name,'strand':ln[6]}
                        if type not in dict[tx_id]:
                            if type == "transcript":
                                dict[tx_id][type]=(start,end)
                    else:
                        if type == "exon":
                            if type not in dict[tx_id]:
                                dict[tx_id][type]=[(start,end)]
                            else:
                                dict[tx_id][type].append((start,end))
            if is_gff > 0:
                if ln[2] == "exon" or ln[2] == "mRNA":
                    chr,type,start,end=ln[0],ln[2],int(ln[3]),int(ln[4])
                    tx_id=ln[-1].split('transcript:')[1].split(';')[0]
                    if ln[2] == "mRNA":
                        type="transcript"
                    if tx_id not in dict:
                        dict[tx_id]={'chr':chr,'strand':ln[6]}
                        if type == "transcript":
                            dict[tx_id][type]=(start,end)
                        if type == "exon":
                            dict[tx_id][type]=[(start,end)]
                    else:
                        if type == "transcript" and type not in dict[tx_id]:
                            dict[tx_id][type]=(start,end)
                        if type == "exon":
                            if type not in dict[tx_id]:
                                dict[tx_id][type]=[(start,end)]
                            else:
                                dict[tx_id][type].append((start,end))
    #convert genomic positions to tx positions
    if is_gff < 0:
        for id in dict:
            tx_pos,tx_start=[],0
            for pair in dict[id]["exon"]:
                tx_end=pair[1]-pair[0]+tx_start
                tx_pos.append((tx_start,tx_end))
                tx_start=tx_end+1
            dict[id]['tx_exon']=tx_pos
    else:
        for id in dict:
            tx_pos,tx_start=[],0
            if dict[id]["strand"] == "-":
                dict[id]["exon"].sort(key=lambda tup: tup[0], reverse=True)
            for pair in dict[id]["exon"]:
                tx_end=pair[1]-pair[0]+tx_start
                tx_pos.append((tx_start,tx_end))
                tx_start=tx_end+1
            dict[id]['tx_exon']=tx_pos
    return dict

def addChrGeneNameToDiffmodTab(annotation_dict,diffmod_table_path,genome,out_dir):
    if genome:
        original_annotation_dict=annotation_dict
        annotation_dict={}
        for tx_id in original_annotation_dict:
            annotation_dict[original_annotation_dict[tx_id]['g_id']]=original_annotation_dict[tx_id]
    file=open(diffmod_table_path,"r")
    header=file.readline()
    entries=file.readlines()
    outfile_path=os.path.join(out_dir,"wChrGName_diffmod.table")
    outfile=open(outfile_path,"w")
    outfile.write('chr,gene_name,'+header)
    header=header.strip().split(',')
    id_ind=header.index('id')
    for ln in entries:
        id=ln.split(',')[id_ind]
        outfile.write(annotation_dict[id]['chr']+','+annotation_dict[id]['g_name']+','+ln)
    outfile.close()
    return outfile_path

def run_postprocessing(diffmod_table_path,out_dir):
    file=open(diffmod_table_path,"r")
    header=file.readline()
    entries=file.readlines()
    outfile_path=os.path.join(out_dir,"majority_direction_kmer_diffmod.table")
    outfile=open(outfile_path,"w")
    outfile.write(header)
    header=header.strip().split(',')
    kmer_ind,dir_ind=header.index('kmer'),header.index('mod_assignment')
    dict={}
    for ln in entries:
        l=ln.strip().split(",")
        if l[kmer_ind] not in dict:
            dict[l[kmer_ind]]={l[dir_ind]:1}
        else:
            if l[dir_ind] not in dict[l[kmer_ind]]:
                dict[l[kmer_ind]][l[dir_ind]]=1
            else:
                dict[l[kmer_ind]][l[dir_ind]]+=1
    for k in dict:
        if len(dict[k]) > 1:  ##consider one modification type per k-mer
            if dict[k]['higher'] <= dict[k]['lower']: ##choose the majority
                dict[k]['choose']='lower'
            else:
                dict[k]['choose']='higher'
        else:
            dict[k]['choose']=list(dict[k].keys())[0]
    for ln in entries:
        l=ln.strip().split(",")
        if l[dir_ind] == dict[l[kmer_ind]]['choose']:
            outfile.write(ln)
    outfile.close()

def postprocessing(args):
    diffmod_dir = args.diffmod_dir
    gtf_or_gff  = args.gtf_or_gff
    genome      = args.genome
    diffmod_table_path = os.path.join(diffmod_dir,"diffmod.table")
    if gtf_or_gff:
        annotation_dict=readAnnotation(gtf_or_gff)
        diffmod_table_path=addChrGeneNameToDiffmodTab(annotation_dict,diffmod_table_path,genome,diffmod_dir)
    run_postprocessing(diffmod_table_path,diffmod_dir)
