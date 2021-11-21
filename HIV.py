# LIBRARIES

import numpy as np
import pandas as pd

# GLOBAL VARIABLES

NUC = ['-', 'A', 'C', 'G', 'T']
PRO = ['-', 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H',
       'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
REF = NUC[0]
CONS_TAG = 'CONSENSUS'
HXB2_TAG = 'B.FR.1983.HXB2-LAI-IIIB-BRU.K03455.19535'
TIME_INDEX = 3
ALPHABET = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ+++++++++++++++++++++++++++'

## Code Ocean directories
#HIV_DIR = '../data/HIV'

# GitHub directories
HIV_DIR = 'data/HIV'


# FUNCTIONS

def index2frame(i):
    """ Return the open reading frames corresponding to a given HXB2 index. """

    frames = []
    
    if ( 790<=i<=2292) or (5041<=i<=5619) or (8379<=i<=8469) or (8797<=i<=9417):
        frames.append(1)
    if (5831<=i<=6045) or (6062<=i<=6310) or (8379<=i<=8653):
        frames.append(2)
    if (2253<=i<=5096) or (5559<=i<=5850) or (5970<=i<=6045) or (6225<=i<=8795):
        frames.append(3)

    return frames


def codon2aa(c, noq=False):
    """ Return the amino acid character corresponding to the input codon. """
    
    # If all nucleotides are missing, return gap
    if c[0]=='-' and c[1]=='-' and c[2]=='-': return '-'
    
    # Else if some nucleotides are missing, return '?'
    elif c[0]=='-' or c[1]=='-' or c[2]=='-':
        if noq: return '-'
        else:   return '?'
    
    # If the first or second nucleotide is ambiguous, AA cannot be determined, return 'X'
    elif c[0] in ['W', 'S', 'M', 'K', 'R', 'Y'] or c[1] in ['W', 'S', 'M', 'K', 'R', 'Y']: return 'X'
    
    # Else go to tree
    elif c[0]=='T':
        if c[1]=='T':
            if    c[2] in ['T', 'C', 'Y']: return 'F'
            elif  c[2] in ['A', 'G', 'R']: return 'L'
            else:                          return 'X'
        elif c[1]=='C':                    return 'S'
        elif c[1]=='A':
            if    c[2] in ['T', 'C', 'Y']: return 'Y'
            elif  c[2] in ['A', 'G', 'R']: return '*'
            else:                          return 'X'
        elif c[1]=='G':
            if    c[2] in ['T', 'C', 'Y']: return 'C'
            elif  c[2]=='A':               return '*'
            elif  c[2]=='G':               return 'W'
            else:                          return 'X'
        else:                              return 'X'
    
    elif c[0]=='C':
        if   c[1]=='T':                    return 'L'
        elif c[1]=='C':                    return 'P'
        elif c[1]=='A':
            if    c[2] in ['T', 'C', 'Y']: return 'H'
            elif  c[2] in ['A', 'G', 'R']: return 'Q'
            else:                          return 'X'
        elif c[1]=='G':                    return 'R'
        else:                              return 'X'
    
    elif c[0]=='A':
        if c[1]=='T':
            if    c[2] in ['T', 'C', 'Y']: return 'I'
            elif  c[2] in ['A', 'M', 'W']: return 'I'
            elif  c[2]=='G':               return 'M'
            else:                          return 'X'
        elif c[1]=='C':                    return 'T'
        elif c[1]=='A':
            if    c[2] in ['T', 'C', 'Y']: return 'N'
            elif  c[2] in ['A', 'G', 'R']: return 'K'
            else:                          return 'X'
        elif c[1]=='G':
            if    c[2] in ['T', 'C', 'Y']: return 'S'
            elif  c[2] in ['A', 'G', 'R']: return 'R'
            else:                          return 'X'
        else:                              return 'X'
    
    elif c[0]=='G':
        if   c[1]=='T':                    return 'V'
        elif c[1]=='C':                    return 'A'
        elif c[1]=='A':
            if    c[2] in ['T', 'C', 'Y']: return 'D'
            elif  c[2] in ['A', 'G', 'R']: return 'E'
            else:                          return 'X'
        elif c[1]=='G':                    return 'G'
        else:                              return 'X'
    
    else:                                  return 'X'


def get_MSA(ref, noArrow=True):
    """Take an input FASTA file and return the multiple sequence alignment, along with corresponding tags. """
    
    temp_msa = [i.split() for i in open(ref).readlines()]
    temp_msa = [i for i in temp_msa if len(i)>0]
    
    msa = []
    tag = []
    
    for i in temp_msa:
        if i[0][0]=='>':
            msa.append('')
            if noArrow: tag.append(i[0][1:])
            else: tag.append(i[0])
        else: msa[-1]+=i[0]
    
    msa = np.array(msa)
    
    return msa, tag


def save_MSA(msa, tag, out, fasta_width=100):
    """ Write a multiple sequence alignment with corresponding tags as a FASTA file. """
    
    f = open(out+'.fasta','w')
    
    for i in range(len(msa)):
        f.write('>'+tag[i]+'\n')
        count=0
        while count<len(msa[i]):
            for j in range(fasta_width):
                f.write(msa[i][count])
                count += 1
                if count==len(msa[i]):
                    break
            f.write('\n')
    f.close()


def filter_sequences(ppt, seq_path, check_time=True):
    """ Filter out sequences without corresponding collection times and re-save them in a readable format. """
    
    # Read in list of accessions and corresponding times
    
    acc2time = {}
    d = [i.split() for i in open('%s/raw/%s-accession2time.dat' % (HIV_DIR, ppt)).readlines()]
    
    for i in d:
        acc2time[i[0]] = int(i[1])
    
    # Read in msa sequences
    
    msa, tag = get_MSA(seq_path, noArrow=True)
    msa, tag = list(msa), list(tag)
    
    HXB2_idx = tag.index(HXB2_TAG)
    HXB2_seq = msa[HXB2_idx]
    del msa[HXB2_idx]
    del tag[HXB2_idx]
    
    cons_idx = tag.index(CONS_TAG)
    cons_seq = msa[cons_idx]
    del msa[cons_idx]
    del tag[cons_idx]

    # Current tag format: (Clade).(Location)  .(Year)      .(Sequence name).(Accession)  .(Patient ID)
    # New tag format:     (Clade).(Patient ID).(Visit time).(Sample time)  .(Sequence ID)

    new_msa = []
    new_tag = []

    for i in range(len(msa)):
        clade = tag[i].split('.')[0]
        acc   = tag[i].split('.')[-2]
        if check_time:
            if acc in acc2time:
                time = acc2time[acc]
                
                if time!=-1:
                    new_msa.append(msa[i])
                    new_tag.append('.'.join([clade, ppt, 'x', str(time), acc]))
            else:
                print('No time found for %s' % acc)
        else:
            new_msa.append(msa[i])
            new_tag.append('.'.join([clade, ppt, 'x', '0', acc]))

    save_MSA([HXB2_seq, cons_seq]+new_msa, [HXB2_TAG, CONS_TAG]+new_tag, HIV_DIR+'/interim/'+seq_path.split('/')[-1].split('.')[0]+'-filtered')


def clip_MSA(HXB2_start, HXB2_end, msa, tag):
    """ Clip the input MSA to the specified range of HXB2 indices and return. """
    
    align_start = 0
    align_end = 0
    HXB2_index = tag.index(HXB2_TAG)
    HXB2_seq = msa[HXB2_index]
    HXB2_count = 0
    for i in range(len(HXB2_seq)):
        if HXB2_seq[i]!='-':
            HXB2_count += 1
        if HXB2_count==HXB2_start:
            align_start = i
        if HXB2_count==HXB2_end+1:
            align_end = i
    return np.array([np.array(list(s[align_start:align_end].upper())) for s in msa])


def filter_excess_gaps(msa, tag, sequence_max_gaps, site_max_gaps, verbose=True):
    """ Remove sequences and sites from the alignment which have excess gaps. """
    
    msa = list(msa)
    tag = list(tag)
    
    HXB2_idx = tag.index(HXB2_TAG)
    HXB2_seq = msa[HXB2_idx]
    del msa[HXB2_idx]
    del tag[HXB2_idx]
    
    cons_idx = tag.index(CONS_TAG)
    cons_seq = msa[cons_idx]
    del msa[cons_idx]
    del tag[cons_idx]
    
    # Remove sequences with too many gaps
    temp_msa = []
    temp_tag = []
    cons_gaps = np.sum(cons_seq=='-')
    for i in range(len(msa)):
        if np.sum(msa[i]=='-')-cons_gaps<sequence_max_gaps:
            temp_msa.append(msa[i])
            temp_tag.append(tag[i])
    temp_msa = np.array(temp_msa)
    if verbose:
        print('\tselected %d of %d sequences with <%d gaps in excess of consensus' %
              (len(temp_msa), len(msa), sequence_max_gaps))
    
    # Drop sites that have too many gaps
    kept_indices = []
    for i in range(len(HXB2_seq)):
        if HXB2_seq[i]!='-' or np.sum(temp_msa[:,i]=='-')/len(temp_msa)<site_max_gaps:
            kept_indices.append(i)
    temp_msa = np.array([HXB2_seq[kept_indices], cons_seq[kept_indices]] + [s[kept_indices] for s in temp_msa])
    temp_tag = [HXB2_TAG, CONS_TAG] + temp_tag
    if verbose:
        print('\tremoved %d of %d sites with >%d%% gaps' %
              (len(msa[0])-len(kept_indices), len(msa[0]), site_max_gaps*100))

    return temp_msa, temp_tag


def get_times(msa, tag, sort=False):
    """Return sequences and times collected from an input MSA and tags (optional: time order them)."""

    times = []
    for i in range(len(tag)):
        if tag[i] not in [HXB2_TAG, CONS_TAG]:
            tsplit = tag[i].split('.')
            times.append(int(tsplit[TIME_INDEX]))
        else:
            times.append(-1)
    
    if sort:
        t_sort = np.argsort(times)
        return np.array(msa)[t_sort], np.array(tag)[t_sort], np.array(times)[t_sort]

    else:
        return np.array(times)


def order_sequences(msa, tag):
    """ Put sequences in time order. """

    msa = list(msa)
    tag = list(tag)
    
    HXB2_idx = tag.index(HXB2_TAG)
    HXB2_seq = msa[HXB2_idx]
    del msa[HXB2_idx]
    del tag[HXB2_idx]
    
    cons_idx = tag.index(CONS_TAG)
    cons_seq = msa[cons_idx]
    del msa[cons_idx]
    del tag[cons_idx]
    
    temp_msa = [HXB2_seq, cons_seq]
    temp_tag = [HXB2_TAG, CONS_TAG]
    msa, tag, temp = get_times(msa, tag, sort=True)
    
    return np.array(temp_msa + list(msa)), np.array(temp_tag + list(tag))


def impute_ambiguous(msa, tag, start_index=0, verbose=True, impute_edge_gaps=False):
    """ Impute ambiguous nucleotides with the most frequently observed ones in the alignment. """
    
    # Impute ambiguous nucleotides
    for i in range(len(msa[0])):
        for j in range(start_index, len(msa)):
            orig = msa[j][i].upper()
            if orig not in NUC:
                avg = [np.sum([msa[k][i]==a for k in range(start_index, len(msa))]) for a in NUC]
                new = NUC[np.argmax(avg)]
                if orig=='R': # A or G
                    if avg[NUC.index('A')]>avg[NUC.index('G')]:
                        new = 'A'
                    else:
                        new = 'G'
                elif orig=='Y': # T or C
                    if avg[NUC.index('T')]>avg[NUC.index('C')]:
                        new = 'T'
                    else:
                        new = 'C'
                elif orig=='K': # G or T
                    if avg[NUC.index('G')]>avg[NUC.index('T')]:
                        new = 'G'
                    else:
                        new = 'T'
                elif orig=='M': # A or C
                    if avg[NUC.index('A')]>avg[NUC.index('C')]:
                        new = 'A'
                    else:
                        new = 'C'
                elif orig=='S': # G or C
                    if avg[NUC.index('G')]>avg[NUC.index('C')]:
                        new = 'G'
                    else:
                        new = 'C'
                elif orig=='W': # A or T
                    if avg[NUC.index('A')]>avg[NUC.index('T')]:
                        new = 'A'
                    else:
                        new = 'T'
                msa[j][i] = new
                if verbose:
                    print('\texchanged %s for %s in sequence %d, site %d' % (new, orig, j, i))
                    
    # Impute leading and trailing gaps
    if impute_edge_gaps:
        for j in range(start_index, len(msa)):
            gap_lead = 0
            gap_trail = 0
            
            gap_idx = 0
            while msa[j][gap_idx]=='-':
                gap_lead += 1
                gap_idx += 1
                
            gap_idx = -1
            while msa[j][gap_idx]=='-':
                gap_trail -= 1
                gap_idx -= 1
                
            for i in range(gap_lead):
                avg = [np.sum([msa[k][i]==a for k in range(start_index, len(msa))]) for a in NUC]
                new = NUC[np.argmax(avg)]
                msa[j][i] = new
                
            for i in range(gap_trail, 0):
                avg = [np.sum([msa[k][i]==a for k in range(start_index, len(msa))]) for a in NUC]
                new = NUC[np.argmax(avg)]
                msa[j][i] = new
                
            if (gap_lead>0) or (gap_trail<0):
                print('\timputed %d leading gaps and %d trailing gaps in sequence %d' % (gap_lead, -1*gap_trail, j))

    return msa


def get_TF(msa, tag, TF_accession, protein=False):
    """ Return the transmitted/founder sequence in an alignment. If there is no known TF sequence,
        return the most frequently observed nucleotide at each site from the earliest available sequences. """
    
    TF_sequence = []

    if TF_accession=='avg':
        idxs = [i for i in range(len(msa)) if tag[i]!=HXB2_TAG and tag[i]!=CONS_TAG]
        temp_msa = np.array(msa)[idxs]
        temp_tag = np.array(tag)[idxs]

        times = get_times(temp_msa, temp_tag, sort=False)
        first_time = np.min(times)
        first_seqs = [temp_msa[i] for i in range(len(temp_msa)) if times[i]==first_time]
        for i in range(len(first_seqs[0])):
            if protein:
                avg = [np.sum([s[i]==a for s in first_seqs]) for a in PRO]
                TF_sequence.append(PRO[np.argmax(avg)])
            else:
                avg = [np.sum([s[i]==a for s in first_seqs]) for a in NUC]
                TF_sequence.append(NUC[np.argmax(avg)])

    else:
        accs = [i.split('.')[-1] for i in tag]
        TF_sequence = msa[accs.index(TF_accession)]

    return TF_sequence


def create_index(msa, tag, TF_seq, cons_seq, HXB2_seq, HXB2_start, min_seqs, max_dt, df_epitope, df_exposed, out_file,
                 return_polymorphic=True, return_truncated=True):
    """ Create a reference to map between site indices for the whole alignment, polymorphic sites only, and HXB2.
        To preserve quality, identify last time point such that all earlier time points have at least min_seqs
        sequences (except time 0, which is allowed to have 1=TF) and maximum time gap of max_dt between samples.
        Include location of known epitopes, flanking residues for those epitopes, and exposed regions of Env.
        Also record the TF and consensus nucleotides at each site. Return the list of polymorphic sites. """

    msa = list(msa)
    tag = list(tag)
    
    HXB2_idx = tag.index(HXB2_TAG)
    HXB2_seq = msa[HXB2_idx]
    del msa[HXB2_idx]
    del tag[HXB2_idx]
    
    cons_idx = tag.index(CONS_TAG)
    cons_seq = msa[cons_idx]
    del msa[cons_idx]
    del tag[cons_idx]
    
    f = open('%s' % out_file, 'w')
    f.write('alignment,polymorphic,HXB2,TF,consensus,epitope,exposed,edge_gap,flanking\n')
    
    # Check for minimum number of sequences/maximum dt to truncate alignment
    temp_msa, temp_tag, times = get_times(msa, tag, sort=True)
    u_times = np.unique(times)
    t_count = [np.sum(times==t) for t in u_times]
    
#    print('\t'.join([str(int(i)) for i in u_times]))
#    print('\t'.join([str(int(i)) for i in t_count]))
#    t_max = 0
#    for i in range(1, len(t_count)):
#        if t_count[i]<min_seqs or u_times[i]-u_times[i-1]>max_dt:
#            break
#        else:
#            t_max += 1
#    t_max = u_times[t_max]
#    temp_msa = temp_msa[times<=t_max]
#    temp_tag = temp_tag[times<=t_max]

    t_allowed = [u_times[0]]
    t_last    = u_times[0]
    for i in range(1, len(t_count)):
        if t_count[i]<min_seqs:
            continue
        elif u_times[i]-t_last>max_dt:
            break
        else:
            t_allowed.append(u_times[i])
            t_last = u_times[i]
    t_max    = t_allowed[-1]
    temp_msa = temp_msa[np.isin(times, t_allowed)]
    temp_tag = temp_tag[np.isin(times, t_allowed)]
    
    HXB2_index = HXB2_start
    polymorphic_index = 0
    polymorphic_sites = []
    for i in range(len(temp_msa[0])):
        
        # Index polymorphic sites
        poly_str = 'NA'
        if np.sum([s[i]==temp_msa[0][i] for s in temp_msa])<len(temp_msa):
            poly_str = '%d' % polymorphic_index
            polymorphic_index += 1
            polymorphic_sites.append(i)
        
        # Index HXB2
        HXB2_str = 'NA'
        if HXB2_seq[i]!='-':
            HXB2_str = '%d' % HXB2_index
            HXB2_index += 1
            HXB2_alpha  = 0
        else:
            HXB2_str = '%d%s' % (HXB2_index-1, ALPHABET[HXB2_alpha])
            HXB2_alpha += 1
        
        # Flag epitope regions
        epitope_str = ''
        flanking = 0
        for epitope_iter, epitope_entry in df_epitope.iterrows():
            if (HXB2_index-1>=epitope_entry.start-15 and HXB2_index-1<epitope_entry.start
                and epitope_entry.detected<=t_max):
                flanking += 1
            elif (HXB2_index-1<=epitope_entry.end+15 and HXB2_index-1>epitope_entry.end
                  and epitope_entry.detected<=t_max):
                flanking += 1
            if (HXB2_index-1>=epitope_entry.start and HXB2_index-1<=epitope_entry.end
                and epitope_entry.detected<=t_max):
                epitope_str = epitope_entry.epitope
            # special case: first 3 AA inserted wrt HXB2
            elif epitope_entry.epitope=='DEPAAVGVG':
                if (i>=3870 and HXB2_index-1<=epitope_entry.end
                    and epitope_entry.detected<=t_max):
                    epitope_str = epitope_entry.epitope
        
        # Flag exposed sites on Env
        exposed = False
        if np.sum((HXB2_index-1>=df_exposed.start) & (HXB2_index-1<=df_exposed.end))>0:
            exposed = True
            
        # Flag edge gaps
        edge_def = 200
        edge_gap = False
        if np.sum(temp_msa[:, i]=='-')>0 and ((i<edge_def) or (len(temp_msa[0])-i<edge_def)):
            gap_seqs = [j for j in range(len(temp_msa)) if temp_msa[j][i]=='-']
            gap_msa = temp_msa[gap_seqs]
            edge_gap = True
            if i<edge_def:
                for s in gap_msa:
                    if np.sum(s[:i]=='-')<i:
                        edge_gap = False
                        break
            else:
                for s in gap_msa:
                    if np.sum(s[i:]=='-')<len(temp_msa[0])-i:
                        edge_gap = False
                        break

        # Save to file
        f.write('%d,%s,%s,%s,%s,%s,%s,%s,%d\n' % (i, poly_str, HXB2_str, TF_seq[i], cons_seq[i], epitope_str, exposed, edge_gap, flanking))
    f.close()

    temp_msa = [HXB2_seq, cons_seq] + list(temp_msa)
    temp_tag = [HXB2_TAG, CONS_TAG] + list(temp_tag)
    
    if return_polymorphic and return_truncated:
        return polymorphic_sites, temp_msa, temp_tag
    elif return_polymorphic:
        return polymorphic_sites
    elif return_truncated:
        return temp_msa, temp_tag


def save_MPL_alignment(msa, tag, out_file, polymorphic_sites=[], return_states=True, protein=False):
    """ Save a nucleotide alignment into MPL-readable form. Optionally return converted states and times. """

    idxs = [i for i in range(len(msa)) if tag[i]!=HXB2_TAG and tag[i]!=CONS_TAG]
    temp_msa = np.array(msa)[idxs]
    temp_tag = np.array(tag)[idxs]
    
    if polymorphic_sites==[]:
        polymorphic_sites = range(len(temp_msa[0]))

    poly_times = get_times(temp_msa, temp_tag, sort=False)

    poly_states = []
    if protein:
        for s in temp_msa:
            poly_states.append([str(PRO.index(a)) for a in s[polymorphic_sites]])
    else:
        for s in temp_msa:
            poly_states.append([str(NUC.index(a)) for a in s[polymorphic_sites]])

    f = open(out_file, 'w')
    for i in range(len(poly_states)):
        f.write('%d\t1\t%s\n' % (poly_times[i], ' '.join(poly_states[i])))
    f.close()

    if return_states:
        return np.array(poly_states, int), np.array(poly_times)


def get_effective_HXB2_index(start, df_index):
    """ Obtain an effective HXB2 index for sites that are inserted relative to HXB2. """

    index = 0
    shift = 0
    found = False
    while not found and shift<start:
        shift += 1
        i = df_index.iloc[start-shift]
        if pd.notnull(i.HXB2):
            found = True
            index = int(i.HXB2) + 1
    if not found:
        print('Never found HXB2 index')

    return index, shift


def get_nonsynonymous(polymorphic_sites, nuc, i, i_HXB2, shift, frames, TF_sequence, match_states, verbose=True):
    """ Return number of reading frames in which the input nucleotide is nonsynonymous in context, compared to T/F. """
    
    ns = 0
    for fr in frames:
        
        pos = int((i_HXB2+shift-fr)%3) # position of the nucleotide in the reading frame
        TF_codon = TF_sequence[i-pos:i-pos+3]
        
        if len(TF_codon)<3 and verbose:
            print('\tmutant at site %d in codon that does not terminate in alignment, assuming syn' % i)
        
        else:
            mut_codon = [a for a in TF_codon]
            mut_codon[pos] = nuc
            replace_indices = [k for k in range(3) if (k+i-pos) in polymorphic_sites and k!=pos]
            
            # If any other sites in the codon are polymorphic, consider mutation in context
            if len(replace_indices)>0:
                is_ns = False
                for s in match_states:
                    TF_codon = TF_sequence[i-pos:i-pos+3]
                    for k in replace_indices:
                        mut_codon[k] = NUC[s[polymorphic_sites.index(k+i-pos)]]
                        TF_codon[k] = NUC[s[polymorphic_sites.index(k+i-pos)]]
                    if codon2aa(mut_codon)!=codon2aa(TF_codon):
                        is_ns = True
                if is_ns:
                    ns += 1
        
            elif codon2aa(mut_codon)!=codon2aa(TF_codon):
                ns += 1

    return ns
    

def get_nonsynonymous_alternate(polymorphic_sites, nuc, i, i_HXB2, shift, frames, TF_sequence, match_states, verbose=True):
    """ Return number of reading frames in which the input nucleotide is nonsynonymous in context, compared to T/F. """
    
    ns        = [-1, -1, -1]
    positions = [-1, -1, -1]
    for fr in frames:
        ns[fr-1] = 0
        
        pos = int((i_HXB2+shift-fr)%3) # position of the nucleotide in the reading frame
        positions[fr-1] = pos
        TF_codon = TF_sequence[i-pos:i-pos+3]
        
        if len(TF_codon)<3 and verbose:
            print('\tmutant at site %d in codon that does not terminate in alignment, assuming syn' % i)
        
        else:
            mut_codon = [a for a in TF_codon]
            mut_codon[pos] = nuc
            replace_indices = [k for k in range(3) if (k+i-pos) in polymorphic_sites and k!=pos]
            
            # If any other sites in the codon are polymorphic, consider mutation in context
            if len(replace_indices)>0:
                is_ns = False
                for s in match_states:
                    TF_codon = TF_sequence[i-pos:i-pos+3]
                    for k in replace_indices:
                        mut_codon[k] = NUC[s[polymorphic_sites.index(k+i-pos)]]
                        TF_codon[k] = NUC[s[polymorphic_sites.index(k+i-pos)]]
                    if codon2aa(mut_codon)!=codon2aa(TF_codon):
                        is_ns = True
                if is_ns:
                    ns[fr-1] = 1
        
            elif codon2aa(mut_codon)!=codon2aa(TF_codon):
                ns[fr-1] = 1

    return ns, positions


def get_glycosylation(polymorphic_sites, nuc, i, i_HXB2, shift, TF_sequence, match_states):
    """ Return net addition/subtraction of N-linked glycosylation motifs (N-X-S/T) in context, compared with T/F. """
    
    notX = ['P', '*', '?', '-']

    if i_HXB2>=6225 and i_HXB2<=8795:
        pos = int((i_HXB2+shift)%3)
        TF_codon = TF_sequence[i-pos:i-pos+3]
        glycan_results = []
        
        if len(TF_codon)==3:
            codon_idx = [k for k in range(3) if (k+i-pos) in polymorphic_sites and k!=pos]
            mut_codon = [a for a in TF_codon]
            mut_codon[pos] = nuc
            
            for s in match_states:
                temp_glycan = 0
                
                # Get codons in context
                for k in codon_idx:
                    mut_codon[k] = NUC[s[polymorphic_sites.index(k+i-pos)]]
                    TF_codon[k] = NUC[s[polymorphic_sites.index(k+i-pos)]]
                mut_aa = codon2aa(mut_codon)
                TF_aa = codon2aa(TF_codon)
                
                if mut_aa!=TF_aa:
                    
                    # Get +/- 2 codons around current site
                    codon_shifts = [-6, -3, 3, 6]
                    match_codons = [list(TF_sequence[i-pos+k:i-pos+3+k]) for k in codon_shifts]
                    for k in range(len(codon_shifts)):
                        for kk in range(3):
                            if (i-pos)+codon_shifts[k]+kk in polymorphic_sites:
                                match_codons[k][kk] = NUC[s[polymorphic_sites.index((i-pos)+codon_shifts[k]+kk)]]
                    
                    # + glycan site
                    
                    if (mut_aa=='N') and len(match_codons[2])==3 and len(match_codons[3])==3:
                        TF_motif = [codon2aa(match_codons[2]), codon2aa(match_codons[3])]
                        if (TF_motif[0] not in notX) and (TF_motif[1] in ['S', 'T']):
                            temp_glycan += 1
                    
                    if (mut_aa not in notX) and (TF_aa in notX) and len(match_codons[1])==3 and len(match_codons[2])==3:
                        TF_motif = [codon2aa(match_codons[1]), codon2aa(match_codons[2])]
                        if (TF_motif[0]=='N') and (TF_motif[1] in ['S', 'T']):
                            temp_glycan += 1
                        
                    if (mut_aa in ['S', 'T']) and (TF_aa not in ['S', 'T']) and len(match_codons[0])==3 and len(match_codons[1])==3:
                        TF_motif = [codon2aa(match_codons[0]), codon2aa(match_codons[1])]
                        if (TF_motif[0]=='N') and (TF_motif[1] not in notX):
                            temp_glycan += 1
                
                    # - glycan site
                    
                    if (TF_aa=='N') and len(match_codons[2])==3 and len(match_codons[3])==3:
                        TF_motif = [codon2aa(match_codons[2]), codon2aa(match_codons[3])]
                        if (TF_motif[0] not in notX) and (TF_motif[1] in ['S', 'T']):
                            temp_glycan -= 1

                    if (TF_aa not in notX) and (mut_aa in notX) and len(match_codons[1])==3 and len(match_codons[2])==3:
                        TF_motif = [codon2aa(match_codons[1]), codon2aa(match_codons[2])]
                        if (TF_motif[0]=='N') and (TF_motif[1] in ['S', 'T']):
                            temp_glycan -= 1

                    if (TF_aa in ['S', 'T']) and (mut_aa not in ['S', 'T']) and len(match_codons[0])==3 and len(match_codons[1])==3:
                        TF_motif = [codon2aa(match_codons[0]), codon2aa(match_codons[1])]
                        if (TF_motif[0]=='N') and (TF_motif[1] not in notX):
                            temp_glycan -= 1
                                    
                glycan_results.append(temp_glycan)
        
        if len(glycan_results)==0:
            return 0
        else:
            return np.median(glycan_results)

    else:
        return 0


def save_trajectories(sites, states, times, TF_sequence, df_index, out_file, protein=False):
    """ Save allele frequency trajectories and supplementary information. """

    index_cols = ['alignment', 'polymorphic', 'HXB2']
    cols = [i for i in list(df_index) if i not in index_cols]
    
    f = open(out_file, 'w')
    f.write('polymorphic_index,alignment_index,HXB2_index,nonsynonymous,nucleotide')
    f.write(',%s,glycan' % (','.join(cols)))
    f.write(',%s\n' % (','.join(['f_at_%d' % t for t in np.unique(times)])))
    
    for i in sites:
        if protein:
            for j in range(len(PRO)):
                traj = []
                for t in np.unique(times):
                    tid = times==t
                    num = np.sum(states[tid].T[sites.index(i)]==j)
                    denom = np.sum(tid)
                    traj.append(num/denom)
                
                # NOTE: treatment of glycans, edge gaps for proteins is incomplete
                if np.sum(traj)!=0:
                    ii = df_index.iloc[i]
                    match_states = states[states.T[sites.index(i)]==j]
                    nonsyn = True
                    CpG = 0
                    glycan = 0
                    f.write('%d,%d,%s,%d,%s' % (sites.index(i), i, str(ii.HXB2), nonsyn, PRO[j]))
                    f.write(',%s,%d,%d' % (','.join([str(ii[c]) if c!='edge_gap' else str(False) for c in cols]), glycan))
                    f.write(',%s\n' % (','.join(['%.4e' % freq for freq in traj])))
    
        else:
            for j in range(len(NUC)):
                traj = []
                for t in np.unique(times):
                    tid = times==t
                    num = np.sum(states[tid].T[sites.index(i)]==j)
                    denom = np.sum(tid)
                    traj.append(num/denom)
                
                if np.sum(traj)!=0:
                    ii = df_index.iloc[i]
                    match_states = states[states.T[sites.index(i)]==j]
                    
                    # Get effective HXB2 index to determine open reading frames
                    eff_HXB2_index = 0
                    shift = 0
                    frames = []
                    try:
                        eff_HXB2_index = int(ii.HXB2)
                        frames = index2frame(eff_HXB2_index)
                    except:
                        eff_HXB2_index = int(ii.HXB2[:-1])
                        shift = ALPHABET.index(ii.HXB2[-1]) + 1
                        frames = index2frame(eff_HXB2_index)
                    
                    # Check whether mutation is nonsynonymous by inserting TF nucleotide in context
                    nonsyn = get_nonsynonymous(sites, NUC[j], i, eff_HXB2_index, shift, frames, TF_sequence, match_states)
            
                    # Flag whether variant is an edge gap
                    edge_gap = False
                    if NUC[j]=='-' and ii.edge_gap==True:
                        edge_gap = True
            
                    # If mutation is in Env, check for modification of N-linked glycosylation site (N-X-S/T motif)
                    glycan = get_glycosylation(sites, NUC[j], i, eff_HXB2_index, shift, TF_sequence, match_states)
                    
                    f.write('%d,%d,%s,%d,%s' % (sites.index(i), i, str(ii.HXB2), nonsyn, NUC[j]))
                    f.write(',%s,%d' % (','.join([str(ii[c]) if c!='edge_gap' else str(edge_gap) for c in cols]), glycan))
                    f.write(',%s\n' % (','.join(['%.4e' % freq for freq in traj])))

    f.close()
    
    
def save_trajectories_alternate(sites, states, times, TF_sequence, df_index, out_file, protein=False):
    """ Save allele frequency trajectories and supplementary information. """

    index_cols = ['alignment', 'polymorphic', 'HXB2']
    cols = [i for i in list(df_index) if i not in index_cols]
    
    f = open(out_file, 'w')
    f.write('polymorphic_index,alignment_index,HXB2_index,nucleotide,frame_1,pos_1,ns_1,frame_2,pos_2,ns_2,frame_3,pos_3,ns_3')
    f.write(',%s,glycan' % (','.join(cols)))
    f.write(',%s\n' % (','.join(['f_at_%d' % t for t in np.unique(times)])))
    
    for i in sites:
        if protein:
            for j in range(len(PRO)):
                traj = []
                for t in np.unique(times):
                    tid = times==t
                    num = np.sum(states[tid].T[sites.index(i)]==j)
                    denom = np.sum(tid)
                    traj.append(num/denom)
                
                # NOTE: treatment of glycans, edge gaps for proteins is incomplete
                if np.sum(traj)!=0:
                    ii = df_index.iloc[i]
                    match_states = states[states.T[sites.index(i)]==j]
                    nonsyn = True
                    CpG = 0
                    glycan = 0
                    f.write('%d,%d,%s,%d,%s' % (sites.index(i), i, str(ii.HXB2), nonsyn, PRO[j]))
                    f.write(',%s,%d,%d' % (','.join([str(ii[c]) if c!='edge_gap' else str(False) for c in cols]), glycan))
                    f.write(',%s\n' % (','.join(['%.4e' % freq for freq in traj])))
    
        else:
            for j in range(len(NUC)):
                traj = []
                for t in np.unique(times):
                    tid = times==t
                    num = np.sum(states[tid].T[sites.index(i)]==j)
                    denom = np.sum(tid)
                    traj.append(num/denom)
                
                if np.sum(traj)!=0:
                    ii = df_index.iloc[i]
                    match_states = states[states.T[sites.index(i)]==j]
                    
                    # Get effective HXB2 index to determine open reading frames
                    eff_HXB2_index = 0
                    shift = 0
                    frames = []
                    
                    try:
                        eff_HXB2_index = int(ii.HXB2)
                        frames = index2frame(eff_HXB2_index)
                    except:
                        eff_HXB2_index = int(ii.HXB2[:-1])
                        shift = ALPHABET.index(ii.HXB2[-1]) + 1
                        frames = index2frame(eff_HXB2_index)
                    
                    # Check whether mutation is nonsynonymous by inserting TF nucleotide in context
                    nonsyn, positions = get_nonsynonymous_alternate(sites, NUC[j], i, eff_HXB2_index, shift, frames, TF_sequence, match_states)
            
                    # Flag whether variant is an edge gap
                    edge_gap = False
                    if NUC[j]=='-' and ii.edge_gap==True:
                        edge_gap = True
            
                    # If mutation is in Env, check for modification of N-linked glycosylation site (N-X-S/T motif)
                    glycan = get_glycosylation(sites, NUC[j], i, eff_HXB2_index, shift, TF_sequence, match_states)
                    
                    f.write('%d,%d,%s,%s' % (sites.index(i), i, str(ii.HXB2), NUC[j]))
                    
                    for k in range(3):
                        f.write(',%d,%d,%d' % (k+1 in frames, positions[k], nonsyn[k]))
                    
                    f.write(',%s,%d' % (','.join([str(ii[c]) if c!='edge_gap' else str(edge_gap) for c in cols]), glycan))
                    f.write(',%s\n' % (','.join(['%.4e' % freq for freq in traj])))

    f.close()


def get_baseline(TF_sequence, df_index, df_poly):
    """ Get baseline number of possible mutations that are:
          - nonsynonymous in known CD8+ T cell epitopes
          - nonsynonymous reversions to clade consensus
          - nonsynonymous reversions to clade consensus not in T cell epitopes
          - synonymous """
    
    n_ns_epitope = 0
    n_ns_rev = 0
    n_ns_rev_nonepi = 0
    n_syn = 0
    
    for i in range(len(df_index)):
        for j in range(1,len(NUC)):
            if TF_sequence[i]!=NUC[j]:
                ii = df_index.iloc[i]
                
                # Get effective HXB2 index to determine open reading frames
                eff_HXB2_index = 0
                shift = 0
                frames = []
                try:
                    eff_HXB2_index = int(ii.HXB2)
                    frames = index2frame(eff_HXB2_index)
                except:
                    eff_HXB2_index = int(ii.HXB2[:-1])
                    shift = ALPHABET.index(ii.HXB2[-1]) + 1
                    frames = index2frame(eff_HXB2_index)
            
                # Check whether mutation is nonsynonymous by inserting mutation in TF background
                is_ns = False
                for fr in frames:
                    pos = int((eff_HXB2_index+shift-fr)%3) # position of the nucleotide in the reading frame
                    TF_codon = TF_sequence[i-pos:i-pos+3]
                    if len(TF_codon)==3:
                        mut_codon = [a for a in TF_codon]
                        mut_codon[pos] = NUC[j]
                        if codon2aa(mut_codon)!=codon2aa(TF_codon):
                            is_ns = True

                # Add to baseline values
                if is_ns and pd.notnull(ii.epitope):
                    n_ns_epitope += 1
                if is_ns and NUC[j]==ii.consensus:
                    n_ns_rev += 1
                if is_ns and NUC[j]==ii.consensus and pd.isnull(ii.epitope):
                    n_ns_rev_nonepi += 1
                if not is_ns:
                    n_syn += 1

    n_ns_epitope    += np.sum((df_poly.nucleotide=='-') & (pd.notnull(df_poly.epitope)))
    n_ns_rev        += np.sum((df_poly.nucleotide=='-') & (df_poly.TF=='-'))
    n_ns_rev_nonepi += np.sum((df_poly.nucleotide=='-') & (df_poly.TF=='-') & (pd.isnull(df_poly.epitope)))

    return n_ns_epitope, n_ns_rev, n_ns_rev_nonepi, n_syn
