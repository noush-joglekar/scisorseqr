import pandas as pd
import pyranges as pr

pd.options.mode.chained_assignment = None

# helper functions - from cerberus (https://github.com/fairliereese/cerberus) - will be uploaded to pip in the future

def get_stranded_gtf_dfs(df):
    """
    Split a GTF df into fwd and rev strands
    Parameters:
        df (pandas DataFrame): DF of gtf
    Returns:
        fwd (pandas DataFrame): DF of all forward-stranded entries from GTF
        rev (pandas DataFrame): DF of all reverse-stranded entries from GTF
    """
    rev = df.loc[df.Strand == '-'].copy(deep=True)
    fwd = df.loc[df.Strand == '+'].copy(deep=True)

    return fwd, rev

def sort_gtf(df):
    """
    Sort a GTF into its proper ordering
    Parameters:
        df (pandas DataFrame): DF of GTF
    Returns:
        df (pandas DataFrame): DF of GTF, sorted
    """
    df['feature_rank'] = df.Feature.map({'gene':0, 'transcript':1, 'exon':2})
    df.feature_rank = df.feature_rank.astype(int)

    fwd, rev = get_stranded_gtf_dfs(df)

    df = pd.DataFrame()
    for temp in [fwd, rev]:
        if len(temp.index) > 0:
            strand = temp.Strand.values.tolist()[0]
            if strand == '+':
                ascending = True
            elif strand == '-':
                ascending = False
            temp.sort_values(by=['gene_id', 'transcript_id', 'feature_rank', 'Start'],
                             ascending=[True, True, True, ascending],
                             na_position='first', inplace=True)

            df = pd.concat([df, temp], ignore_index=True)
    df.drop('feature_rank', axis=1, inplace=True)
    return df

def get_update_ends_settings(strand, mode):
    """
    Returns which columns to refer to and which min/max function
    to use depending on looking at forward / rev strand or
    tss / tes
    Parameters:
        strand (str): {'+', '-'}
        mode (str): {'tss', 'tes'}
    Returns:
        old_end (str): Name of column to modify; {'Start', 'End'}
        new_end (str): Name of column to pull new value from; {'Start_end', 'End_end'}
        gene_func (str): What function to apply to new_end; {'min', 'max'}
    """
    if mode == 'tss':
        if strand == '+':
            old_end = 'Start'
            new_end = 'Start_end'
            gene_func = 'min'
        elif strand == '-':
            old_end = 'End'
            new_end = 'End_end'
            gene_func = 'max'
    elif mode == 'tes':
        if strand == '+':
            old_end = 'End'
            new_end = 'End_end'
            gene_func = 'max'
        elif strand == '-':
            old_end = 'Start'
            new_end  = 'Start_end'
            gene_func = 'min'

    return old_end, new_end, gene_func

def make_hier_entry(df, how='t'):
    """
    Generate entries for each gene and transcript using 
    the exons that belong to each of them.
    
    Parameters:
        df (pandas DataFrame): GTF DF with exons only
        how (str): {'t', 'g'}; for transcripts or genes
            respectively
    
    Returns:
        t_df (pandas DataFrame): DF with transcripts
            or genes added
    """
    agg_dict = {'min_coord': 'min', 'max_coord': 'max'}
    t_df = df.copy(deep=True)
    t_df['min_coord'] = t_df[['Start', 'End']].min(axis=1)
    t_df['max_coord'] = t_df[['Start', 'End']].max(axis=1)
    if how == 't':
        gb_cols = ['Chromosome', 'Strand', 'gene_name',
                   'gene_id', 'transcript_id', 'transcript_name',
                   'tss_id', 'tes_id',
                   'new_transcript_id', 'original_transcript_id',
                   'original_transcript_name', 'ag1', 'ag2']
        gb_cols = list(set(gb_cols)&(set(t_df.columns)))
    elif how == 'g':
        gb_cols = ['Chromosome', 'Strand', 'gene_name',
                   'gene_id']
        gb_cols = list(set(gb_cols)&(set(t_df.columns)))

    cols = gb_cols + ['min_coord', 'max_coord']
    t_df = t_df[cols]
    t_df = t_df.groupby(gb_cols).agg(agg_dict).reset_index()
    t_df.rename({'min_coord': 'Start', 'max_coord': 'End'}, axis=1, inplace=True)
    if how == 't':
        t_df['Feature'] = 'transcript'
    elif how == 'g':
        t_df['Feature'] = 'gene'

    return t_df

def make_gtf_from_exons(df):
    """
    Turn a GTF with just exon entries into a GTF with 
    exons, transcripts, and genes
    
    Parameters:
        df (pandas DataFrame): DF representation of a GTF
            with just exon entries
    
    Returns:
        df (pandas DataFrame): DF representation of a GTF with 
            exon, transcript, and gene entries
    """

    # make transcript entries
    t_df = make_hier_entry(df, how='t')

    # make gene entries
    g_df = make_hier_entry(df, how='g')

    # concat everything and sort by gene id, transcript id, feature rank (gene =0, t =1, exon=2), then start coords
    df = pd.concat([df, t_df, g_df])
    df = sort_gtf(df)
    return df

def tsv_to_gtf(tsv,
               ofile):
    """
    Convert a scisorseq-style list of TSS regions, polyA regions, and intron chains 
    into a GTF.
    
    Parameters:
        tsv (str): Path to file to convert
        gtf (str): Path to output GTF to save
    """
    df = pd.read_csv(tsv, sep='\t', header=None,
                names=['transcript_id', 'original_transcript_id',
                       'tss', 'ic', 'tes'])

    # pull chromosome and strand out from tss
    df[['Chromosome', 'Strand']] = df.tss.str.split('_', expand=True)[[0,3]]

    # get intron chain coordinates
    df['coords'] = df.ic.str.split(';%;')
    df.drop('ic', axis=1, inplace=True)
    df = df.explode('coords', ignore_index=True)

    # remove weird blank entries
    df = df.loc[df.coords != '']

    # get start and stop of each intron
    df[['Start', 'End']] = df.coords.str.split('_', expand=True)[[1,2]]
    df[['Start', 'End']] = df[['Start', 'End']].astype(int)
    df.drop('coords', axis=1, inplace=True)

    # increment / decrement coords to turn to exons
    df.Start = df.Start-1
    df.End = df.End+1

    # melt to get coords individually
    df = df.melt(id_vars=['transcript_id', 'original_transcript_id',
                     'tss', 'tes', 'Chromosome', 'Strand'],
                     value_vars=['Start', 'End'],
                     value_name='Coord')
    df.drop('variable', axis=1, inplace=True)
    df.Coord = df.Coord.astype(int)

    # sort based on strand
    fwd, rev = get_stranded_gtf_dfs(df)
    fwd.sort_values(by=['Chromosome', 'transcript_id', 'Coord'],
                    ascending=[True, True, True], inplace=True)
    rev.sort_values(by=['Chromosome', 'transcript_id', 'Coord'],
                    ascending=[True, True, False], inplace=True)

    # get start / end coords of tss / tes regions
    df[['tss_Start', 'tss_End']] = df.tss.str.split('_', expand=True)[[1,2]].astype(int)
    df[['tes_Start', 'tes_End']] = df.tes.str.split('_', expand=True)[[1,2]].astype(int)

    fwd, rev = get_stranded_gtf_dfs(df)

    # loop through strands and ends
    ends = pd.DataFrame()
    for mode in ['tss', 'tes']:
        for strand, temp in zip(['+', '-'], [fwd, rev]):
            end_cols = ['{}_Start'.format(mode),
                        '{}_End'.format(mode)]
            keep_cols = ['transcript_id', 'original_transcript_id',
                         'Chromosome', 'Strand']
            keep_cols += end_cols
            temp = temp[keep_cols]
            old_end, new_end, func = get_update_ends_settings(strand, mode)
            temp['Coord'] = temp[end_cols].apply(func, axis=1)
            temp.drop(end_cols, axis=1, inplace=True)
            temp.drop_duplicates(inplace=True)
            ends = pd.concat([ends, temp])


    # drop unnecessary columns
    df.drop(['tss', 'tes', 'tss_Start',
             'tes_Start', 'tss_End', 'tes_End'],
            axis=1, inplace=True)

    # concatenate and then sort based on strand
    df = pd.concat([df, ends])
    fwd, rev = get_stranded_gtf_dfs(df)
    fwd.sort_values(by=['Chromosome', 'transcript_id', 'Coord'],
                    ascending=[True, True, True], inplace=True)
    rev.sort_values(by=['Chromosome', 'transcript_id', 'Coord'],
                    ascending=[True, True, False], inplace=True)
    df = pd.concat([fwd, rev])

    # add gene id
    df['gene_id'] = df.transcript_id.str.split('_', expand=True)[0]

    # pivot from individual coordinates to pairs of coordinates for each exon
    fwd, rev = get_stranded_gtf_dfs(df)
    df = pd.DataFrame()
    for strand, temp in zip(['+', '-'], [fwd, rev]):
        if strand == '+':
            start_ind = 0
            end_ind = 1
        elif strand == '-':
            start_ind = 1
            end_ind = 0

        exons = pd.DataFrame({'Start':temp['Coord'].iloc[start_ind::2].values,
                              'Start_tid':temp['transcript_id'].iloc[start_ind::2].values,
                              'End':temp['Coord'].iloc[end_ind::2].values,
                              'End_tid':temp['transcript_id'].iloc[end_ind::2].values})
        temp = temp.iloc[::2]
        temp.drop('Coord', axis=1, inplace=True)
        exons.drop(['Start_tid','End_tid'], axis=1, inplace=True)

        temp['Start'] = exons['Start'].tolist()
        temp['End'] = exons['End'].tolist()

        df = pd.concat([df, temp])


    # add exon feature tag
    df['Feature'] = 'exon'

    # add gene and transcript entries
    df = make_gtf_from_exons(df)

    # decrement the start coordinate to convert bed style coordinates
    # to gtf style coordinates
    df['Start'] = df['Start']-1

    # write
    df = pr.PyRanges(df)
    df.to_gtf(ofile)
