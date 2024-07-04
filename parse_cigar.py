#!/usr/bin/env python
from cigar import Cigar
import pandas as pd
import sys

cigar_file = sys.argv[1]
f = open(cigar_file)
aln = f.readlines()[0]
cigar = list(Cigar(aln).items())

blocks = sys.argv[2]
synteny_transformed = pd.read_csv(blocks, header=None, sep = '\t')
synteny_transformed.columns = ["Synteny Block in Transformed Sequence", "Start Vertex in Transformed Graph", "End Vertex in Transformed Graph", "Start Pos in Transformed Sequence", "End Pos in Transformed Sequence", "Include Left Vertex", "Include Right Vertex", "Is Short Version"]

ref_len = 0
query_len = 0
query2ref = {}
for seq_len, flag in cigar:
    if flag == 'M':
        for i in range(seq_len):
            query2ref[query_len] = ref_len
            ref_len += 1
            query_len += 1
    elif flag == 'X':
        for i in range(seq_len):
            query2ref[query_len] = ref_len
            ref_len += 1
            query_len += 1
    elif flag == 'D':
        for i in range(seq_len):
            query2ref[query_len] = ref_len
            ref_len += 1
    elif flag == 'I':
        for i in range(seq_len):
            query2ref[query_len] = ref_len
            query_len += 1
    else:
        print("Unrecognized ciggar flag")

query2ref[query_len] = ref_len
# print(ref_len, query_len)

def similarity(cigar, start, end):
    ref_len = 0
    query_len = 0
    matches = 0
    mismatches = 0
    for seq_len, flag in cigar:
        if flag == 'M':
            for i in range(seq_len):
                query2ref[query_len] = ref_len
                if query_len >= start and query_len < end:
                    matches += 1
                ref_len += 1
                query_len += 1
        elif flag == 'X':
            for i in range(seq_len):
                query2ref[query_len] = ref_len
                if query_len >= start and query_len < end:
                    mismatches += 1
                ref_len += 1
                query_len += 1
        elif flag == 'D':
            for i in range(seq_len):
                query2ref[query_len] = ref_len
                if query_len >= start and query_len < end:
                    mismatches += 1
                ref_len += 1
        elif flag == 'I':
            for i in range(seq_len):
                query2ref[query_len] = ref_len
                if query_len >= start and query_len < end:
                    mismatches += 1
                query_len += 1
        else:
            print("Unrecognized ciggar flag")
    #return matches/(matches + mismatches)
    return matches

def count_J(start, end, df):
    df_new = df.loc[df["class"] == 'J']
    df_new = df_new.loc[(df_new["start"] >= start) & (df_new["end"] <= end)]
    return len(df_new)

def count_D(start, end, df):
    df_new = df.loc[df["class"] == 'D']
    df_new = df_new.loc[(df_new["start"] >= start) & (df_new["end"] <= end)]
    return len(df_new)

def count_V(start, end, df):
    df_new = df.loc[df["class"] == 'V']
    df_new = df_new.loc[(df_new["start"] >= start) & (df_new["end"] <= end)]
    return len(df_new)

synteny_transformed["Start Pos in Original Sequence"] = synteny_transformed["Start Pos in Transformed Sequence"].apply(lambda x: query2ref[x])
synteny_transformed["End Pos in Original Sequence"] = synteny_transformed["End Pos in Transformed Sequence"].apply(lambda x: query2ref[x])
synteny_transformed["Similarity of Transformed Block and Original Block"] = synteny_transformed.apply(lambda x: similarity(cigar, x["Start Pos in Transformed Sequence"], x["End Pos in Transformed Sequence"])/max(x["End Pos in Transformed Sequence"]-x["Start Pos in Transformed Sequence"],x["End Pos in Original Sequence"]-x["Start Pos in Original Sequence"], 1), axis=1)
synteny_transformed["Similarity of Transformed Block and Original Block"] = synteny_transformed["Similarity of Transformed Block and Original Block"].round(2)
synteny_transformed["Length in Original Sequence"] = synteny_transformed["End Pos in Original Sequence"] - synteny_transformed["Start Pos in Original Sequence"]

if len(sys.argv) == 5:
    gene = sys.argv[4]
    gene = pd.read_csv(gene)
    synteny_transformed["V genes"] = synteny_transformed.apply(lambda x: count_V(x["Start Pos in Original Sequence"], x["End Pos in Original Sequence"], gene), axis=1)
    synteny_transformed["D genes"] = synteny_transformed.apply(lambda x: count_D(x["Start Pos in Original Sequence"], x["End Pos in Original Sequence"], gene), axis=1)
    synteny_transformed["J genes"] = synteny_transformed.apply(lambda x: count_J(x["Start Pos in Original Sequence"], x["End Pos in Original Sequence"], gene), axis=1)


identifier = sys.argv[3]
synteny_transformed.to_csv(blocks + "." + identifier + ".csv", index=None, sep = '\t')

if len(sys.argv) == 5:
    def get_contained(start, end, df):
        df_new = df.loc[(df["Start Pos in Original Sequence"]<=start) & (df["End Pos in Original Sequence"] >= end)]
        if len(df_new):
            return df_new["Synteny Block in Transformed Sequence"].to_list()[0]
        else:
            return ""
    gene["block"] = gene.apply(lambda x: get_contained(x["start"], x["end"], synteny_transformed), axis=1)
    len(gene.loc[gene["block"] == ""])

    def get_border(start, end, df):
        df_new = df.loc[((df["Start Pos in Original Sequence"]>start) & (df["Start Pos in Original Sequence"] <= end)) | ((df["End Pos in Original Sequence"]>=start) & (df["End Pos in Original Sequence"] < end))]
        if len(df_new):
            return df_new["Synteny Block in Transformed Sequence"].to_list()
        else:
            return ""
    gene["border"] = gene.apply(lambda x: get_border(x["start"], x["end"], synteny_transformed), axis=1)
    len(gene.loc[gene["border"] != ""])

    gene.to_csv(blocks + "." + identifier + ".gene2blocks.csv", index=None, sep = '\t')