#!/usr/bin/env python
import pandas as pd
import sys

blocks1 = sys.argv[1]
blocks2 = "null"
if len(sys.argv) >= 3:
    blocks2 = sys.argv[2]

df1 = pd.read_csv(blocks1, sep = '\t')
df1["Synteny Block in Transformed Sequence"] = df1["Synteny Block in Transformed Sequence"].apply(lambda x: int(x[:x.find('(')]))
df1 = df1.loc[df1["End Pos in Original Sequence"] - df1["Start Pos in Original Sequence"] >= 2000]
df1.reset_index(drop=True, inplace=True)
df2 = pd.DataFrame()
if blocks2 != "null":
    df2 = pd.read_csv(blocks2, sep = '\t')
    df2["Synteny Block in Transformed Sequence"] = df2["Synteny Block in Transformed Sequence"].apply(lambda x: int(x[:x.find('(')]))
    df2 = df2.loc[df2["End Pos in Original Sequence"] - df2["Start Pos in Original Sequence"] >= 2000]
    df2.reset_index(drop=True, inplace=True)

block2times1 = {}
block2times2 = {}
for block in df1["Synteny Block in Transformed Sequence"]:
    if abs(block) in block2times1:
        block2times1[abs(block)] += 1
    else:
        block2times1[abs(block)] = 1
if blocks2 != "null":
    for block in df2["Synteny Block in Transformed Sequence"]:
        if abs(block) in block2times2:
            block2times2[abs(block)] += 1
        else:
            block2times2[abs(block)] = 1

length_duplication = 0
duplicons = []
for i in range(len(df1)):
    if block2times1[abs(df1["Synteny Block in Transformed Sequence"][i])] > 1:
        if df1["Synteny Block in Transformed Sequence"][i] not in duplicons:
            duplicons.append(df1["Synteny Block in Transformed Sequence"][i])
        length_duplication += (df1["End Pos in Original Sequence"][i] - df1["Start Pos in Original Sequence"][i])

length_duplication = 0
duplicons = []
for i in range(len(df2)):
    if block2times2[abs(df2["Synteny Block in Transformed Sequence"][i])] > 1:
        if df2["Synteny Block in Transformed Sequence"][i] not in duplicons:
            duplicons.append(df2["Synteny Block in Transformed Sequence"][i])
        length_duplication += (df2["End Pos in Original Sequence"][i] - df2["Start Pos in Original Sequence"][i])

duplicons1 = []
lengths1 = []
genes1 = []
times1 = []
sims1 = []
starts1 = []
ends1 = []
duplicon2letter = {}
string = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
letter_index = 0
for i in range(len(df1)):
    if abs(df1["Synteny Block in Transformed Sequence"][i]) not in block2times2:
        block2times2[abs(df1["Synteny Block in Transformed Sequence"][i])] = 0
    if block2times1[abs(df1["Synteny Block in Transformed Sequence"][i])] > 1 or block2times2[abs(df1["Synteny Block in Transformed Sequence"][i])] > 1:
        duplicons1.append(df1["Synteny Block in Transformed Sequence"][i])
        lengths1.append((df1["End Pos in Original Sequence"][i] - df1["Start Pos in Original Sequence"][i])//1000)
        if "Genes" in df1.columns:
            genes1.append(df1["Genes"][i])
        elif "V genes" in df1.columns:
            genes1.append(df1["V genes"][i] + df1["D genes"][i] + df1["J genes"][i])
        times1.append(block2times1[abs(df1["Synteny Block in Transformed Sequence"][i])])
        starts1.append(df1["Start Pos in Original Sequence"][i])
        ends1.append(df1["End Pos in Original Sequence"][i])
        if abs(df1["Synteny Block in Transformed Sequence"][i]) not in duplicon2letter:
            duplicon2letter[abs(df1["Synteny Block in Transformed Sequence"][i])] = string[letter_index]
            letter_index += 1

if len(genes1):
    for i in genes1:
        print(f"{i}\t", end='')  # Fixed width 4
    print("#Genes")

for i in range(len(starts1)):
    if i == 0:
        print(f"{starts1[i]//1000}\t", end='')
    if i < len(starts1)-1:
        print(f"{(starts1[i+1] - ends1[i])//1000 if starts1[i+1] > ends1[i] else 0}\t", end='')
print("Distance (Kb)")

for i in lengths1:
    print(f"{i}\t", end='')  # Numbers left-justified, fixed width 4
print("Length (Kb)")

for i, dup in enumerate(duplicons1):
    letter = duplicon2letter[abs(dup)]
    if dup < 0:
        letter = '-' + letter
    if times1[i] > 1:
        print(f"{letter}\t", end='')  # Upper case, fixed width 4
    else:
        print(f"{letter.lower()}\t", end='')  # Lower case, fixed width 4
print("Genome 1")
# print()

if blocks2 != "null":
    duplicons2 = []
    lengths2 = []
    #duplicon2letter = {}
    genes2 = []
    times2 = []
    starts2= []
    ends2 = []
    sims2 = []
    # letter_index = 16
    for i in range(len(df2)):
        if abs(df2["Synteny Block in Transformed Sequence"][i]) not in block2times1:
            block2times1[abs(df2["Synteny Block in Transformed Sequence"][i])] = 0
        if block2times2[abs(df2["Synteny Block in Transformed Sequence"][i])] > 1 or block2times1[abs(df2["Synteny Block in Transformed Sequence"][i])] > 1:
            duplicons2.append(df2["Synteny Block in Transformed Sequence"][i])
            lengths2.append((df2["End Pos in Original Sequence"][i] - df2["Start Pos in Original Sequence"][i])//1000)
            if "Genes" in df2.columns:
                genes2.append(df2["Genes"][i])
            elif "V genes" in df2.columns:
                genes2.append(df2["V genes"][i] + df2["D genes"][i] + df2["J genes"][i])
            times2.append(block2times2[abs(df2["Synteny Block in Transformed Sequence"][i])])
            starts2.append(df2["Start Pos in Original Sequence"][i])
            ends2.append(df2["End Pos in Original Sequence"][i])
            if abs(df2["Synteny Block in Transformed Sequence"][i]) not in duplicon2letter:
                duplicon2letter[abs(df2["Synteny Block in Transformed Sequence"][i])] = string[letter_index]
                letter_index += 1
    for i, dup in enumerate(duplicons2):
        letter = duplicon2letter[abs(dup)]
        if dup < 0:
            letter = '-' + letter
        if times2[i] > 1:
            print(f"{letter}\t", end='')  # Upper case, fixed width 4
        else:
            print(f"{letter.lower()}\t", end='')  # Lower case, fixed width 4

    print("Genome 2")
    # print()

    for i in lengths2:
        print(f"{i}\t", end='')  # Numbers left-justified, fixed width 4
    print("Length (Kb)")

    for i in range(len(starts2)):
        if i == 0:
            print(f"{starts2[i]//1000}\t", end='')
        if i < len(starts2)-1:
            print(f"{(starts2[i+1] - ends2[i])//1000 if starts2[i+1] > ends2[i] else 0}\t", end='')
    print("Distance (Kb)")

    if len(genes2):
        for i in genes2:
            print(f"{i}\t", end='')  # Fixed width 4
        print("#Genes")