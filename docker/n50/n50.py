#!/usr/bin/env python3

"""
Author:     Jonathan Stites
Date:       03 April 2015
Company:    Dovetail Genomics


Takes as input a fasta file or gzipped fasta file and outputs statistics
in json format about the contiguity of the assembly contigs and assembly 
scaffolds.

Options:
-a, --assembly       fasta or gzipped-fasta
-o, --output         output file
-m, --minlen         minimum sequence length to count towards n50, etc.
-k                   discount for kmer length of k (1 gives default every-base-counts)

Usage:
assembly_statistics.py -a hg19.fa -o hg19_stats.json

Notes:
Detects gzip fasta by extension. For contig statistics,
breaks scaffolds at runs of 10 or more consecutive Ns. By default,
statistics are based on contigs/scaffolds with >= [min_length] bases.
"""


import argh
import json
import gzip
from collections import defaultdict


def main(assembly=None, output="/dev/stdout", minlen=1000, genome=0, brief=False, k=1):
    stats = {}
    gc, at, n, lengths, con_gc, con_at, con_n, con_lengths = read_assembly(assembly, minlen)
    if (k > 1):
        lengths = [(length + 1 - k) for length in lengths]
        con_lengths = [(length + 1 - k) for length in con_lengths]
    stats["contigs"] = get_stats(con_lengths, con_gc, con_at, con_n, minlen, genome)
    if (lengths != con_lengths and n != con_n):
        stats["scaffolds"] = get_stats(lengths, gc, at, n, minlen, genome)
    else:
        stats["scaffolds"] = stats["contigs"]

    with open(output, 'w') as out:
        if brief:
            for key in stats.keys():
                print(assembly, key,
                      stats[key]["# contigs"], stats[key]["Total length"], 
                      stats[key]["N50"], stats[key]["L50"],
                      file=out)
        else:
            json.dump(stats, out, indent=4, sort_keys=True)

def read_assembly(assembly, minlen):
    gc, at, n = 0, 0, 0
    con_gc, con_at, con_n = 0, 0, 0
    lengths = []
    con_lengths = []
    for header, seq in read_fasta(assembly):
        gc, at, n = update_gc(seq, gc, at, n, minlen)
        lengths.append(len(seq))
        for contig in break_scaffold(seq):
            con_gc, con_at, con_n = update_gc(contig, con_gc, con_at, con_n, minlen, scaffolds=False)
            con_lengths.append(len(contig))

    return gc, at, n, lengths, con_gc, con_at, con_n, con_lengths

def get_stats(lengths, gc, at, n, minlen, genome):
    stats = {}
    stats["# contigs (>= 0 bp)"] = len(lengths)
    stats["# contigs (>= 1000 bp)"] = sum(1 for i in lengths if i >= 1000)
    stats["Total length (>= 0 bp)"] = sum(lengths)
    stats["Total length (>= 1000 bp)"] = sum(i for i in lengths if i >= 1000)
    stats["# contigs"] = sum(1 for i in lengths if i >= minlen)
    stats["Largest contig"] = max(lengths)
    stats["Total length"] = sum(i for i in lengths if i >= minlen)
    stats["Estimated genome size"] = genome
    try:
        stats["GC (%)"] = round(gc*100/(gc+at), 2)
    except ZeroDivisionError:
        stats["GC (%)"] = 0.0
    n50, l50 = nx(lengths, .5, minlen, genome)
    n90, l90 = nx(lengths, .90, minlen, genome)
    g = 'G' if genome else ''
    stats[f"N{g}50"] = n50
    stats[f"N{g}90"] = n90
    stats[f"L{g}50"] = l50
    stats[f"L{g}90"] = l90
    try:
        stats["# N's per 100 kbp"] = round(n*100000/ sum([i for i in lengths if i >= minlen]), 2)
    except ZeroDivisionError:
        stats["# N's per 100 kbp"] = 0.0
    return stats

def read_fasta(assembly):
    if assembly.endswith(".gz"):
        read_mode = "rb"
        open_cmd = gzip.open
    else:
        read_mode = "r"
        open_cmd = open
        
    with open_cmd(assembly, mode=read_mode) as assembly_handle:
        header = None
        for line in assembly_handle:
            if read_mode == "rb":
                line = line.decode("utf-8")
            if line.startswith(">"):
                if header:
                    yield header, ''.join(seq)
                header = line.rstrip()
                seq = []
            else:
                seq.append(line.rstrip())
        yield header, ''.join(seq)


def update_gc(seq, gc, at, n, minlen, scaffolds=True):
    gc_bases = set(["G", "g", "C", "c"])
    at_bases = set(["A", "a", "T", "t"])
    n_bases = set(["N", "n"])
    if scaffolds and len(seq) < minlen:
        return gc, at, n

    counts = get_counts(seq)
    for base in counts.keys():
        if base in gc_bases:
            gc += counts[base]
        elif base in at_bases:
            at += counts[base]
        else:
            n += counts[base]
    return gc, at, n

def get_counts(seq):
    counts = defaultdict(int)
    for base in seq:
        counts[base] += 1
    return counts

def nx(lengths, fraction, minlen, genome=0):
    sorted_lengths = sorted(lengths, reverse=True)
    mid = (sum([i for i in sorted_lengths if i >= minlen]
               ) if not genome else genome) *fraction
    total = 0
    for i, length in enumerate(sorted_lengths):
        total += length
        if total >= mid:
            return length, i + 1
    return 0, 0

def break_scaffold(seq):
    i = 0
    cur_contig_start = 0

    while (i < len(seq)) and (seq.find("N", i) != -1):
        start = seq.find("N", i)
        end = start + 1
        while (end != len(seq)) and (seq[end] == "N"):
            end +=1
        
        i = end + 1
        if (end - start) >= 10:
            if cur_contig_start != start:
                yield seq[cur_contig_start: start]
            cur_contig_start = end
    if len(seq[cur_contig_start: ]) != 0:
        yield seq[cur_contig_start: ]


class InvalidBase(Exception):
    pass

if __name__ == "__main__":
    argh.dispatch_command(main)

