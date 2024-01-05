#! /usr/bin/envpython

def est_base(x):
    if x in "ATCG":
        return True
    else :
        return False

def est_adn(x):
    i = 0
    while i <= len(x)-1:
        if x[i] in "ATCG":
            i += 1
        else :
            return False
    return True

def arn(x):
    if est_adn(x):
        arn = x.replace("T","U")
        return arn
    else :
        return None 

def arn_to_codons(x):
    tab = []
    for i in range(0, len(x), 3) :
        codon = x[i:i+3]
        tab.append(codon)
    if len(x) % 3 == 0:
        return tab
    else :
        x = len(tab)
        del(tab[x-1])
        return tab

from json import *
def load_dico_codons_aa(x):
    fichier = open (x,"r")
    strjson = fichier.read()
    fichier.close()
    codonaa = loads(strjson)
    return codonaa

from itertools import product

def codons_stop(x):
    codons = []
    codonsstop = []
    t = list(x)
    L = ["A", "U", "G", "C"]
    for i in product(L, repeat=3):
        codons.append(i)
    k = 0
    while k < len(codons):
        s = ""
        j = 0
        while j < len(codons[k]):
            s += codons[k][j]
            j += 1
        if s not in t:
            codonsstop.append(s)
        k += 1
    return codonsstop

def codons_to_aa(codons, codonaa):
    tab = []
    i = 0
    while i < len(codons) and codons[i] in codonaa :
            aa = codonaa[codons[i]]
            tab.append(aa)
            i += 1
    return tab

def nextIndice(tab, ind, elements):
    i = ind
    while i < len(tab) :
        if tab[i] in elements:
            return i
        i += 1
    return len(tab)

def decoupe_sequence(seq, start, stop):
    tab = []
    i = 0
    while i < len(seq):
        y = []
        if seq[i] in start :
            i += 1
            while seq[i] not in stop :
                y.append(seq[i])
                i += 1
            tab.append(y)
        else :
            i += 1
    return tab

def codons_to_seq_codantes(seq, codonaa):
    codonsstop = codons_stop(codonaa)
    tab = decoupe_sequence(seq, ["AUG"], codonsstop)
    return tab

def seq_codantes_to_seq_aas(seq, codonaa):
    tab = []
    i = 0
    while i < len(seq):
        k = seq[i]
        j = 0
        y = []
        while j < len(k):
            x = codonaa[k[j]]
            y.append(x)
            j += 1
        tab.append(y)
        i += 1
    return tab

def adn_encode_molecule(brin, codonaa, seq):
    if est_adn(brin):
        ARN = arn(brin)
        codon = arn_to_codons(ARN)
        sequence = codons_to_seq_codantes(codon, codonaa)
        aa = seq_codantes_to_seq_aas(sequence, codonaa)
        if aa == seq :
            return True
        else :
            return False
