#! /usr/bin/envpython
from biology import *
def test_est_base():
    assert biology.est_base("A") == True
    assert biology.est_base("T") == True
    assert biology.est_base("Z") == False
    print("test ok")
    
def test_est_adn():
    assert biology.est_adn("ATGTCAAA") == True
    assert biology.est_adn("ATBOAATG") == False
    print("test ok")

def test_arn():
    assert biology.arn("ATTGCA")=="AUUGCA"
    assert biology.arn("ATBOAATG")== None
    print("test ok")

def test_arn_to_codons():
    assert biology.arn_to_codons("CGUUAGGGG") == ["CGU", "UAG", "GGG"]
    assert biology.arn_to_codons("CGUAAU") == ["CGU", "AAU"]
    assert biology.arn_to_codons("CGUAAUGC") == ["CGU", "AAU"]
    print("test ok")

def test_codons_to_aa():
    assert biology.codons_to_aa(["CGU", "AAU", "UAA", "GGG", "CGU"], codonaa) == ["Arginine", "Asparagine"] 
    assert biology.codons_to_aa(["UUC","GUU","UAG","UAU","UCC"],codonaa) == ["Phenylalanine","Valine"]
    print("test ok")

def test_nextIndice():
    assert nextIndice(["bonjour", "hello", "buongiorno", "ciao", "bye"], 3, ["hello", "bye"]) == 4
    assert nextIndice(["hi", "yes"], 2, ["no", "how"]) == 2
    assert nextIndice(["bonjour", "hello", "buongiorno", "ciao", "bye"], 0, ["hello", "bye"]) == 1
    print("test ok")

def test_decoupe_sequence():
    assert decoupe_sequence(["val1", "début", "val2", "val3", "end", "val4", "fin", "begin", "val5", "fin", "val6"], ["début", "begin"], ["fin", "end"]) == [    ["val2", "val3"],
    ["val5"]]
    assert decoupe_sequence(["a", "b", "start", "c", "end"], ["start", "deb"], ["end", "fin"]) == [["c"]]
    print("test ok")

def test_codons_to_seq_codantes():
    assert codons_to_seq_codantes(["CGU", "UUU", "AUG", "CGU", "AUG", "AAU", "UAA", "AUG", "GGG", "CCC",  "CGU", "UAG", "GGG"], codonaa) == [["CGU", "AUG", "AAU"],["GGG", "CCC", "CGU"]]
    print("test ok")

def test_seq_codantes_to_seq_aas():
    assert seq_codantes_to_seq_aas([["CGU", "AUG", "AAU"],["GGG", "CCC", "CGU"]], codonaa) == [["Arginine", "Methionine", "Asparagine"],["Glycine", "Proline", "Arginine"]]
    print("test ok")

def test_adn_encode_molecule():
    assert adn_encode_molecule("CGTTTTATGCGTATGAATTAAATGGGGCCCCGTTAGGGG", codonaa, [["Arginine", "Methionine", "Asparagine"],["Glycine", "Proline", "Arginine"]]) == True
    assert adn_encode_molecule("CGTTTTATGCGTATGAATTAAATGGGGCCCCGTTAGGGG", codonaa, ["Glycine", "Proline", "Arginine"]) == False
    print("test ok")


