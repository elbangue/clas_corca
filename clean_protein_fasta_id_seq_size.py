#!/usr/bin/env python
# coding: utf-8

# In[9]:


from Bio import SeqIO
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC



print ("********************************************************************************************************************")
print ("Clean Bulk of Protein Sequences in fasta format considering ID, AA sequence and User arbitrary size")
print ("********************************************************************************************************************")
print ("University of Saarland| 2020 | elbangue@gmail.com")
print ("")

#Entrada de datos a analizar
entrada_fasta= input("Path and filename of file in fasta format: ")
size_cut=100
size_cut = int(input("Protein abitrary AA inferior limit (100aa default): "))
print ("Processing.....", entrada_fasta)

mylist = []
# Read entries, fasta: sequence
for ident in SeqIO.parse(entrada_fasta, "fasta"):
    floatlist=[]
    
    #Do list considering identical fasta ID and sequence
    for ident2 in SeqIO.parse(entrada_fasta, "fasta"):
        if ident.id == ident2.id or ident.seq ==ident2.seq: #if only want selected ID or sequence then delete one
            floatlist.append (ident2) 
            
    #Initial values of accesories variables        
    y=0 
    d=0
    
    #New fasta list without redundant identical ID and sequences
    if len(floatlist) > 1:
        for x in floatlist:
            if len(x.seq)>=y:
                y = len(x.seq)
                new =d
            d=d+1    
        if mylist != []:
            condi = True
            for coge in mylist:
                if coge.id == floatlist[new].id:
                    condi = False
            if condi == True:
                if len(floatlist[new].seq)>size_cut:
                    mylist.append(floatlist[new])
        
    else:
        if len(ident.seq)>size_cut:
            mylist.append(ident)
            
y=0 #Initial value of accesorie variable

for char in entrada_fasta:
    y=y+1
    if char == "." :
        d = y-1
        
#Write cleaned fasta to file 
SeqIO.write(mylist, entrada_fasta[0:d] + '_proccess_id&seq&size' + str(size_cut) +'aa.fasta', 'fasta')   
print ("Copy file to:", entrada_fasta[0:d] + '__proccess_id&seq&size' + str(size_cut) +'aa.fasta')

