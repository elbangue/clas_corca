
# coding: utf-8

# In[20]:


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.Restriction import*

print ("********************************************************************************************************************")
print ("Classification of C.cassiicola isolates by the in silico restriction polymorphism analysis of ga4 and caa5 sequences")
print ("********************************************************************************************************************")
print ("Instituto Biológico de Sao Paulo (IB)| 2019 | elbangue@gmail.com")
print ("")

#Entrada de datos a analizar
entrada_caa5= raw_input ("Path and filename of caa5 sequences in fasta format: ")
entrada_ga4= raw_input ("Path and filename of ga4 sequences in fasta format: ")
print entrada_caa5
print entrada_ga4
#Entradas y Validación 
lista_d_validacion_ga4 = []
lista_d_validacion_caa5 = []

#Validando ga4 vs caa5
for seq_record_ga4 in SeqIO.parse(entrada_ga4, "fasta"):
    lista_d_validacion_ga4.append(seq_record_ga4.id)
for seq_record_caa5 in SeqIO.parse(entrada_caa5, "fasta"):
    lista_d_validacion_caa5.append(seq_record_caa5.id)
lista_d_validacion_ga4.sort()
lista_d_validacion_caa5.sort()
misalida =[]
misalida2 =[]
print ("********************************************************************************************************************")

if len(lista_d_validacion_ga4)!= len(lista_d_validacion_caa5):
    print "Ga4=>", len(lista_d_validacion_ga4), "/ caa5 =>", len(lista_d_validacion_caa5), ", some isolates wont be analyzed, edit sequences in files !!!!"
    print ("")
else:
    print "Same number of ga4 and caa5 fasta sequences to analyse!! /", len (lista_d_validacion_caa5)
    print ("")
for a in lista_d_validacion_caa5:
    if a in lista_d_validacion_ga4:
        misalida2.append(a)
    else:
        misalida.append(a)
if misalida !=[]:
    print "Problems with some ID!, edit please! /", misalida
    print ("")
    
#FIN VALIDACIONES
        
narrador=Seq("")

PhL1_ga4=RestrictionBatch([BsrI, BstZ17I])
PhL2_ga4=RestrictionBatch([RsaI,TaqI])
PhL3_ga4=RestrictionBatch([AluI,BstZ17I, TaqI])
PhL4_ga4=RestrictionBatch([AluI, TaqI])
PhL5_ga4=RestrictionBatch([BstZ17I, TaqI])
PhL6_ga4=RestrictionBatch([BsrI, BstZ17I, TaqI])
PhL7_ga4=RestrictionBatch([BsrI, BstZ17I, TaqI])
PhL8_ga4=RestrictionBatch([BstZ17I, TaqI])

PhL1_ga4_group=[]
PhL2_ga4_group=[]
PhL3_ga4_group=[]
PhL4_ga4_group=[]
PhL5_ga4_group=[]
PhL6_ga4_group=[]
PhL7_ga4_group=[]
PhL8_ga4_group=[]

for seq_record_ga4 in SeqIO.parse(entrada_ga4, "fasta"):
    ga4_inicio = seq_record_ga4.seq.find("CCCAACGTCATGACTC")
    ga4_final = seq_record_ga4.seq.find("TCAGCCATCAACC")
   
    if (ga4_inicio > ga4_final) or (ga4_inicio == ga4_final):
        seq_record_ga4.seq = seq_record_ga4.seq.reverse_complement()
        ga4_final = seq_record_ga4.seq.find("GAGTCATGACGTTGGG")
        ga4_inicio = seq_record_ga4.seq.find("GGTTGATGGCTGA")
        
        
       
    
    rb=RestrictionBatch([RsaI, BsrI, TaqI, BstZ17I, AluI])
    moment_isolate = RestrictionBatch([])
    busqueda=rb.search (seq_record_ga4.seq)
    if busqueda!=[]:
        for x in busqueda:
            if busqueda.get(x)!= []: 
                for y in busqueda.get(x):
                    if y >= ga4_inicio and y <= ga4_final:
                        moment_isolate.add (x)
                        
 
    
  
    if moment_isolate == PhL1_ga4:
        PhL1_ga4_group.append (seq_record_ga4)
    if moment_isolate == PhL2_ga4:
        PhL2_ga4_group.append (seq_record_ga4)
    if moment_isolate == PhL3_ga4:
        PhL3_ga4_group.append (seq_record_ga4)
    if moment_isolate == PhL4_ga4:
        PhL4_ga4_group.append (seq_record_ga4)
    if moment_isolate == PhL5_ga4:
        PhL5_ga4_group.append (seq_record_ga4)
    if moment_isolate == PhL6_ga4:
        PhL6_ga4_group.append (seq_record_ga4)
    if moment_isolate == PhL7_ga4:
        PhL7_ga4_group.append (seq_record_ga4)
    if moment_isolate == PhL8_ga4:
        PhL8_ga4_group.append (seq_record_ga4)
    
        
PhL1_caa5=RestrictionBatch([])
PhL2_caa5=RestrictionBatch([BstEII, SalI])
PhL3_caa5=RestrictionBatch([])
PhL4_caa5=RestrictionBatch([NaeI])
PhL5_caa5=RestrictionBatch([])
PhL6_caa5=RestrictionBatch([BstZI])
PhL7_caa5=RestrictionBatch([])
PhL8_caa5=RestrictionBatch([BstEII])

PhL1_caa5_group=[]
PhL2_caa5_group=[]
PhL3_caa5_group=[]
PhL4_caa5_group=[]
PhL5_caa5_group=[]
PhL6_caa5_group=[]
PhL7_caa5_group=[]
PhL8_caa5_group=[]

PhL1=[]
PhL2=[]
PhL3=[]
PhL4=[]
PhL5=[]
PhL6=[]
PhL7=[]
PhL8=[]
                               
for seq_record_caa5 in SeqIO.parse(entrada_caa5, "fasta"):   
     
   
    rb=RestrictionBatch([BstEII, BstZI, NaeI, SalI])
    moment_isolate = RestrictionBatch([])
    busqueda=rb.search (seq_record_caa5.seq)
    if busqueda!=[]:
        for x in busqueda:
            if busqueda.get(x) != []:
                for y in busqueda.get(x):
               # print busqueda.get(x)
                    if y!=0:
                        moment_isolate.add (x)
                        
  
   
    if moment_isolate == PhL1_caa5 or moment_isolate == PhL8_caa5:#este agrego se debe a algunas exepciones observadas
        PhL1_caa5_group.append (seq_record_caa5)
    if moment_isolate == PhL2_caa5:
        PhL2_caa5_group.append (seq_record_caa5)
    if moment_isolate == PhL3_caa5:
        PhL3_caa5_group.append (seq_record_caa5)
    if moment_isolate == PhL4_caa5:
        PhL4_caa5_group.append (seq_record_caa5)
    if moment_isolate == PhL5_caa5:
        PhL5_caa5_group.append (seq_record_caa5)
    if moment_isolate == PhL6_caa5:
        PhL6_caa5_group.append (seq_record_caa5)
    if moment_isolate == PhL7_caa5:
        PhL7_caa5_group.append (seq_record_caa5)
    if moment_isolate == PhL8_caa5:
        PhL8_caa5_group.append (seq_record_caa5)

print ("***********************PhL1***********************")

for ga4 in PhL1_ga4_group:
    for caa5 in PhL1_caa5_group:
        if ga4.id == caa5.id:
            PhL1.append(ga4.id)
            break
print("PhL1 #/"), len(PhL1)
print "PhL1 isolates ID ==>", PhL1

print ("***********************PhL2***********************")            
for ga4 in PhL2_ga4_group:
    for caa5 in PhL2_caa5_group:
        if ga4.id == caa5.id:
            PhL2.append(ga4.id)
            break
print("PhL2 #/"), len(PhL2)
print "PhL2 isolates ID ==> ", PhL2

print ("***********************PhL3***********************")            
for ga4 in PhL3_ga4_group:
    for caa5 in PhL3_caa5_group:
        if ga4.id == caa5.id:
            PhL3.append(ga4.id)
            break
print("PhL3 #/"), len(PhL3)
print "PhL3 isolates ID ==> ", PhL3

print ("***********************PhL4***********************")            
for ga4 in PhL4_ga4_group:
    for caa5 in PhL4_caa5_group:
        if ga4.id == caa5.id:
            PhL4.append(ga4.id)
            break
print("PhL4 #/"), len(PhL4)
print "PhL4 isolates ID ==>", PhL4

print ("***********************PhL5***********************")            
for ga4 in PhL5_ga4_group:
    for caa5 in PhL5_caa5_group:
        if ga4.id == caa5.id:
            PhL5.append(ga4.id)
            break
print("PhL5 #/"), len(PhL5)
print "PhL5 isolates ID ==> ", PhL5

print ("***********************PhL6***********************")            
for ga4 in PhL6_ga4_group:
    for caa5 in PhL6_caa5_group:
        if ga4.id == caa5.id:
            PhL6.append(ga4.id)
            break
print("PhL6 #/"), len(PhL6)
print "PhL6 isolates ID ==>", PhL6

print ("***********************PhL7***********************")            
for ga4 in PhL7_ga4_group:
    for caa5 in PhL7_caa5_group:
        if ga4.id == caa5.id:
            PhL7.append(ga4.id)
            break
print("PhL7 #/"), len(PhL7)
print "PhL7 isolates ID ==>", PhL7

print ("***********************PhL8***********************")            
for ga4 in PhL8_ga4_group:
    for caa5 in PhL8_caa5_group:
        if ga4.id == caa5.id:
            PhL8.append(ga4.id)
            break
print("PhL8 #/"), len(PhL8)
print "PhL8 isolates ID ==>", PhL8

lista_de_no_clas =[]
clas_isolates = PhL1 + PhL2 +PhL3+PhL4 + PhL5 +PhL6+PhL7+PhL8
print ("All classified isolates /"), clas_isolates
print "Total of analyzed isolates /", len(misalida2)
print ("Number of classified isolates /"), len (clas_isolates)
print ("Number of unclassified isolates /"), len(misalida2) -len(clas_isolates)
for a in misalida2:
    if not a in clas_isolates:
        lista_de_no_clas.append (a)
if lista_de_no_clas !=[]:
    print ("Some isolates were impossible to classify :( !!!!")
    print lista_de_no_clas, misalida
else: 
    print ("All isolates were succefully classified!:) !!!!")

