# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 15:03:53 2022

@author: Dell
"""
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import xlwt
from xlwt import Workbook
  
def Print_Count(records): #print the percentage content and average content for A,C,T,G,CG for a list of SeqRecord
    consensus_seq=""
    A=[]
    T=[]
    C=[]
    G=[]
    CG=[]
    
    for i in range(len(records)):               #loop on each record
        A_count = records[i].seq.count('A')     #Count A in one seq
        C_count = records[i].seq.count('C')     #Count C in one seq
        G_count = records[i].seq.count('G')     #Count G in one seq
        T_count = records[i].seq.count('T')     #Count T in one seq
        
        #calculate percentage content
        A_percent=(float(A_count) / len(records[i].seq)) * 100
        C_percent=(float(C_count) / len(records[i].seq)) * 100
        T_percent=(float(T_count) / len(records[i].seq)) * 100
        G_percent=(float(G_count) / len(records[i].seq)) * 100
        CG_percent=((float(G_count)+float(C_count)) / len(records[i].seq)) * 100
        
        #store the percentage content in a list
        A.append(A_percent)
        C.append(C_percent)
        T.append(T_percent)
        G.append(G_percent)
        CG.append(CG_percent)
        
    #compute the percentage average of content
    A_avg=sum(A)/len(A)
    C_avg=sum(C)/len(C)
    T_avg=sum(T)/len(T)
    G_avg=sum(G)/len(G)
    CG_avg=sum(CG)/len(CG)    
   
    print("")
    print("Percentage Average of A content: ",A_avg)
    print("Percentage Average of C content: ",C_avg)
    print("Percentage Average of T content: ",T_avg)
    print("Percentage Average of G content: ",G_avg)
    print("Percentage Average of CG content: ",CG_avg)
    print("")


def Get_Consensus(seqs): #creat a consensus sequence from the given list of SeqRecord and save it in a fasta file
    
    consensus_seq=""
    
    for i in range(0,len(seqs[0].seq)): #loop on each index of the seqs
        char=""                      #string to hold the char at index i for all seqs
        for j in range(0,len(seqs)): #loop on seqs
            char+=seqs[j].seq[i]     #add the char at i of seq[j]
          
        res=Counter(char) #create an iterable hashtable
        consensus_seq +=max(res,key=res.get) #add to consesnsus_seq the most repeated char in char
   
    record = SeqRecord(        #creat a seq record containing the consensus seq
    Seq(consensus_seq),
    id="1",
    name="consensus",
    description="consensus sequense of 10 delta sequences",
    )
    SeqIO.write(record, "consensus.fasta", "fasta") #save the consensus seq in a fasta file(consensus.fasta)

def compare_column(f,seqs):  #return true if characters at index f in the list seqs are similar and false otherwise
        char=""                      #string to hold the char at index f for all seqs
        for j in range(0,len(seqs)): #loop on seqs
            char+=seqs[j].seq[f]     #add the char at f of seq[j]
        if(char.count(char[0])==len(char)): #check if all characters of char are similar to the first character
            return True
        else:
            return False
def Get_Dissimilar_Region(seqs,consensus): #print the indexes of dissimilar regions between seqs and consensus by comparing similar regions in seqs with tha same region in consensus seq 
    i=0 #pointer1
    k=0 #pointer2
    similar_regions=[]  #indexes of start and end of similar regions in case seqs
    dissimilar_regions=[]  #indexes of dissimilar regions between case seqs and consensus seq
    
    #obtain similair regions in case sequences
    while(i<len(seqs[0].seq)):  
        k=i+1
        if(compare_column(i,seqs)): #start of a probable similar region
            similar_column=True
            while(similar_column and k<len(seqs[0].seq)): #loop untill reaching a column where it is not similar
                if(compare_column(k,seqs)):
                    k+=1
                else:
                    similar_column=False
        if(k>i+1): #similar region is found
            similar_regions.append(i)
            similar_regions.append(k-1)
        i=k+1 
   
    
    # #make lenght of consensus seq equal to case seqs by adding "E" at the end only if lenght of consensus seq<case seqs
    # if(len(consensus[0].seq)<len(seqs[0].seq)):
    #     consensus[0].seq+="E"*(len(seqs[0].seq)-len(consensus))
    
    #find dissimilar regions by comparing similar regions defined in similar_regions with tha same region in consensus seq 
    n=0
    row=1
    #create an excel sheet
    wb = Workbook()
    sheet1 = wb.add_sheet('Sheet 1')
    sheet1.write(0,0, 'index')
    sheet1.write(0,1, 'case sequence')
    sheet1.write(0,2, 'refrence sequence') 
    print("")
    print("Dissimilar regions at idexes:")
   #compare between obtained similar regions in case sequense and the consensus sequence
    while (n<len(similar_regions)):
       
        if(consensus[0].seq[similar_regions[n]:similar_regions[n+1]+1])!=(seqs[0].seq[similar_regions[n]:similar_regions[n+1]+1]):
            
            dissimilar_regions.append(similar_regions[n])
            dissimilar_regions.append(similar_regions[n+1])
            print(f"{similar_regions[n]} to {similar_regions[n+1]}")
            sheet1.write(row,0, f"{similar_regions[n]} to {similar_regions[n+1]}")
            sheet1.write(row,1, f"{seqs[0].seq[similar_regions[n]:similar_regions[n+1]+1]}")
            sheet1.write(row,2, f"{consensus[0].seq[similar_regions[n]:similar_regions[n+1]+1]}")
        row+=1
        n=n+2
    wb.save('Dissimilar regions.xls')
    
#Main
refrence_seqs=list(SeqIO.parse("delta_align.fasta", "fasta")) 
Get_Consensus(refrence_seqs) #get consensus sequence

omicron = list(SeqIO.parse("Omicron10.fasta", "fasta"))
print("case sequences")
Print_Count(omicron) 

consensus= list(SeqIO.parse("consensus.fasta", "fasta"))  
print("refrence sequence")
Print_Count(consensus)
omicron_aligned_with_consensus=list(SeqIO.parse("aligned_omicron_with_consensus.fasta", "fasta")) 
aligned_consensus=list(SeqIO.parse("aligned_consensus.fasta","fasta"))
Get_Dissimilar_Region(omicron_aligned_with_consensus,aligned_consensus)