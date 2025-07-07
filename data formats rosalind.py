from Bio import Entrez
from Bio import SeqIO
#store the given id as a var
given_id = ["JX317622 JX308813 NM_001194889 JX462669 NM_001270868 JX308815 JN573266 JX393321 JX462666 NM_001081821"]
Entrez.email = "enter your email@gmail.com"

#grab the ids and relevant information in fasta format as a records variable
handle= Entrez.efetch(db="nucleotide", id=given_id, rettype="fasta")
records = list (SeqIO.parse(handle, "fasta"))

#create two lists the count list and the final shortest count
#the shortest count will be used as our output 
#the count list will contain counts of each sequence in the records list
shortest_count =[]
count_list =[]

#create a for loop for each record in the records
for i in range(0,len(records)):
    #then append the counts of each to the count list
    count_list.append(len(records[i].seq))
    

#since the for loop went through each record 
#the indices of the counts should be the same
#now we can find the shortest sequence count indicies 
ind_short = count_list.index(min(count_list))

#then since we know the indices index we can apply it to the records
#to have the records with the lowest sequence count
shortest_count = records[ind_short]

fasta_string = f">{shortest_count.description}\n{shortest_count.seq}"
print(fasta_string)

