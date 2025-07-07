import numpy as np

#first we will have to define the algorithm as a function
#of the two sequences given and given scores 
def needleman_wunsch(seq1, seq2, match_score=1, mismatch_score=-1, gap_penalty=-1):
    
    #now lets define the length of the sequences
    n = len(seq1)
    m = len(seq2)
    
    #Initialize the scoring matrix
    score_matrix = np.zeros((n+1, m+1), dtype=int)
    
    #Initialize the first row and column
    for i in range(n+1):
        score_matrix[i][0] = gap_penalty * i
    
    for j in range(m+1):
        score_matrix[0][j] = gap_penalty * j
    
    #fill the scoring matrix
    for i in range(1, n+1):
        for j in range(1, m+1):
            match= score_matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score)
            delete = score_matrix[i-1][j] + gap_penalty
            insert = score_matrix[i][j-1] + gap_penalty
            score_matrix[i][j] = max(match, delete, insert)
            
    aligned_seq1 = ""
    aligned_seq2 = ""
    i = n
    j = m
    
    while i > 0 or j > 0:
        if i > 0 and j > 0 and score_matrix[i][j] == score_matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score):
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
        elif i > 0 and score_matrix[i][j] == score_matrix[i-1][j] + gap_penalty:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        else:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1
    
    return aligned_seq1, aligned_seq2, score_matrix[n][m]

#now lets try an example
seq1= "ATTACGCTA"
seq2= "GTTACGCTT"   

aligned_seq1, aligned_seq2, score = needleman_wunsch(seq1, seq2)

print(aligned_seq1)