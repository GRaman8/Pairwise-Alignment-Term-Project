import numpy as np

def get_scoring_function():
    score_map = {}                                                 # create an empty array
    print("Enter scoring function values (format: A C 2 means match A-C scores 2)")
    while True:
        line = input("Scoring function (or 'done' to finish): ")
        if line.lower() == 'done':                                 # looks for whether the user is done with entering the scoring functions inputs.
            break
        parts = line.split()
        if len(parts) != 3:                                        #checks for the input format stuff.
            print("Invalid format. Use: <char1> <char2> <score>")
            continue
        score_map[(parts[0], parts[1])] = int(parts[2])
        score_map[(parts[1], parts[0])] = int(parts[2])  # Ensure symmetry
    return score_map

def affine_gap_alignment(seq1, seq2, scoring, gap_open, gap_extend):
    m, n = len(seq1), len(seq2)
    
    M = np.zeros((m+1, n+1))  # Match/mismatch matrix
    I = np.zeros((m+1, n+1))  # Insertion matrix
    D = np.zeros((m+1, n+1))  # Deletion matrix
    
    # Initialize matrices
    for i in range(1, m+1):
        I[i][0] = gap_open + (i-1) * gap_extend
        M[i][0] = I[i][0]
    for j in range(1, n+1):
        D[0][j] = gap_open + (j-1) * gap_extend
        M[0][j] = D[0][j]
    
    # Fill matrices
    for i in range(1, m+1):
        for j in range(1, n+1):
            match_score = scoring.get((seq1[i-1], seq2[j-1]), -1)  # Default -1 if no score given
            M[i][j] = max(M[i-1][j-1] + match_score, I[i-1][j-1] + match_score, D[i-1][j-1] + match_score)
            I[i][j] = max(M[i][j-1] + gap_open, I[i][j-1] + gap_extend)
            D[i][j] = max(M[i-1][j] + gap_open, D[i-1][j] + gap_extend)
    
    # Backtrack to get alignment
    alignment1, alignment2 = "", ""
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0 and M[i][j] == M[i-1][j-1] + scoring.get((seq1[i-1], seq2[j-1]), -1):
            alignment1 = seq1[i-1] + alignment1
            alignment2 = seq2[j-1] + alignment2
            i -= 1
            j -= 1
        elif j > 0 and M[i][j] == I[i][j]:
            alignment1 = "-" + alignment1
            alignment2 = seq2[j-1] + alignment2
            j -= 1
        else:
            alignment1 = seq1[i-1] + alignment1
            alignment2 = "-" + alignment2
            i -= 1
    
    return alignment1, alignment2, M[m][n]

if __name__ == "__main__":
    seq1 = input("Enter first DNA sequence: ")
    seq2 = input("Enter second DNA sequence: ")
    
    scoring = get_scoring_function()
    gap_open = int(input("Enter gap opening penalty: "))
    gap_extend = int(input("Enter gap extension penalty: "))
    
    alignment1, alignment2, score = affine_gap_alignment(seq1, seq2, scoring, gap_open, gap_extend)
    print("\nOptimal Alignment:")
    print(alignment1)
    print(alignment2)
    print(f"Alignment Score: {score}")
