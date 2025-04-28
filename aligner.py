# Version-3
print("Version-3")
print("DNA Sequence Alignment Tool with Affine Gap Penalty")
import numpy as np
import re

def get_scoring_function():                                                                  # Manual Inputs for the Scoring function
    """Get custom scoring function from user input"""
    score_map = {}
    print("Enter scoring function values (format: A C 2 means match A-C scores 2). Type 'done' to finish.")
    while True:
        line = input("Scoring function: ")
        if line.lower() == 'done':
            break
        parts = line.strip().split()
        if len(parts) != 3 or not parts[2].lstrip('-').isdigit():
            print("Invalid format. Use: [base1] [base2] [score]")
            continue
        score_map[(parts[0], parts[1])] = int(parts[2])
        score_map[(parts[1], parts[0])] = int(parts[2])  # Ensure symmetry
    return score_map

def is_valid_dna(seq):                                                                        # Input validation function
    """Check if sequence contains only valid DNA characters"""
    return re.fullmatch(r"[ACGT]+", seq.upper()) is not None

def print_matrices(M, I, D, seq1, seq2):                                                      # Visualization of the Work that has been done in the form of matrices 
    """Print dynamic programming matrices for debugging"""
    print("\nMatch/Mismatch Matrix (M):")
    print_matrix(M, seq1, seq2)
    print("\nInsertion Matrix (I):")
    print_matrix(I, seq1, seq2)
    print("\nDeletion Matrix (D):")
    print_matrix(D, seq1, seq2)

def print_matrix(matrix, seq1, seq2):
    """Helper function to print a single matrix with sequence labels"""                       # Actual Visualization work for the matrices.     
    # Print header
    print("    ", end="")
    print("   ", end="")
    for j in range(len(seq2)):
        print(f"{seq2[j]:4}", end="")
    print()
    
    # Print rows with labels
    for i in range(matrix.shape[0]):
        if i == 0:
            print("  ", end="")
        else:
            print(f"{seq1[i-1]} ", end="")
        
        for j in range(matrix.shape[1]):
            if matrix[i][j] == float('-inf'):
                print("-∞  ", end="")
            else:
                print(f"{int(matrix[i][j]):3} ", end="")
        print()

def visualize_alignment(alignment1, alignment2):
    """Visualize the alignment with match/mismatch indicators"""
    match_line = ""
    for i in range(len(alignment1)):
        if alignment1[i] == alignment2[i]:
            match_line += "|"
        elif alignment1[i] == "-" or alignment2[i] == "-":
            match_line += " "
        else:
            match_line += "x"
    
    print(alignment1)
    print(match_line)
    print(alignment2)

def affine_gap_alignment(seq1, seq2, scoring, gap_open, gap_extend, debug=False):
    """
    Perform pairwise sequence alignment with affine gap penalty
    
    Args:
        seq1, seq2: DNA sequences to align
        scoring: Dictionary with scoring function values
        gap_open: Penalty for opening a gap
        gap_extend: Penalty for extending a gap
        debug: If True, print matrices for debugging
        
    Returns:
        Tuple of (alignment1, alignment2, score)
    """
    # Check for empty sequences
    if not seq1 or not seq2:
        print("Error: Sequences cannot be empty")
        return "", "", 0
        
    m, n = len(seq1), len(seq2)

    # Initialize matrices
    M = np.full((m+1, n+1), float('-inf'))  # Match/mismatch matrix
    I = np.full((m+1, n+1), float('-inf'))  # Insertion matrix (gap in seq1)
    D = np.full((m+1, n+1), float('-inf'))  # Deletion matrix (gap in seq2)
    
    # Initialize traceback matrices
    trace_M = np.full((m+1, n+1), '', dtype=object)
    trace_I = np.full((m+1, n+1), '', dtype=object)
    trace_D = np.full((m+1, n+1), '', dtype=object)

    # Initialize first cell
    M[0][0] = 0
    I[0][0] = float('-inf')
    D[0][0] = float('-inf')

    # Initialize first column
    for i in range(1, m+1):
        M[i][0] = float('-inf')
        I[i][0] = float('-inf')
        D[i][0] = gap_open + (i-1) * gap_extend
        trace_D[i][0] = 'D'  # Gap in seq2
    
    # Initialize first row
    for j in range(1, n+1):
        M[j][0] = float('-inf')
        I[0][j] = gap_open + (j-1) * gap_extend
        D[0][j] = float('-inf')
        trace_I[0][j] = 'I'  # Gap in seq1

    # Fill matrices
    for i in range(1, m+1):
        for j in range(1, n+1):
            # Match/Mismatch matrix
            match_score = scoring.get((seq1[i-1], seq2[j-1]), -1)
            
            # Calculate scores for M coming from each matrix
            from_M = M[i-1][j-1] + match_score
            from_I = I[i-1][j-1] + match_score
            from_D = D[i-1][j-1] + match_score
            
            # Find max and update traceback
            if from_M >= from_I and from_M >= from_D:
                M[i][j] = from_M
                trace_M[i][j] = 'M'
            elif from_I >= from_M and from_I >= from_D:
                M[i][j] = from_I
                trace_M[i][j] = 'I'
            else:
                M[i][j] = from_D
                trace_M[i][j] = 'D'
            
            # Insertion matrix (gap in seq1)
            from_M_to_I = M[i][j-1] + gap_open
            from_I_to_I = I[i][j-1] + gap_extend
            
            if from_M_to_I >= from_I_to_I:
                I[i][j] = from_M_to_I
                trace_I[i][j] = 'M'
            else:
                I[i][j] = from_I_to_I
                trace_I[i][j] = 'I'
            
            # Deletion matrix (gap in seq2)
            from_M_to_D = M[i-1][j] + gap_open
            from_D_to_D = D[i-1][j] + gap_extend
            
            if from_M_to_D >= from_D_to_D:
                D[i][j] = from_M_to_D
                trace_D[i][j] = 'M'
            else:
                D[i][j] = from_D_to_D
                trace_D[i][j] = 'D'

    if debug:
        print_matrices(M, I, D, seq1, seq2)

    # Find the maximum score and matrix to start traceback from
    max_score = max(M[m][n], I[m][n], D[m][n])
    if M[m][n] == max_score:
        current_matrix = 'M'
    elif I[m][n] == max_score:
        current_matrix = 'I'
    else:
        current_matrix = 'D'

    # Traceback to reconstruct alignment
    alignment1, alignment2 = "", ""
    i, j = m, n

    while i > 0 or j > 0:
        if current_matrix == 'M':
            if i > 0 and j > 0:
                if trace_M[i][j] == 'M':
                    next_matrix = 'M'
                elif trace_M[i][j] == 'I':
                    next_matrix = 'I'
                else:  # 'D'
                    next_matrix = 'D'
                    
                alignment1 = seq1[i-1] + alignment1
                alignment2 = seq2[j-1] + alignment2
                i -= 1
                j -= 1
                current_matrix = next_matrix
            else:
                # Should not happen with proper initialization
                break
                
        elif current_matrix == 'I':
            if j > 0:
                if trace_I[i][j] == 'M':
                    next_matrix = 'M'
                else:  # 'I'
                    next_matrix = 'I'
                    
                alignment1 = '-' + alignment1
                alignment2 = seq2[j-1] + alignment2
                j -= 1
                current_matrix = next_matrix
            else:
                # Should not happen with proper initialization
                break
                
        elif current_matrix == 'D':
            if i > 0:
                if trace_D[i][j] == 'M':
                    next_matrix = 'M'
                else:  # 'D'
                    next_matrix = 'D'
                    
                alignment1 = seq1[i-1] + alignment1
                alignment2 = '-' + alignment2
                i -= 1
                current_matrix = next_matrix
            else:
                # Should not happen with proper initialization
                break

    return alignment1, alignment2, max_score

def calculate_alignment_score(alignment1, alignment2, scoring, gap_open, gap_extend):
    """Calculate and verify the score of an alignment"""
    score = 0
    gap_count1 = 0  # Count consecutive gaps in seq1
    gap_count2 = 0  # Count consecutive gaps in seq2
    
    for i in range(len(alignment1)):
        if alignment1[i] == '-':
            if gap_count1 == 0:  # Opening a new gap
                score += gap_open
            else:  # Extending existing gap
                score += gap_extend
            gap_count1 += 1
            gap_count2 = 0  # Reset other gap counter
        elif alignment2[i] == '-':
            if gap_count2 == 0:  # Opening a new gap
                score += gap_open
            else:  # Extending existing gap
                score += gap_extend
            gap_count2 += 1
            gap_count1 = 0  # Reset other gap counter
        else:  # Match or mismatch
            score += scoring.get((alignment1[i], alignment2[i]), -1)
            gap_count1 = 0
            gap_count2 = 0
            
    return score

def run_test_cases(scoring, gap_open, gap_extend):
    """Run predefined test cases to verify algorithm correctness"""
    test_cases = [
        # (seq1, seq2, expected_score)
        ("ACGT", "AGT", None),  # Simple case with one gap
        ("AGTC", "AGC", None),  # Simple case with one gap
        ("ACGTACGT", "ACGACGT", None),  # Case with one gap
        ("ACGTGTCAGT", "ACGTCAGT", None),  # More complex
        ("AGCTAGCT", "AGCT", None)  # Multiple gaps
    ]
    
    print("\nRunning test cases...")
    for i, (seq1, seq2, _) in enumerate(test_cases):
        print(f"\nTest Case {i+1}: {seq1} vs {seq2}")
        alignment1, alignment2, score = affine_gap_alignment(seq1, seq2, scoring, gap_open, gap_extend)
        
        print("Alignment:")
        visualize_alignment(alignment1, alignment2)
        print(f"Score: {score}")
        
        # Verify score by manual calculation
        verified_score = calculate_alignment_score(alignment1, alignment2, scoring, gap_open, gap_extend)
        if score == verified_score:
            print(f" Score verification: algorithm score {score} matches calculated score {verified_score}")
        else:
            print(f"Score mismatch: algorithm score {score} differs from calculated score {verified_score}")

if __name__ == "__main__":
    print("Pairwise Sequence Alignment Tool (with Affine Gap Penalty)")
    print("--------------------------------------------------------")

    # Get input sequences
    seq1 = input("Enter first DNA sequence: ").upper()
    while not seq1 or not is_valid_dna(seq1):
        if not seq1:
            print("Error: Sequence cannot be empty")
        else:
            print("Invalid sequence. Only A, C, G, T are allowed.")
        seq1 = input("Enter first DNA sequence: ").upper()

    seq2 = input("Enter second DNA sequence: ").upper()
    while not seq2 or not is_valid_dna(seq2):
        if not seq2:
            print("Error: Sequence cannot be empty")
        else:
            print("Invalid sequence. Only A, C, G, T are allowed.")
        seq2 = input("Enter second DNA sequence: ").upper()

    # Define default scoring
    default_scoring = {
        ('A', 'A'): 2, ('C', 'C'): 2, ('G', 'G'): 2, ('T', 'T'): 2,
        ('A', 'C'): -1, ('A', 'G'): -1, ('A', 'T'): -1,
        ('C', 'G'): -1, ('C', 'T'): -1, ('G', 'T'): -1
    }
    
    # Show default scoring to user
    print("\nDefault scoring matrix:")
    print("  | A  | C  | G  | T  |")
    print("--------------------- ")
    for base1 in "ACGT":
        row = f"{base1} |"
        for base2 in "ACGT":
            row += f" {default_scoring.get((base1, base2), 0):2d} |"
        print(row)
    print()
    
    use_default = input("Use default scoring? (Y/n): ").strip().lower()
    if use_default == 'n':
        scoring = get_scoring_function()
    else:
        scoring = default_scoring.copy()
        # Ensure symmetry
        for (a, b), v in list(scoring.items()):
            scoring[(b, a)] = v

    gap_open = int(input("Enter gap opening penalty (e.g., -2): "))
    gap_extend = int(input("Enter gap extension penalty (e.g., -1): "))

    # Option for debugging mode
    debug_mode = input("Enable debug mode to show matrices? (y/N): ").strip().lower() == 'y'
    
    # Perform alignment
    alignment1, alignment2, score = affine_gap_alignment(seq1, seq2, scoring, gap_open, gap_extend, debug=debug_mode)

    print("\nOptimal Alignment:")
    visualize_alignment(alignment1, alignment2)
    print(f"Alignment Score: {score}")
    
    # Verify the score by recalculating
    verified_score = calculate_alignment_score(alignment1, alignment2, scoring, gap_open, gap_extend)
    print(f"Verified Score: {verified_score}")
    
    # Ask if user wants to run test cases
    run_tests = input("\nRun predefined test cases? (y/N): ").strip().lower() == 'y'
    if run_tests:
        run_test_cases(scoring, gap_open, gap_extend)