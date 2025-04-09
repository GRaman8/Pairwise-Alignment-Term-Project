# Pairwise Sequence Alignment using Affine Gap Penalty (Python CLI Tool)

## Overview
This project implements a **command-line tool** in Python to perform **pairwise alignment** of two DNA sequences using an **affine gap penalty model**. The tool supports custom scoring for character matches/mismatches and handles penalties for both **gap opening** and **gap extension**, providing more biologically realistic alignments.

This was created as a project solution for a bioinformatics or algorithms course, and it supports both basic and bonus requirements described in the assignment PDF.

## Features

- Align two DNA sequences with customizable scoring parameters
- Support for affine gap penalties (different costs for opening vs extending gaps)
- Visualize the alignment with match/mismatch indicators
- Debug mode to view dynamic programming matrices
- Built-in score verification
- Test cases to verify implementation correctness

## Algorithm Details

The implementation uses a three-matrix approach to handle affine gap penalties:

1. **Match/Mismatch Matrix (M)**: Tracks alignments ending with matched/mismatched characters
2. **Insertion Matrix (I)**: Tracks alignments ending with a gap in sequence 1
3. **Deletion Matrix (D)**: Tracks alignments ending with a gap in sequence 2

### Dynamic Programming Approach

For each position (i,j) in the alignment:

1. **M[i][j]**: Max score of alignments ending with x[i] aligned to y[j]
   - From M: M[i-1][j-1] + match_score(x[i], y[j])
   - From I: I[i-1][j-1] + match_score(x[i], y[j])
   - From D: D[i-1][j-1] + match_score(x[i], y[j])

2. **I[i][j]**: Max score of alignments ending with a gap in x aligned to y[j]
   - From M: M[i][j-1] + gap_open
   - From I: I[i][j-1] + gap_extend

3. **D[i][j]**: Max score of alignments ending with x[i] aligned to a gap in y
   - From M: M[i-1][j] + gap_open
   - From D: D[i-1][j] + gap_extend

### Traceback

After filling the matrices, the algorithm traces back from the maximum score to reconstruct the optimal alignment, transitioning between matrices based on the stored traceback information.

## Implementation

The code is implemented in Python using NumPy for efficient matrix operations. Key components include:

- **Matrix Initialization**: Special handling for first row/column to account for gap penalties
- **Traceback Matrices**: Tracking the source of each optimal score to reconstruct the alignment
- **Affine Gap Model**: Separate gap opening (-2 by default) and gap extension (-1 by default) penalties
- **Scoring Function**: Customizable scoring for matches, mismatches, and gaps

## How to Run
```bash
git clone https://github.com/GRaman8/Pairwise-Alignment-Term-Project.git
cd Pairwise-Alignment-Term-Project
python3 aligner.py
```

### Input

1. Enter two DNA sequences (only A, C, G, T characters allowed)
2. Choose to use the default scoring matrix or provide custom values
3. Specify gap opening and extension penalties
4. Choose whether to enable debug mode

### Example

```
Enter first DNA sequence: ACGTTGA
Enter second DNA sequence: GCTAG

Default scoring matrix:
  | A  | C  | G  | T  |
--------------------- 
A |  2 | -1 | -1 | -1 |
C |  0 |  2 | -1 | -1 |
G |  0 |  0 |  2 | -1 |
T |  0 |  0 |  0 |  2 |

Use default scoring? (Y/n): Y
Enter gap opening penalty (e.g., -2): -2
Enter gap extension penalty (e.g., -1): -1
Enable debug mode to show matrices? (y/N): y

Optimal Alignment:
ACGTTGA
x| |x| 
GC-TAG-
Alignment Score: 0
Verified Score: 0
```

## Working Process

1. **Initialization**: 
   - Set up three matrices (M, I, D) for tracking different alignment states
   - Initialize boundaries with appropriate gap penalties

2. **Fill Matrices**:
   - For each position, calculate optimal scores for:
     - Aligning two characters (match/mismatch)
     - Inserting a gap in sequence 1
     - Inserting a gap in sequence 2
   - Store the path taken in traceback matrices

3. **Traceback**:
   - Begin at the cell with maximum score
   - Follow the traceback matrices to reconstruct the alignment
   - Transition between matrices (M, I, D) as needed

4. **Verification**:
   - Recalculate the score from the alignment to verify correctness
   - Compare with the score from the dynamic programming process

## Testing

The tool includes built-in test cases to verify correctness:
- Simple cases with one gap
- Cases with multiple gaps
- More complex alignments

Each test case shows the alignment, score, and verification status.

## Understanding the Output

- The alignment shows both sequences with gaps inserted as needed
- The middle line indicates the relationship between aligned positions:
  - `|` for matching characters
  - `x` for mismatches
  - ` ` (space) for positions with gaps

## Time and Space Complexity

- **Time Complexity**: O(mn) where m and n are the lengths of the sequences
- **Space Complexity**: O(mn) for storing the dynamic programming matrices

## Affine Gap Model

The affine gap model uses a different penalty for opening a gap versus extending an existing gap. This is biologically more accurate as it's often more likely for a gap to be extended than for a new gap to be created.

Formula: gap_penalty = gap_open + (gap_length - 1) * gap_extend

## Requirements
- Python 3.x
- Numpy
 
## File Structure
```
sequence-aligner/
├── aligner.py           # Main Python script
├── Project.pdf          # Project description           
├── README.md            # This file
├── .gitignore           # Ignore virtual env, __pycache__, etc.
├──  requirements.txt    # Contains python packages
├──  LICENSE             # Contains Licensing information
```

## License
This project is licensed under the MIT License.

## References

- Needleman, S.B. and Wunsch, C.D. (1970) A general method applicable to the search for similarities in the amino acid sequence of two proteins. Journal of Molecular Biology, 48, 443-453.
- Gotoh, O. (1982) An improved algorithm for matching biological sequences. Journal of Molecular Biology, 162, 705-708.