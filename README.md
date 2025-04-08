# Pairwise Sequence Alignment using Affine Gap Penalty (Python CLI Tool)

## 📌 Project Overview
This project implements a **command-line tool** in Python to perform **pairwise alignment** of two DNA sequences using an **affine gap penalty model**. The tool supports custom scoring for character matches/mismatches and handles penalties for both **gap opening** and **gap extension**, providing more biologically realistic alignments.

This was created as a project solution for a bioinformatics or algorithms course, and it supports both basic and bonus requirements described in the assignment PDF.


## 🚀 Features
- Accepts two DNA sequences as user input
- Accepts a **custom scoring matrix** (e.g., A-C = 2)
- Accepts **gap opening** and **gap extension penalties**
- Implements the **affine gap alignment algorithm** using dynamic programming
- Displays the optimal alignment and final score


## 🧠 How It Works

### 1. **User Input**
The CLI tool prompts the user to:
- Enter two sequences (sequence1 and sequence2)
- Input scoring values (e.g., A A 2 for a match, A C -1 for a mismatch)
- Enter gap opening and gap extension penalties

### 2. **Scoring Matrix**
A dictionary is created to store symmetric match/mismatch scores (e.g., `scoring['AC'] = -1`, also assigns `scoring['CA'] = -1`).

### 3. **Affine Gap Algorithm**
We use **three matrices**:
- `M[i][j]`: optimal alignment score ending with match/mismatch
- `I[i][j]`: alignment ending in a gap in sequence 1 (insertion)
- `D[i][j]`: alignment ending in a gap in sequence 2 (deletion)

Each cell is filled based on recursive relationships:
```python
M[i][j] = max(M[i-1][j-1], I[i-1][j-1], D[i-1][j-1]) + match_score
I[i][j] = max(M[i][j-1] + gap_open, I[i][j-1] + gap_extend)
D[i][j] = max(M[i-1][j] + gap_open, D[i-1][j] + gap_extend)
```

### 4. **Traceback**
After filling the matrices, we backtrack from the bottom-right cell to reconstruct the aligned sequences by tracing which move led to each score (match/mismatch, gap in seq1, gap in seq2).

### 5. **Output**
The tool prints:
- Aligned sequence 1
- Aligned sequence 2
- Final alignment score



## 🧪 Sample Run
```
Enter first DNA sequence: ACGT
Enter second DNA sequence: AGT
Enter scoring function (format: A C 2). Type 'done' to finish:
A A 2
C C 2
G G 2
T T 2
A G -1
A C -1
C G -1
G T -1
done
Enter gap opening penalty: -2
Enter gap extension penalty: -1

Optimal Alignment:
ACGT
A-GT
Alignment Score: 3
```



## 🛠️ Requirements
- Python 3.x
- Numpy



## 🏗️ How to Run
```bash
git clone https://github.com/GRaman8/Pairwise-Alignment-Term-Project.git
cd Pairwise-Alignment-Term-Project
python3 aligner.py
```



## 📁 File Structure
```
sequence-aligner/
├── aligner.py           # Main Python script
├── Project.pdf          # Project description           
├── README.md            # This file
├── .gitignore           # Ignore virtual env, __pycache__, etc.
├──  requirements.txt    # Contains python packages
```




## 📜 License
This project is licensed under the MIT License.

