```markdown
# SNPsnap

**SNPsnap** is a Windows desktop application for performing multiple sequence alignment (MSA) using the MAFFT FFT-NS-1 algorithm, with a focus on discovering single-nucleotide polymorphisms (SNPs). It provides a user-friendly graphical interface (built with PyQt6) and generates a color-coded Excel report highlighting SNPs, insertions, deletions, and conserved regions.

> **Note:** This application is designed for **Windows only**. It uses bundled MAFFT Windows binaries and a batch file wrapper that relies on Windows paths and commands.

## Features

- **Multiple input methods**: Load reference and sample sequences from FASTA files or paste directly into the application.
- **Alignment engine**: Uses the powerful MAFFT aligner (FFT-NS-1) via bundled Windows binaries.
- **Interactive results**:
  - Alignment table with color-coded modifications (SNPs, insertions, deletions).
  - Clustal format output displayed in a dedicated tab.
  - Real-time logging of application events.
- **Export options**:
  - **Excel report**: Formatted worksheet with compressed conserved regions, modification types, and SNP counts per sample.
  - **Clustal file**: Save the raw alignment in Clustal format.
- **Summary statistics**: Shows mismatches, gaps, and alignment lengths for each sample.

## System Requirements

- **Operating System**: Windows 7, 8, 10, or 11 (64-bit recommended)
- **Python** 3.8 or newer (if running from source)
- **RAM**: 4 GB minimum (8+ GB recommended for large alignments)
- **Disk space**: 500 MB for installation

## Installation Options

### Option 1: Run from Source (for developers)

#### 1. Install Python dependencies
```bash
pip install PyQt6 pandas numpy xlsxwriter biopython pyinstaller
```

#### 2. Clone or download the repository
Ensure the following structure is preserved:
```
snpsnap/
├── mafft.bat
├── mafft.py
├── document.py
├── pyqt.py
├── usr/
│   ├── bin/
│   │   └── bash.exe
│   └── lib/
│       └── mafft/
│           └── ... (MAFFT scripts)
```

#### 3. Launch the application
```bash
python pyqt.py
```

### Option 2: Build a Standalone Executable (Recommended for end users)

Package SNPsnap into a single `.exe` file using PyInstaller. This creates a portable application that doesn't require Python installation.

#### Build steps:

1. **Install PyInstaller** (if not already installed):
   ```bash
   pip install pyinstaller
   ```

2. **Run the build command** from the project root:
   ```bash
   pyinstaller --onefile --windowed --name SNPsnap ^
     --add-data "mafft.bat;." ^
     --add-data "usr;usr" ^
     pyqt.py
   ```

3. **Locate the executable**:
   - The built `SNPsnap.exe` will be in the `dist/` folder
   - You can move this `.exe` anywhere on your system - it's completely portable
   - Double-click `SNPsnap.exe` to run the application

#### Build options explained:
- `--onefile`: Creates a single executable file
- `--windowed`: Prevents a console window from appearing (clean GUI experience)
- `--add-data`: Includes the MAFFT binaries and batch file in the executable
- `--name SNPsnap`: Names the output file `SNPsnap.exe`

### Option 3: Use the Pre-built Executable

If a pre-built `SNPsnap.exe` is provided in the releases section:
1. Download the `.exe` file
2. Place it in any folder (e.g., Desktop or Documents)
3. Double-click to run - no installation needed

## Usage

1. **Launch SNPsnap** by double-clicking `SNPsnap.exe` (or running `python pyqt.py` from source)

2. **Choose input mode** using the tabs at the top:
   - **File Input**: Select FASTA files from your computer
   - **Text Input**: Paste FASTA content directly

3. **Provide reference and sample sequences**:
   - **File mode**: 
     - Click "Browse..." next to Reference FASTA to select your reference file **only accepts files in .fasta .fa .fna formats!**
     - Click "Browse..." next to Sample FASTA(s) to select one or more sample files (hold Ctrl to select multiple)
   - **Text mode**:
     - Paste reference FASTA content in the top text area
     - Paste sample FASTA content in the bottom text area (can include multiple sequences)

4. Click **Run Alignment**. The application will:
   - Invoke MAFFT to perform alignment
   - Process the results (may take a few seconds depending on sequence length)
   - Display results automatically when complete

5. **Explore results** in three tabs:
   - **Alignment Table**: 
     - Conserved regions are compressed into ranges (e.g., "1-50")
     - SNPs appear in red text
     - Deletions have light red background
     - Insertions have light yellow background
     - A summary row at the bottom shows total SNPs per sample
   - **Clustal Output**: Raw alignment in standard Clustal format
   - **Logs**: Technical messages and any errors

6. **Save your work**:
   - Click **Save Excel** to generate a formatted Excel file (`.xlsx`) with all formatting preserved
   - Click **Save Clustal** to export the alignment as a text file (`.aln`)

## Input File Format

SNPsnap accepts standard FASTA format:

```
>sequence_name
ATCGATCGATCG
>another_sequence
GCTAGCTAGCTA
```

**Valid characters**: A, T, C, G, U, N, R, Y, K, M, S, W, B, D, H, V, and - (gap)

## Troubleshooting

| Problem | Solution |
|---------|----------|
| **"MAFFT not found" error** | Ensure `usr/` folder is in the same directory as the executable. When building with PyInstaller, verify `--add-data` includes all files. |
| **Excel export fails** | Install xlsxwriter: `pip install xlsxwriter`. If using standalone executable, rebuild with proper dependencies. |
| **Invalid characters in sequences** | Check your FASTA files for non-standard characters. Only IUPAC nucleotide codes and '-' are allowed. |
| **Application crashes on startup** | Ensure all bundled files are present. Try running from source to see detailed error messages. |
| **Alignment takes too long** | Reduce sequence length or number of samples. The default timeout is 5 minutes. |
| **Anti-virus blocking execution** | Add SNPsnap to your antivirus exceptions. The app uses MAFFT which may trigger false positives. |

## File Structure

```
snpsnap/
├── pyqt.py              # Main GUI application
├── mafft.py             # MAFFT wrapper and SequenceComparator class
├── document.py          # Excel report generator
├── mafft.bat            # Windows batch file to launch MAFFT
├── usr/                 # MAFFT Windows binaries
│   ├── bin/
│   │   └── bash.exe     # MSYS bash (required for MAFFT)
│   └── lib/
│       └── mafft/       # MAFFT Perl scripts and modules
└── README.md            # This file
```

## License

This project is licensed with the **MIT license**

SNPsnap incorporates **MAFFT** (version 7.505), which is also GPL-licensed. The full MAFFT license can be found at: https://mafft.cbrc.jp/alignment/software/license.html

## Acknowledgments

- **MAFFT** by Kazutaka Katoh: [https://mafft.cbrc.jp/alignment/software/](https://mafft.cbrc.jp/alignment/software/)
- **PyQt6** for the GUI framework
- **Biopython** for sequence handling
- **pandas** and **xlsxwriter** for Excel export
- **MSYS** for providing the Windows environment for MAFFT

## Support

For issues, questions, or contributions:
- Open an issue on the GitHub repository
- Include the error message from the Logs tab when reporting bugs
- Provide example input files if possible (or a minimal example that reproduces the issue)

---

*A fun and useful university project of mine built with ❤️ for the people who are also tired of relying on conventional bioinformatics tools for SNP discovery like me*
```
