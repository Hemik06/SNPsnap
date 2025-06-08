# SNPSnap/mafft.py
import glob
import subprocess
import copy
from io import StringIO
from Bio import AlignIO, SeqIO
import pandas as pd
from Bio.SeqRecord import SeqRecord
from typing import List, Optional
import os
import sys

def get_mafft_path_and_cwd():
    """
    Determines the path to the MAFFT executable and the required CWD.
    If running in a PyInstaller bundle, points to the bundled mafft.bat
    and sets CWD to the bundle's root for mafft.bat to find its resources.
    Otherwise, assumes 'mafft.bat' is in PATH
    """
    if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
        mafft_exe = os.path.join(sys._MEIPASS, 'mafft.bat')
        mafft_cwd = sys._MEIPASS
        return mafft_exe, mafft_cwd
    return 'mafft.bat', None


class MafftError(Exception):
    """Base class for MAFFT-related exceptions"""
    pass


class InputFileError(MafftError):
    """Error related to input files"""
    pass


class MafftExecutionError(MafftError):
    """Error during MAFFT execution"""
    pass


class SequenceComparator:
    THREADS = 4
    VALID_NUCLEOTIDES = set('ATCGUNRYKMSWBDHVNatcgunrykmbdhvn-')

    def __init__(self, reference: SeqRecord, samples: List[SeqRecord]):
        try:
            self.reference = copy.deepcopy(reference)
            self.reference.id = "REFERENCE"
            self.samples = [copy.deepcopy(s) for s in samples]
            self._validate_sequences()
            self.alignment, self.raw_clustal = self._perform_alignment()
        except MafftError as e:
            raise e
        except Exception as e:
            raise MafftError(f"Unexpected error during initialization: {str(e)}") from e

    def _validate_sequences(self):
        all_records = [self.reference] + self.samples
        for record in all_records:
            if not record.seq:
                raise InputFileError(f"Sequence for '{record.id}' is empty.")
            sequence_str = str(record.seq)
            invalid_chars = set(char for char in sequence_str if char not in self.VALID_NUCLEOTIDES)
            if invalid_chars:
                raise InputFileError(
                    f"Invalid characters {sorted(list(invalid_chars))} found in sequence '{record.id}'. "
                    f"Allowed characters: {''.join(sorted(list(self.VALID_NUCLEOTIDES)))}."
                )

    @classmethod
    def from_text(cls, reference_text: str, samples_text: str):
        try:
            ref_records = list(SeqIO.parse(StringIO(reference_text), "fasta"))
            if not ref_records:
                raise InputFileError("Reference FASTA text is empty or invalid.")
            if len(ref_records) != 1:
                raise InputFileError("Reference must contain exactly one sequence.")
            ref = ref_records[0]
            samples = list(SeqIO.parse(StringIO(samples_text), "fasta"))
            if not samples:
                raise InputFileError("No valid sample sequences found in FASTA text.")
            return cls(ref, samples)
        except InputFileError:
            raise
        except Exception as e:
            raise InputFileError(f"Invalid FASTA text format: {str(e)}") from e

    @classmethod
    def from_files(cls, ref_path: str, sample_paths: List[str]):
        try:
            if not os.path.exists(ref_path):
                raise InputFileError(f"Reference file not found: {ref_path}")
            if not os.path.isfile(ref_path):
                raise InputFileError(f"Reference path is not a file: {ref_path}")
            try:
                ref = SeqIO.read(ref_path, "fasta")
            except ValueError as e:
                 if "More than one record found" in str(e):
                     raise InputFileError(f"Reference file '{ref_path}' must contain exactly one sequence.") from e
                 raise InputFileError(f"Invalid reference FASTA file ('{ref_path}'): {str(e)}") from e
            except Exception as e:
                raise InputFileError(f"Error reading reference FASTA file ('{ref_path}'): {str(e)}") from e

            samples = []
            if not sample_paths:
                raise InputFileError("No sample file paths provided.")
            for path in sample_paths:
                if not os.path.exists(path):
                    raise InputFileError(f"Sample file not found: {path}")
                if not os.path.isfile(path):
                    raise InputFileError(f"Sample path is not a file: {path}")
                try:
                    current_sample_records = list(SeqIO.parse(path, "fasta"))
                    if not current_sample_records:
                        raise InputFileError(f"Sample FASTA file '{path}' is empty or invalid.")
                    if len(current_sample_records) > 1:
                         raise InputFileError(f"Sample file '{path}' contains multiple sequences. Each sample file should contain one sequence.")
                    samples.append(current_sample_records[0])
                except InputFileError:
                    raise
                except Exception as e:
                    raise InputFileError(f"Invalid sample FASTA file ('{path}'): {str(e)}") from e
            if not samples:
                raise InputFileError("No valid sample sequences could be read from the provided files.")
            return cls(ref, samples)
        except MafftError:
            raise
        except Exception as e:
            raise MafftError(f"File loading failed: {str(e)}") from e

    def _perform_alignment(self) -> tuple[AlignIO.MultipleSeqAlignment, str]:
        if not self.samples:
            raise InputFileError("Cannot perform alignment: No sample sequences loaded.")

        mafft_executable, mafft_cwd_for_subprocess = get_mafft_path_and_cwd()

        try:
            all_sequences = [self.reference] + self.samples
            with StringIO() as buffer:
                SeqIO.write(all_sequences, buffer, "fasta")
                input_data = buffer.getvalue()

            if not input_data.strip():
                 raise MafftExecutionError("Internal error: Generated FASTA input for MAFFT is empty.")

            cmd = [mafft_executable, "--auto", "--thread", str(self.THREADS), "--clustalout", "-"]

            # START MODIFICATION: Hide console window on Windows
            creation_flags = 0
            if sys.platform == "win32":
                # 0x08000000 means CREATE_NO_WINDOW
                # This flag is only used when the process is created directly.
                # If mafft.bat itself calls other console programs, those might still appear briefly
                # unless mafft.bat itself is modified or those tools are also launched with no window.
                # However, for mafft.bat itself, this should hide its main window.
                creation_flags = subprocess.CREATE_NO_WINDOW
            # END MODIFICATION

            result = subprocess.run(
                cmd,
                input=input_data,
                text=True,
                capture_output=True,
                check=False,
                timeout=300,
                cwd=mafft_cwd_for_subprocess,
                creationflags=creation_flags # MODIFIED: Added creationflags
            )

            if result.returncode != 0:
                error_msg = f"MAFFT failed (exit code {result.returncode}):\nSTDERR:\n{result.stderr}\nSTDOUT:\n{result.stdout}"
                if "Empty input" in result.stderr or "empty" in result.stderr.lower():
                    error_msg = "MAFFT reported empty input. Check your FASTA files for valid sequences."
                elif "contains duplicated identical sequences" in result.stderr:
                     error_msg = f"MAFFT error: Input contains duplicated identical sequences. {result.stderr}"
                raise MafftExecutionError(error_msg)

            if "WARNING" in result.stderr:
                if "contains invalid characters" in result.stderr:
                    raise MafftExecutionError("MAFFT warning: Input contains invalid nucleotide characters. Please check sequence data.")
                if "too short" in result.stderr:
                    raise MafftExecutionError("MAFFT warning: Input sequences might be too short for reliable alignment.")

            clustal_output = result.stdout
            if not clustal_output.strip():
                raise MafftExecutionError("MAFFT produced empty CLUSTAL output despite success code. stderr: " + result.stderr)

            try:
                alignment = AlignIO.read(StringIO(clustal_output), "clustal")
            except ValueError as e:
                raise MafftExecutionError(f"Failed to parse MAFFT's CLUSTAL output: {str(e)}\nOutput:\n{clustal_output[:1000]}") from e

            if len(alignment) != len(all_sequences):
                raise MafftExecutionError(
                    f"Alignment record count ({len(alignment)}) does not match input sequence count ({len(all_sequences)}). "
                    "This may indicate an issue with MAFFT processing or sequence IDs."
                )

            return alignment, clustal_output
        except subprocess.TimeoutExpired:
            raise MafftExecutionError("MAFFT alignment timed out after 5 minutes. Sequences may be too long or too numerous.")
        except FileNotFoundError:
            if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
                raise MafftExecutionError(f"MAFFT executable ({mafft_executable}) not found or not executable within the application bundle. Ensure it was correctly bundled and permissions are set.")
            else:
                raise MafftExecutionError("MAFFT executable ('mafft.bat') not found. Ensure MAFFT is installed and 'mafft.bat' is in your system's PATH.")
        except MafftError:
            raise
        except Exception as e:
            raise MafftExecutionError(f"An unexpected error occurred during MAFFT alignment: {str(e)}") from e

    def write_clustal(self, filename: str) -> None:
        try:
            with open(filename, 'w') as f:
                f.write(self.raw_clustal)
        except IOError as e:
            raise MafftError(f"Error writing Clustal file '{filename}': {e}") from e

    @property
    def nucleotide_matrix(self) -> pd.DataFrame:
        if not self.alignment:
             raise MafftError("Alignment has not been performed or failed.")
        aligned_records = {rec.id: str(rec.seq) for rec in self.alignment}
        if "REFERENCE" not in aligned_records:
            raise MafftError("Reference sequence not found in alignment output. Expected ID 'REFERENCE'.")
        ref_seq_str = aligned_records['REFERENCE']
        data_for_df = {}
        data_for_df['REFERENCE'] = list(ref_seq_str)
        for rec_id, seq_str in aligned_records.items():
            if rec_id != 'REFERENCE':
                data_for_df[rec_id] = list(seq_str)
        positions: List[Optional[int]] = []
        counter = 0
        for base in ref_seq_str:
            if base == '-':
                positions.append(None)
            else:
                counter += 1
                positions.append(counter)
        df = pd.DataFrame(
            data_for_df,
            index=pd.RangeIndex(1, len(ref_seq_str) + 1, name='AlignmentIndex')
        )
        df.insert(0, 'Position', positions)
        df['Position'] = pd.to_numeric(df['Position'], errors='coerce')
        return df

    @property
    def summary(self) -> pd.DataFrame:
        if not self.alignment:
             raise MafftError("Alignment has not been performed or failed.")
        ref_record_aligned = None
        for rec in self.alignment:
            if rec.id == "REFERENCE":
                ref_record_aligned = rec
                break
        if ref_record_aligned is None:
            raise MafftError("Reference sequence (ID 'REFERENCE') not found in the alignment for summary generation.")
        ref_str = str(ref_record_aligned.seq)
        stats = []
        for rec in self.alignment:
            mismatches = 0
            gaps_in_comparison = 0 # Initialize for reference case
            if rec.id == "REFERENCE":
                pass # mismatches and gaps_in_comparison already 0 or N/A
            else:
                seq_str = str(rec.seq)
                mismatches = sum(1 for r_base, s_base in zip(ref_str, seq_str) if r_base != s_base and r_base != '-' and s_base != '-')
                gaps_in_comparison = sum(1 for r_base, s_base in zip(ref_str, seq_str) if r_base == '-' or s_base == '-')
            stats.append({
                'ID': rec.id,
                'Mismatches_vs_Ref': mismatches if rec.id != "REFERENCE" else 'N/A',
                'Total_Gaps_In_Seq': str(rec.seq).count('-'),
                'Gaps_vs_Ref': gaps_in_comparison if rec.id != "REFERENCE" else 'N/A',
                'Aligned_Sequence_Length': len(rec.seq),
            })
        return pd.DataFrame(stats)


if __name__ == '__main__':
    
    try:
        print("Attempting alignment with valid files...")
        
        reference_file = "agr14.fasta" # Or your reference fasta file
        sample_files_pattern = "samples/*.fasta" # Or your pattern

        if not os.path.exists(reference_file):
            print(f"Warning: Reference file '{reference_file}' not found for example.")
        else:
            sample_files = glob.glob(sample_files_pattern)
            if not sample_files:
                 print(f"Warning: No sample files found matching '{sample_files_pattern}' for example.")
            else:
                comparator = SequenceComparator.from_files(
                    reference_file,
                    sample_files
                )
                matrix = comparator.nucleotide_matrix
                print("\nNucleotide Matrix:")
                print(matrix.head()) # Print head to keep it concise
                print("\nSummary Statistics:")
                print(comparator.summary)
                comparator.write_clustal("example_alignment.aln")
                print("\nWrote example_alignment.aln")

        # Example of trying with an invalid file (if _validate_sequences is working)
        # Create a dummy 'reference.fasta' and 'sample_invalid_chars.fasta' for this test
        # e.g., sample_invalid_chars.fasta: >invalid\nACGTX
        if os.path.exists("reference.fasta") and os.path.exists("sample_invalid_chars.fasta"):
            print("\nAttempting alignment with an invalid character file (should fail)...")
            try:
                SequenceComparator.from_files("reference.fasta", ["sample_invalid_chars.fasta"])
            except InputFileError as e:
                print(f"Caught expected InputFileError: {e}")
        else:
            print("\nSkipping invalid char test: 'reference.fasta' or 'sample_invalid_chars.fasta' not found.")

    except MafftError as e:
        print(f"MAFFT Error in example: {e}")
    except Exception as e:
        print(f"Unexpected error in example: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()