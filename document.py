import logging
import sys
from typing import List, Optional

import glob
import numpy as np
import pandas as pd
try:
    import xlsxwriter
except ImportError:
    print("Error: 'xlsxwriter' library is required for Excel export. "
          "Please install it using 'pip install xlsxwriter'")
    sys.exit(1)
# Assuming mafft SequenceComparator exists and works as intended
from mafft import SequenceComparator


# Setup basic logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')


class AlignmentExporter:
    """
    Exports a nucleotide alignment DataFrame to an Excel file.

    Highlights differences (SNPs, Insertions, Deletions) compared to a
    reference sequence and compresses conserved regions for readability.
    """

    DEFAULT_REF_COL = 'REFERENCE'
    DEFAULT_POS_COL = 'Position'
    MOD_TYPE_COL = 'Modification Type'

    def __init__(
        self,
        matrix: pd.DataFrame,
        reference_column: str = DEFAULT_REF_COL,
        position_column: str = DEFAULT_POS_COL
    ):
        """
        Initializes the exporter with the alignment data.

        Args:
            matrix: DataFrame containing alignment data. Must include position,
                    reference, and sample columns.
            reference_column: Name of the column containing reference bases.
            position_column: Name of the column containing sequence positions.

        Raises:
            TypeError: If matrix is not a pandas DataFrame.
            ValueError: If matrix is empty, or required columns are missing.
        """
        if not isinstance(matrix, pd.DataFrame):
            raise TypeError("Input 'matrix' must be a pandas DataFrame.")
        if matrix.empty:
            raise ValueError("Input 'matrix' cannot be empty.")

        self.df = matrix.copy()
        self.ref_col = reference_column

        # Standardize position column name, allowing for legacy name
        if 'BiologicalPosition' in self.df.columns and position_column == self.DEFAULT_POS_COL:
            logging.info("Renaming 'BiologicalPosition' column to '%s'", self.DEFAULT_POS_COL)
            self.df.rename(columns={'BiologicalPosition': self.DEFAULT_POS_COL}, inplace=True)
        self.pos_col = self.DEFAULT_POS_COL # Use standardized name internally

        # Validate required columns exist
        if self.ref_col not in self.df.columns:
            raise ValueError(f"Reference column '{self.ref_col}' not found in DataFrame.")
        if self.pos_col not in self.df.columns:
            raise ValueError(f"Position column '{self.pos_col}' not found in DataFrame.")

        self.sample_cols = [
            c for c in self.df.columns if c not in (self.pos_col, self.ref_col, self.MOD_TYPE_COL)
        ]
        if not self.sample_cols:
            logging.warning("No sample columns found in the input DataFrame.")

        # Ensure base columns are string type and handle potential NaNs
        self.df[self.ref_col] = self.df[self.ref_col].fillna('-').astype(str).str.upper()
        for col in self.sample_cols:
            self.df[col] = self.df[col].fillna('-').astype(str).str.upper()
        # Ensure position column is suitable for numeric operations later
        self.df[self.pos_col] = pd.to_numeric(self.df[self.pos_col], errors='coerce')


        # Determine modification types using a vectorized approach
        self.df[self.MOD_TYPE_COL] = self._determine_modification_types()

        # Blank out sample bases that exactly match the reference base,
        # unless the reference base itself is an alignment gap ('-').
        for col in self.sample_cols:
            mask = (self.df[col] == self.df[self.ref_col]) & (self.df[self.ref_col] != '-')
            self.df.loc[mask, col] = ''

    def _determine_modification_types(self) -> pd.Series:
        """
        Determines the modification type for each alignment position
        using vectorized operations. Priority: Insertion > Deletion > SNP > Conserved.
        """
        ref = self.df[self.ref_col]
        # Handle case with no sample columns gracefully
        if not self.sample_cols:
             is_insertion = (ref == '-')
             mod_types = pd.Series('Conserved', index=self.df.index)
             mod_types.loc[is_insertion] = 'Insertion'
             return mod_types
             
        samples_df = self.df[self.sample_cols]

        is_insertion = (ref == '-')

        # Deletion: any sample has '-' where ref does not
        has_deletion = ((samples_df == '-').any(axis=1)) & (~is_insertion)

        # SNP: any sample differs from ref (''), excluding '-' for both ref and sample
        differs_mask = samples_df.ne(ref, axis=0)
        is_valid_sample_base = (samples_df != '') & (samples_df != '-')
        has_snp = (differs_mask & is_valid_sample_base).any(axis=1) & (~is_insertion) & (~has_deletion)

        # Assign types based on priority
        mod_types = pd.Series('Conserved', index=self.df.index)
        mod_types.loc[has_snp] = 'SNP'
        mod_types.loc[has_deletion] = 'Deletion'
        mod_types.loc[is_insertion] = 'Insertion' # Highest priority

        return mod_types

    @staticmethod
    def _compress_regions(indices: List[int]) -> List[str]:
        """Compresses consecutive integer indices into ranges."""
        if not indices:
            return []
        indices = sorted(list(set(indices))) # Ensure unique sorted list
        if not indices:
            return []
            
        cuts = np.where(np.diff(indices) != 1)[0] + 1
        groups = np.split(np.array(indices), cuts)
        return [
            f"{g[0]}-{g[-1]}" if len(g) > 1 else str(g[0])
            for g in groups if len(g) > 0 # Ensure group is not empty
        ]

    def _write_conserved_buffer(
        self,
        ws: 'xlsxwriter.worksheet.Worksheet',
        buf_indices: List[int],
        current_row_idx: int,
        fmt_cons: 'xlsxwriter.format.Format'
    ) -> tuple[int, List[list]]:
        """Writes compressed conserved regions from the buffer to the worksheet."""
        if not buf_indices:
            return current_row_idx, []

        # Get valid numeric positions from the buffer indices
        pos_vals = [
            int(x) for x in self.df.loc[buf_indices, self.pos_col]
            if pd.notna(x)
        ]
        
        output_rows_segment = []
        compressed_ranges = self._compress_regions(pos_vals)

        for rng in compressed_ranges:
            new_row_data = [rng, '', 'Conserved'] + [''] * len(self.sample_cols)
            output_rows_segment.append(new_row_data)
            ws.write_row(current_row_idx, 0, new_row_data, fmt_cons)
            current_row_idx += 1

        return current_row_idx, output_rows_segment

    def export_excel(self, filename: str) -> Optional[pd.DataFrame]:
        """
        Exports the processed alignment data to an Excel file.

        Args:
            filename: Path to the output Excel file (.xlsx).

        Returns:
            A DataFrame mirroring the structure written to Excel, or None if
            an error occurs during Excel writing.
        """
        try:
            with pd.ExcelWriter(filename, engine='xlsxwriter') as writer:
                wb = writer.book
                ws = wb.add_worksheet('Alignment')

                # Define formats
                fmt_hdr = wb.add_format({'bold': True, 'bg_color': '#DDEBF7', 'border': 1})
                fmt_snp = wb.add_format({'font_color': '#FF0000'})  # Red text for SNPs
                fmt_deletion = wb.add_format({'bg_color': '#FFC7CE'})  # Light red background for deletions
                fmt_deletion_blue = wb.add_format({'bg_color': '#FFC7CE', 'font_color': '#0000FF'})  # Blue text on red for deletions
                fmt_insertion = wb.add_format({'bg_color': '#FFFFCC'})  # Light yellow background for insertions
                fmt_cons = wb.add_format({'bg_color': '#E2EFDA'})  # Light green for conserved regions

                # Write headers
                headers = [self.pos_col, self.ref_col, self.MOD_TYPE_COL] + self.sample_cols
                ws.write_row(0, 0, headers, fmt_hdr)

                row_idx = 1  # Excel row index starts from 1 (0 is header)
                conserved_buffer = []
                all_output_rows = []  # Collect all rows for the final DataFrame

                # Iterate through DataFrame rows by index
                for i in self.df.index:
                    mtype = self.df.at[i, self.MOD_TYPE_COL]

                    if mtype == 'Conserved':
                        conserved_buffer.append(i)
                        continue

                    # Process conserved buffer if we hit a non-conserved row
                    if conserved_buffer:
                        row_idx, added_rows = self._write_conserved_buffer(
                            ws, conserved_buffer, row_idx, fmt_cons
                        )
                        all_output_rows.extend(added_rows)
                        conserved_buffer = []  # Clear buffer

                    # Process the current non-conserved row
                    pos_val = self.df.at[i, self.pos_col]
                    pos_display = str(int(pos_val)) if pd.notna(pos_val) else ''
                    ref_display = self.df.at[i, self.ref_col]
                    sample_values = [self.df.at[i, s] for s in self.sample_cols]
                    current_row_data = [pos_display, ref_display, mtype] + sample_values
                    all_output_rows.append(current_row_data)

                    # Determine row format based on modification type
                    if mtype == 'Deletion':
                        row_fmt = fmt_deletion
                    elif mtype == 'Insertion':
                        row_fmt = fmt_insertion
                    else:
                        row_fmt = None  # No background for SNPs

                    # Write position, reference, and modification type with row format
                    ws.write(row_idx, 0, pos_display, row_fmt)
                    ws.write(row_idx, 1, ref_display, row_fmt)
                    ws.write(row_idx, 2, mtype, row_fmt)

                    # Write sample data with appropriate formatting
                    for col_offset, cell_value in enumerate(sample_values):
                        col_idx = 3 + col_offset
                        cell_fmt = None

                        if mtype == 'Deletion':
                            if cell_value == '-':
                                cell_fmt = fmt_deletion_blue
                            else:
                                cell_fmt = row_fmt
                        elif mtype == 'Insertion':
                            cell_fmt = row_fmt
                        elif mtype == 'SNP' and cell_value != '':
                            cell_fmt = fmt_snp

                        ws.write(row_idx, col_idx, cell_value, cell_fmt)

                    row_idx += 1

                # Process any remaining conserved buffer after the loop
                if conserved_buffer:
                    row_idx, added_rows = self._write_conserved_buffer(
                        ws, conserved_buffer, row_idx, fmt_cons
                    )
                    all_output_rows.extend(added_rows)

                # Set column widths for better readability
                ws.set_column(0, 0, 12)  # Position
                ws.set_column(1, 1, 10)  # Reference
                ws.set_column(2, 2, 18)  # Modification Type
                if self.sample_cols:
                    ws.set_column(3, 2 + len(self.sample_cols), 8)  # Sample columns

                # Auto-filter
                ws.autofilter(0, 0, row_idx - 1, len(headers) - 1)

                # Freeze top row
                ws.freeze_panes(1, 0)

                # Add sum of SNPs row
                if self.sample_cols:
                    # Calculate SNP counts for each sample
                    snp_counts = {}
                    for sample in self.sample_cols:
                        mask = (self.df[self.MOD_TYPE_COL] == 'SNP') & (self.df[sample] != '')
                        snp_counts[sample] = mask.sum()
                    # Prepare sum row data
                    sum_row_data = ["sum of SNPs", "Ref", ""] + [snp_counts.get(s, 0) for s in self.sample_cols]
                    # Create format for the summary row
                    fmt_summary = wb.add_format({'bg_color': '#B7DEE8'})  # Pastel blue background
                    # Write the row
                    ws.write_row(row_idx, 0, sum_row_data, fmt_summary)
                    all_output_rows.append(sum_row_data)
                    row_idx += 1

                logging.info("Successfully exported alignment report to %s", filename)
                return pd.DataFrame(all_output_rows, columns=headers)

        except Exception as e:
            logging.error("Failed to export alignment to Excel: %s", e, exc_info=True)
            return None

def export_alignment(matrix: pd.DataFrame, filename: str) -> Optional[pd.DataFrame]:
    """
    Convenience function to create an AlignmentExporter and export the matrix.

    Args:
        matrix: DataFrame containing alignment data.
        filename: Path to the output Excel file (.xlsx).

    Returns:
        DataFrame mirroring the Excel output, or None on failure.
    """
    try:
        exporter = AlignmentExporter(matrix)
        return exporter.export_excel(filename)
    except (ValueError, TypeError) as e:
        logging.error("Error initializing AlignmentExporter: %s", e)
        return None
    except Exception as e: # Catch potential issues during export not caught inside
        logging.error("An unexpected error occurred during alignment export: %s", e)
        return None


# Example Usage (requires placeholder mafft and fasta files):
if __name__ == "__main__":
    # --- Original example using hypothetical SequenceComparator ---
    # try:
    # Placeholder: Replace with actual file paths if running
    # Ensure 'reference.fasta' and 'sample*.fasta' exist
    sample_files = glob.glob("sample[1-4].fasta") # Use valid samples
        
    comparator = SequenceComparator.from_files("reference.fasta", sample_files)
    alignment_matrix = comparator.nucleotide_matrix

    # Using dummy data instead for demonstration
    alignment_matrix = alignment_matrix

    print("Exporting alignment report...")
    export_alignment(alignment_matrix, "alignment_report_from_comparator.xlsx")

    # Example of writing Clustal (assuming comparator has this method)
    print("Writing CLUSTAL file...")
    comparator.write_clustal("alignment.clustal")
    
    # except FileNotFoundError:
    #     print("Error: Fasta files for SequenceComparator not found.")
    # except NameError:
    #     print("Error: SequenceComparator class not available.")
    # except Exception as e:
    #     print(f"An error occurred: {e}")