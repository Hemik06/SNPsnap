# SNPSnap/pyqt.py
import sys
import logging
import os # Ensure os is imported
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
                             QPushButton, QFileDialog, QLabel, QLineEdit, QTabWidget,
                             QTableWidget, QTableWidgetItem, QTextEdit, QMessageBox,
                             QPlainTextEdit)
from PyQt6.QtCore import QThread, pyqtSignal, Qt
from PyQt6.QtGui import QFont

# --- Placeholder for external modules (for testing if not available) ---
try:
    from document import export_alignment
    from mafft import SequenceComparator, InputFileError, MafftExecutionError, MafftError
except ImportError:
    # Create stubs if the actual modules are not available
    logging.warning("Using placeholder stubs for 'document' and 'mafft' modules.")
    import pandas as pd

    def export_alignment(matrix_data, output_path):
        if isinstance(matrix_data, list) and len(matrix_data) > 1:
            df = pd.DataFrame(matrix_data[1:], columns=matrix_data[0])
        elif isinstance(matrix_data, pd.DataFrame):
            df = matrix_data # In our case, matrix_data from SequenceComparator is already a DataFrame
        else:
            df = pd.DataFrame(matrix_data or [])
        
        try:
            # The stub here is simplified. The real function does more.
            df.to_excel(output_path, index=False) 
            logging.info(f"Stub: Matrix exported to {output_path}")
        except Exception as e:
            logging.error(f"Stub: Error exporting matrix to {output_path}: {e}")
        return df # Crucially, it returns the DataFrame

    class MafftError(Exception): pass
    class InputFileError(MafftError): pass
    class MafftExecutionError(MafftError): pass

    class SequenceComparator:
        def __init__(self, ref_seq_record, sample_seq_records): # Adjusted for SeqRecord inputs
            self.ref_seq = str(ref_seq_record.seq)
            self.sample_seqs = [str(s.seq) for s in sample_seq_records]
            
            if not self.ref_seq or not self.sample_seqs:
                 raise InputFileError("Stub: Empty reference or sample sequences.")
            if "INPUT_ERROR" in self.ref_seq: raise InputFileError("Stub: Simulated input error")
            if "EXEC_ERROR" in self.ref_seq: raise MafftExecutionError("Stub: Simulated exec error")
            if "MAFFT_ERROR" in self.ref_seq: raise MafftError("Stub: Simulated mafft error")

            self._raw_clustal = f">REFERENCE\n{self.ref_seq}\n"
            for i, s_text in enumerate(self.sample_seqs):
                self._raw_clustal += f">sample{i+1}\n{s_text}\n"
            
            # Simulate nucleotide_matrix more closely to actual output
            cols = ['Position', 'REFERENCE'] + [f'sample{i+1}' for i in range(len(self.sample_seqs))]
            # This is a simplified matrix for stubbing
            data = {'Position': list(range(1, len(self.ref_seq) + 1)),
                    'REFERENCE': list(self.ref_seq)}
            for i, s_text in enumerate(self.sample_seqs):
                data[f'sample{i+1}'] = list(s_text.ljust(len(self.ref_seq),'-')[:len(self.ref_seq)])
            
            self._nucleotide_matrix_df = pd.DataFrame(data, columns=cols)


        @property
        def raw_clustal(self): return self._raw_clustal
        
        @property
        def nucleotide_matrix(self): # This is what's passed to export_alignment
            return self._nucleotide_matrix_df 

        @classmethod
        def _parse_fasta_text_to_seqrecords(cls, text_content, is_ref=False):
            from io import StringIO
            from Bio import SeqIO
            from Bio.SeqRecord import SeqRecord

            records = list(SeqIO.parse(StringIO(text_content), "fasta"))
            if is_ref:
                if not records: raise InputFileError("Stub: No sequence found in reference text.")
                if len(records) !=1: raise InputFileError("Stub: Reference must be single sequence.")
                return records[0]
            if not records: raise InputFileError("Stub: No sequences found in sample text.")
            return records


        @classmethod
        def from_files(cls, ref_file_path, sample_file_paths):
            from Bio.SeqRecord import SeqRecord # ensure SeqRecord is available for stub
            try:
                with open(ref_file_path, 'r') as f: ref_text = f.read()
                ref_rec = cls._parse_fasta_text_to_seqrecords(ref_text, is_ref=True)
            except Exception as e: raise InputFileError(f"Stub: Error reading ref file {ref_file_path}: {e}")
            
            sample_recs = []
            for p in sample_file_paths:
                try:
                    with open(p, 'r') as f: sample_text = f.read()
                    sample_recs.extend(cls._parse_fasta_text_to_seqrecords(sample_text))
                except Exception as e: raise InputFileError(f"Stub: Error reading sample file {p}: {e}")
            
            if not sample_recs: raise InputFileError("Stub: No sequences found in sample files.")
            return cls(ref_rec, sample_recs)

        @classmethod
        def from_text(cls, ref_data_text, sample_data_text):
            from Bio.SeqRecord import SeqRecord # ensure SeqRecord is available for stub
            ref_rec = cls._parse_fasta_text_to_seqrecords(ref_data_text, is_ref=True)
            sample_recs = cls._parse_fasta_text_to_seqrecords(sample_data_text)
            return cls(ref_rec, sample_recs)

# --- End of Placeholder ---


class AlignmentThread(QThread):
    finished = pyqtSignal(object, object, str)  # (original_matrix_df, display_df, clustal_output)
    error = pyqtSignal(str, str)  # (error_type, message)

    def __init__(self, input_type: str, ref_data: str, sample_data: str, parent=None):
        super().__init__(parent)
        self.input_type = input_type
        self.ref_data = ref_data
        self.sample_data = sample_data

    def run(self):
        temp_excel_path = "temp_AlignmentSheet.xlsx" # Define path for temp file
        original_matrix_df = None # To store the matrix from comparator

        try:
            if self.input_type == 'files':
                sample_paths = [p.strip() for p in self.sample_data.split(";") if p.strip()]
                if not sample_paths:
                    raise InputFileError("No valid sample file paths provided.")
                comparator = SequenceComparator.from_files(self.ref_data, sample_paths)
            else: # 'text' input
                comparator = SequenceComparator.from_text(self.ref_data, self.sample_data)

            clustal_output = comparator.raw_clustal
            original_matrix_df = comparator.nucleotide_matrix # This is a DataFrame

            # export_alignment formats original_matrix_df for Excel display and returns it as display_df
            # It also writes to temp_excel_path as a side effect.
            display_df = export_alignment(original_matrix_df, temp_excel_path)
            
            if display_df is None:
                # This means export_alignment failed, which should be logged by export_alignment itself.
                # emit an error or handle this appropriately.
                # For now, assume if display_df is None, it's an error handled by on_alignment_finished checking df.
                 raise MafftError("Failed to generate the formatted alignment table for display (export_alignment returned None).")


            self.finished.emit(original_matrix_df, display_df, clustal_output)

        except (InputFileError, MafftExecutionError, MafftError) as e:
            self.error.emit(type(e).__name__, str(e))
        except Exception as e:
            logging.exception("Unexpected error in AlignmentThread")
            self.error.emit("UnknownError", f"Unexpected error: {str(e)}")
        finally:
            # Ensure the temporary file is deleted if it exists
            if os.path.exists(temp_excel_path):
                try:
                    os.remove(temp_excel_path)
                    logging.info(f"Temporary file {temp_excel_path} deleted.")
                except OSError as e_remove:
                    # Log error if deletion fails, but don't let it crash the thread
                    logging.error(f"Error deleting temporary file {temp_excel_path}: {e_remove}")


class LogHandler(logging.Handler):
    def __init__(self, text_edit):
        super().__init__()
        self.text_edit = text_edit

    def emit(self, record):
        msg = self.format(record)
        self.text_edit.appendPlainText(msg)


class FastaComparatorApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.alignment_thread = None
        self.current_data = {
            'original_matrix_df': None, # Will restore raw nucleotide matrix data 
            'display_df': None,         # Will store the DataFrame formatted for Excel/QTableWidget
            'clustal': ""
        }
        self.init_ui()
        self.setup_logging()
        self._on_input_changed()
        self.toggle_processing_ui(True)

    def init_ui(self):
        self.setWindowTitle("SNPsnap")
        self.setGeometry(100, 100, 900, 700)

        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        layout = QVBoxLayout(main_widget)

        self.input_tabs = QTabWidget()
        self._setup_file_input_tab()
        self._setup_text_input_tab()
        layout.addWidget(self.input_tabs)
        self.input_tabs.currentChanged.connect(self._on_input_changed)

        self._setup_control_buttons(layout)
        self._setup_output_tabs(layout)

    def _setup_file_input_tab(self):
        file_tab = QWidget()
        layout = QVBoxLayout(file_tab)
        ref_layout = QHBoxLayout()
        self.ref_edit = QLineEdit()
        self.ref_edit.setPlaceholderText("Path to reference FASTA file")
        btn_browse_ref = QPushButton("Browse...")
        btn_browse_ref.clicked.connect(lambda: self._browse_file(self.ref_edit))
        ref_layout.addWidget(QLabel("Reference FASTA:"))
        ref_layout.addWidget(self.ref_edit)
        ref_layout.addWidget(btn_browse_ref)
        sample_layout = QHBoxLayout()
        self.sample_edit = QLineEdit()
        self.sample_edit.setPlaceholderText("Path(s) to sample FASTA files, separated by ;")
        btn_browse_samples = QPushButton("Browse...")
        btn_browse_samples.clicked.connect(self._browse_samples)
        sample_layout.addWidget(QLabel("Sample FASTA(s):"))
        sample_layout.addWidget(self.sample_edit)
        sample_layout.addWidget(btn_browse_samples)
        layout.addLayout(ref_layout)
        layout.addLayout(sample_layout)
        self.input_tabs.addTab(file_tab, "File Input")
        self.ref_edit.textChanged.connect(self._on_input_changed)
        self.sample_edit.textChanged.connect(self._on_input_changed)

    def _setup_text_input_tab(self):
        text_tab = QWidget()
        layout = QVBoxLayout(text_tab)
        self.ref_textedit = QTextEdit()
        self.ref_textedit.setPlaceholderText("Paste reference FASTA content here...")
        self.samples_textedit = QTextEdit()
        self.samples_textedit.setPlaceholderText("Paste one or more sample FASTA records here...")
        layout.addWidget(QLabel("Reference FASTA:"))
        layout.addWidget(self.ref_textedit, stretch=1)
        layout.addWidget(QLabel("Sample FASTA(s):"))
        layout.addWidget(self.samples_textedit, stretch=1)
        self.input_tabs.addTab(text_tab, "Text Input")
        self.ref_textedit.textChanged.connect(self._on_input_changed)
        self.samples_textedit.textChanged.connect(self._on_input_changed)

    def _setup_control_buttons(self, layout):
        btn_layout = QHBoxLayout()
        self.btn_run = QPushButton("Run Alignment")
        self.btn_run.clicked.connect(self.run_alignment)
        self.btn_save_excel = QPushButton("Save Excel")
        self.btn_save_excel.clicked.connect(lambda: self.save_results('Excel'))
        self.btn_save_clustal = QPushButton("Save Clustal")
        self.btn_save_clustal.clicked.connect(lambda: self.save_results('Clustal'))
        btn_layout.addWidget(self.btn_run)
        btn_layout.addStretch(1)
        btn_layout.addWidget(self.btn_save_excel)
        btn_layout.addWidget(self.btn_save_clustal)
        layout.addLayout(btn_layout)

    def _setup_output_tabs(self, layout):
        self.output_tabs = QTabWidget()
        self.table_widget = QTableWidget()
        self.clustal_textedit = QTextEdit()
        self.clustal_textedit.setFont(QFont("Courier New", 10))
        self.clustal_textedit.setLineWrapMode(QTextEdit.LineWrapMode.NoWrap)
        self.clustal_textedit.setReadOnly(True)
        self.log_textedit = QPlainTextEdit()
        self.log_textedit.setReadOnly(True)
        self.log_textedit.setFont(QFont("Consolas", 9))
        self.log_textedit.setMaximumBlockCount(1000)
        self.output_tabs.addTab(self.table_widget, "Alignment Table")
        self.output_tabs.addTab(self.clustal_textedit, "Clustal Output")
        self.output_tabs.addTab(self.log_textedit, "Logs")
        layout.addWidget(self.output_tabs, stretch=1)

    def setup_logging(self):
        log_handler = LogHandler(self.log_textedit)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%H:%M:%S')
        log_handler.setFormatter(formatter)
        logger = logging.getLogger()
        logger.addHandler(log_handler)
        logger.setLevel(logging.INFO)
        logging.info("Application started. Logging initialized.")

    def _on_input_changed(self):
        if self.alignment_thread and self.alignment_thread.isRunning():
            self.btn_run.setEnabled(False)
        else:
            self.validate_inputs()

    def validate_inputs(self) -> bool:
        valid = False
        current_tab_index = self.input_tabs.currentIndex()
        try:
            if current_tab_index == 0:
                ref_path = self.ref_edit.text().strip()
                sample_paths_str = self.sample_edit.text().strip()
                if ref_path and sample_paths_str:
                    actual_sample_paths = [p.strip() for p in sample_paths_str.split(";") if p.strip()]
                    if os.path.exists(ref_path) and actual_sample_paths and \
                       all(os.path.exists(p) for p in actual_sample_paths):
                        valid = True
            else:
                if self.ref_textedit.toPlainText().strip() and \
                   self.samples_textedit.toPlainText().strip():
                    valid = True
        except Exception as e:
            logging.error(f"Error during input validation: {e}")
            valid = False
        self.btn_run.setEnabled(valid)
        return valid

    def toggle_processing_ui(self, processing_finished: bool):
        self.input_tabs.setEnabled(processing_finished)
        if processing_finished:
            self.btn_run.setText("Run Alignment")
            self._on_input_changed()
        else:
            self.btn_run.setText("Processing...")
            self.btn_run.setEnabled(False)
        # Enable/disable save buttons based on whether 'original_matrix_df' (for Excel) or 'clustal' is available
        can_save_excel = processing_finished and self.current_data.get('original_matrix_df') is not None
        can_save_clustal = processing_finished and bool(self.current_data.get('clustal'))
        self.btn_save_excel.setEnabled(can_save_excel)
        self.btn_save_clustal.setEnabled(can_save_clustal)

    def _browse_file(self, line_edit_widget):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select FASTA File", "", "FASTA Files (*.fasta *.fa *.fna);;All Files (*)"
        )
        if path:
            line_edit_widget.setText(path)
            logging.info(f"Selected file: {path}")
            self._on_input_changed()

    def _browse_samples(self):
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Select Sample FASTA Files", "", "FASTA Files (*.fasta *.fa *.fna);;All Files (*)"
        )
        if paths:
            self.sample_edit.setText(";".join(paths))
            logging.info(f"Selected {len(paths)} sample files: {';'.join(paths)}")
            self._on_input_changed()

    def run_alignment(self):
        if not self.validate_inputs():
            QMessageBox.warning(self, "Input Invalid", "Please ensure all inputs are valid before running.")
            return
        input_type = 'files' if self.input_tabs.currentIndex() == 0 else 'text'
        try:
            if input_type == 'files':
                ref_data = self.ref_edit.text().strip()
                sample_data = self.sample_edit.text().strip()
            else:
                ref_data = self.ref_textedit.toPlainText().strip()
                sample_data = self.samples_textedit.toPlainText().strip()

            self.table_widget.setRowCount(0)
            self.table_widget.setColumnCount(0)
            self.clustal_textedit.clear()
            self.current_data = {'original_matrix_df': None, 'display_df': None, 'clustal': ""}

            self.alignment_thread = AlignmentThread(input_type, ref_data, sample_data, parent=self)
            self.alignment_thread.finished.connect(self.on_alignment_finished)
            self.alignment_thread.error.connect(self.on_alignment_error)
            self.alignment_thread.start()
            self.toggle_processing_ui(False)
            logging.info(f"Alignment process started (type: {input_type})...")
        except Exception as e:
            logging.exception("Error preparing for alignment")
            self.show_error_message("SetupError", str(e))

    # MODIFIED: on_alignment_finished now receives original_matrix_df and display_df
    def on_alignment_finished(self, original_matrix_df, display_df, clustal_output):
        self.current_data.update({
            'original_matrix_df': original_matrix_df, # Store the raw matrix from comparator
            'display_df': display_df,                 # Store the formatted DataFrame for display
            'clustal': clustal_output
        })
        try:
            # Use display_df for the QTableWidget
            if display_df is not None and not display_df.empty:
                self.table_widget.setRowCount(display_df.shape[0])
                self.table_widget.setColumnCount(display_df.shape[1])
                self.table_widget.setHorizontalHeaderLabels(display_df.columns.astype(str).tolist())
                for i in range(display_df.shape[0]):
                    for j in range(display_df.shape[1]):
                        item_value = str(display_df.iat[i, j])
                        self.table_widget.setItem(i, j, QTableWidgetItem(item_value))
                self.table_widget.resizeColumnsToContents()
            else:
                logging.warning("Alignment finished but display DataFrame is None or empty.")
                self.table_widget.setRowCount(0)
                self.table_widget.setColumnCount(0)

            self.clustal_textedit.setText(clustal_output)
            logging.info("Alignment completed successfully. Results displayed.")
        except Exception as e:
            logging.exception("Error displaying alignment results")
            self.show_error_message("DisplayError", f"Error displaying results: {str(e)}")
        finally:
            self.toggle_processing_ui(True)

    def on_alignment_error(self, error_type, message):
        self.show_error_message(error_type, message)
        self.toggle_processing_ui(True)

    def show_error_message(self, error_type, message):
        error_map = {
            'InputFileError': 'Input Error',
            'MafftExecutionError': 'Alignment Execution Error',
            'MafftError': 'MAFFT Logic Error',
            'UnknownError': 'Unexpected Error',
            'SetupError': 'Application Setup Error',
            'DisplayError': 'Results Display Error'
        }
        title = error_map.get(error_type, 'Error')
        logging.error(f"{title} ({error_type}): {message}")
        QMessageBox.critical(self, title, f"{message}")

    def save_results(self, file_type: str):
        data_to_save = None
        default_suffix = ""
        file_filter = ""

        if file_type == 'Excel':
            # Use the original_matrix_df for saving, as export_alignment will re-process it.
            data_to_save = self.current_data.get('original_matrix_df')
            default_suffix = ".xlsx"
            file_filter = f"Excel Files (*{default_suffix})"
        elif file_type == 'Clustal':
            data_to_save = self.current_data.get('clustal')
            default_suffix = ".aln"
            file_filter = f"Clustal Alignment Files (*.aln);;Text Files (*.txt)"
        else:
            logging.error(f"Unknown file type for saving: {file_type}")
            QMessageBox.warning(self, "Save Error", f"Unknown file type: {file_type}")
            return

        if data_to_save is None:
            QMessageBox.warning(self, "No Data", f"No {file_type} data available to save.")
            return

        default_filename = f"AlignmentSheet{default_suffix}"
        initial_selected_filter = file_filter.split(";;")[0]
        path, _ = QFileDialog.getSaveFileName(
            self, f"Save {file_type} File", default_filename, file_filter, initial_selected_filter
        )

        if not path:
            return
        
        try:
            if file_type == 'Excel':
                # export_alignment takes the original_matrix_df and the target path,
                # then processes it into the Excel format and saves it.
                returned_df = export_alignment(data_to_save, path)
                if returned_df is None: # Check if export_alignment failed
                    raise Exception("export_alignment function failed to save the Excel file.")
            else:  # Clustal (text)
                with open(path, 'w') as f:
                    f.write(data_to_save)
            
            logging.info(f"Saved {file_type} file: {path}")
            QMessageBox.information(self, "Success", f"{file_type} file saved successfully to\n{path}")
        except Exception as e:
            logging.exception(f"Failed to save {file_type} file")
            QMessageBox.critical(self, "Save Error", f"Failed to save {file_type} file to\n{path}\n\nError: {str(e)}")

    def closeEvent(self, event):
        if self.alignment_thread and self.alignment_thread.isRunning():
            logging.info("Attempting to stop alignment thread on close...")
            if not self.alignment_thread.wait(1000):
                logging.warning("Alignment thread did not finish, terminating.")
                self.alignment_thread.terminate()
                self.alignment_thread.wait()
        super().closeEvent(event)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    app.setApplicationName("SNPsnap")
    app.setOrganizationName("MyOrg")
    window = FastaComparatorApp()
    window.show()
    sys.exit(app.exec())