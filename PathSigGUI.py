"""
PathSigGUI.py

"""
import os
import shutil
import sys
import datetime
from PySide6 import QtWidgets, QtGui, QtCore
from PySide6.QtWidgets import QApplication
import NumberCruncher
from libraries.StyleSheetLoader import load_stylesheet
import libraries.Tool_Box as ToolBox

__version__ = "0.0.2"
__author__ = "Dennis A. Simpson"
__copyright__ = "Copyright 2027, University of North Carolina at Chapel Hill"
__email__ = "dennis@email.unc.edu"

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.notice_color = "#8abaed"
        self.warning_color = "#b0ed8a"
        self.critical_error_color = "#f73e8e"

        self.setWindowTitle("PathSigGUI v" + __version__)
        self.setStyleSheet(load_stylesheet())

        self.path_to_tsv = None
        self.output_folder = None
        self.temp_tsv_path = "C:{0}Users{0}{1}{0}Documents{0}TempTSV.tsv".format(os.sep, os.getlogin())
        self.console = QtWidgets.QTextEdit()

        self.input_edit = QtWidgets.QLineEdit()
        self.input_btn = QtWidgets.QPushButton("Browse Input")
        self.input_label = QtWidgets.QLabel("TSV Data File:")

        self.output_folder_edit = QtWidgets.QLineEdit()
        self.output_btn = QtWidgets.QPushButton("Browse Output")
        self.output_label = QtWidgets.QLabel("Output folder:")

        self.fam_label = QtWidgets.QLabel("FAM Label:")
        self.fam_label_edit = QtWidgets.QLineEdit()
        self.fam_label_edit.setPlaceholderText("e.g., NHEJ, TMEJ-AAG, etc.")

        self.hex_label = QtWidgets.QLabel("HEX Label:")
        self.hex_label_edit = QtWidgets.QLineEdit()
        self.hex_label_edit.setPlaceholderText("e.g., NHEJ, TMEJ-AAG, etc.")

        self.cy5_label = QtWidgets.QLabel("CY5 Label:")
        self.cy5_label_edit = QtWidgets.QLineEdit()
        self.cy5_label_edit.setPlaceholderText("e.g., NHEJ, TMEJ-AAG, etc.")

        self.cy55_label = QtWidgets.QLabel("CY5.5 Label:")
        self.cy55_label_edit = QtWidgets.QLineEdit()
        self.cy55_label_edit.setPlaceholderText("e.g., NHEJ, TMEJ-AAG, etc.")

        self.rox_label = QtWidgets.QLabel("Rox Label:")
        self.rox_label_edit = QtWidgets.QLineEdit()
        self.rox_label_edit.setPlaceholderText("e.g., NHEJ, TMEJ-AAG, etc.")

        self.atto590_label = QtWidgets.QLabel("Atto590 Label:")
        self.atto590_label_edit = QtWidgets.QLineEdit()
        self.atto590_label_edit.setPlaceholderText("e.g., NHEJ, TMEJ-AAG, Copies, etc.")

        self.run_btn = QtWidgets.QPushButton("▶ Process Data")
        self.exit_btn = QtWidgets.QPushButton("✖ EXIT")
        self.main_widget = QtWidgets.QWidget()

        self.message_box = QtWidgets.QMessageBox()

        self.initialize_ui()

    def initialize_ui(self):
        # main_widget = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(self.main_widget)
        layout.setContentsMargins(20, 20, 20, 20)
        layout.setSpacing(15)

        self.setCentralWidget(self.main_widget)
        self.center_window(self.main_widget)

        # Header/Console Label
        self.console.setReadOnly(True)
        self.console.setPlaceholderText("System logs will appear here...")
        layout.addWidget(self.console)

        # Input path
        self.input_edit.setPlaceholderText("Select Input TSV File")
        self.input_btn.setObjectName("browse_btn")
        self.input_btn.clicked.connect(self.select_tsv_file)

        # Output path
        self.output_folder_edit.setPlaceholderText("Select output directory… (will default to TSV file folder if left blank)")
        self.output_btn.setObjectName("browse_btn")
        self.output_btn.clicked.connect(self.select_output_folder)

        label_width = 150

        for label in [self.input_label, self.output_label]:
            label.setFixedWidth(label_width)
            label.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter)

        grid = QtWidgets.QGridLayout()
        grid.addWidget(self.input_label, 0, 0)
        grid.addWidget(self.input_edit, 0, 1)
        grid.addWidget(self.input_btn, 0, 2)

        grid.addWidget(self.output_label, 1, 0)
        grid.addWidget(self.output_folder_edit, 1, 1)
        grid.addWidget(self.output_btn, 1, 2)

        grid.setColumnStretch(1, 1)
        layout.addLayout(grid)

        # ---- Channel Settings ----
        channel_box = QtWidgets.QGroupBox("Channel Settings")
        channel_layout = QtWidgets.QGridLayout()

        channel_label_list = \
            [self.fam_label, self.hex_label, self.cy5_label, self.cy55_label, self.rox_label, self.atto590_label]

        for label in channel_label_list:
            label.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter)
            label.setFixedWidth(label_width)

        channel_layout.addWidget(self.fam_label, 0, 0); channel_layout.addWidget(self.fam_label_edit, 0, 1)
        channel_layout.addWidget(self.hex_label, 1, 0); channel_layout.addWidget(self.hex_label_edit, 1, 1)
        channel_layout.addWidget(self.cy5_label, 0, 2); channel_layout.addWidget(self.cy5_label_edit, 1, 3)
        channel_layout.addWidget(self.cy55_label, 1, 2); channel_layout.addWidget(self.cy55_label_edit, 0, 3)
        channel_layout.addWidget(self.rox_label, 0, 4); channel_layout.addWidget(self.rox_label_edit, 0, 5)
        channel_layout.addWidget(self.atto590_label, 1, 4); channel_layout.addWidget(self.atto590_label_edit, 1, 5)

        # ---- Settings container (horizontal) ----
        settings_container = QtWidgets.QHBoxLayout()
        channel_box.setLayout(channel_layout)
        settings_container.addWidget(channel_box, 3)

        # Add the horizontal settings container to the main layout
        layout.addLayout(settings_container)
        channel_box.setLayout(channel_layout)

        # Button Row
        button_layout = QtWidgets.QHBoxLayout()
        button_layout.setSpacing(10)

        self.run_btn.setObjectName("run_btn")
        self.run_btn.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        self.run_btn.clicked.connect(self.process_data)

        self.exit_btn.setObjectName("exit_btn")
        self.exit_btn.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        self.exit_btn.clicked.connect(self.close)

        button_layout.addWidget(self.run_btn)
        button_layout.addWidget(self.exit_btn)

        layout.addLayout(button_layout)
        self.setCentralWidget(self.main_widget)



        #QtWidgets.QMessageBox.critical(self, "You Forgot Something", "Enter a valid TSV file and output folder.")
    def process_data(self):
        self.log_to_console(
            "Passing Channel Labels.  If any missing labels are intentional, it is safe to ignore the warnings",
            self.notice_color)

        channel_label_list = []
        if self.fam_label_edit.text():
            channel_label_list.append(self.fam_label_edit.text())
        else:
            channel_label_list.append("")
            self.log_to_console("No FAM Channel Label Entered.", self.warning_color)

        if self.hex_label_edit.text():
            channel_label_list.append(self.hex_label_edit.text())
        else:
            channel_label_list.append("")
            self.log_to_console("No HEX Channel Label Entered.", self.warning_color)

        if self.cy5_label_edit.text():
            channel_label_list.append(self.cy5_label_edit.text())
        else:
            channel_label_list.append("")
            self.log_to_console("No CY5 Channel Label Entered.", self.warning_color)

        if self.cy55_label_edit.text():
            channel_label_list.append(self.cy55_label_edit.text())
        else:
            channel_label_list.append("")
            self.log_to_console("No CY5.5 Channel Label Entered.", self.warning_color)

        if self.rox_label_edit.text():
            channel_label_list.append(self.rox_label_edit.text())
        else:
            channel_label_list.append("")
            self.log_to_console("No ROX Channel Label Entered.", self.warning_color)

        if self.atto590_label_edit.text():
            channel_label_list.append(self.atto590_label_edit.text())
        else:
            channel_label_list.append("")
            self.log_to_console("No ATTO590 Channel Label Entered.", self.warning_color)

        # Question popup with Yes/No buttons
        reply = self.message_box.question(
            self,
            "Take Heed",
            "Are Channel Labels Correct?",
            QtWidgets.QMessageBox.StandardButton.Yes | QtWidgets.QMessageBox.StandardButton.No,
            QtWidgets.QMessageBox.StandardButton.No  # Default focused button
        )

        if reply == QtWidgets.QMessageBox.StandardButton.Yes:
            if self.path_to_tsv and self.output_folder:
                NumberCruncher.DataProcessing(self.path_to_tsv, self.output_folder, channel_label_list)
            else:
                self.log_to_console("Please select a TSV file and/or output folder.", self.warning_color)
                self.message_box.critical(self, "You Forgot Something", "Enter a valid TSV file and output folder.")

    def select_tsv_file(self) -> None:

        self.path_to_tsv, _ = \
            QtWidgets.QFileDialog.getOpenFileName(self, self.tr("File Select"),
                                                  self.tr("C:{0}Users{0}{1}{0}Documents{0}"
                                                          .format(os.sep, os.getlogin())))

        if  os.path.basename(self.path_to_tsv) == "TempTSV.tsv":
            self.log_to_console("Cannot use file name TempTSV.tsv as input file.", self.warning_color)
        elif self.path_to_tsv:
            shutil.copyfile(self.path_to_tsv, self.temp_tsv_path)
        else:
            self.log_to_console("TSV File Not Selected.", self.warning_color)

        if self.path_to_tsv:
            self.input_edit.setText(self.path_to_tsv)

        if not self.output_folder:
            self.output_folder = os.path.dirname(self.path_to_tsv) + os.sep
            self.output_folder_edit.setText(self.output_folder)

    def select_output_folder(self):
        output_folder = QtWidgets.QFileDialog.getExistingDirectory(self, "Select Output Folder")

        if output_folder:
            self.output_folder_edit.setText(output_folder+os.sep)

    def log_to_console(self, text, color=None):
        if not color:
            color = self.notice_color
        self.console.moveCursor(QtGui.QTextCursor.MoveOperation.End)
        # always precede printed line with datetime stamp
        self.console.insertHtml(f"<span style='color: {color};'>[{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {text}</span><br>")
        self.console.ensureCursorVisible()

    def center_window(self, central_widget):
        """
        These settings are a bit empirical.  They are designed to center the window on the screen and set the size to make
        everything readable.

        :param central_widget:
        :return:
        """
        # screen = self.screen()
        screen_geometry = QtGui.QGuiApplication.primaryScreen().availableGeometry()
        central_widget.setWindowTitle("PathSig Data Processing")
        w_scale = 0.55
        h_scale = 0.85

        if screen_geometry.width() > 2000:
            w_scale = 0.70
        elif screen_geometry.width() < 1500:
            w_scale = 0.79
        if screen_geometry.height() > 1200:
            h_scale = 0.70
        # print(screen_geometry.width(), screen_geometry.height())
        width = int(screen_geometry.width() * w_scale)
        height = int(screen_geometry.height() * h_scale)
        self.resize(width, height)
        """
        central_widget.resize(screen_geometry.width() * w_scale, screen_geometry.height() * h_scale)
        central_widget.setGeometry(QtWidgets.QStyle.alignedRect(QtCore.Qt.LayoutDirection.LeftToRight,
                                                                QtCore.Qt.AlignmentFlag.AlignCenter,
                                                                central_widget.size(), screen_geometry, ), )
        """


if __name__ == "__main__":

    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    window = MainWindow()
    window.show()
    sys.exit(app.exec())