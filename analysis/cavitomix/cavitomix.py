# -*- coding: utf-8 -*-
"""
CavitOmiX
main class definitions of the PyMOL plugin
"""
import bz2
import os
import pickle
import re
from datetime import datetime
from logging import getLogger
from typing import Dict

from PyQt5.QtCore import QSize, Qt, QThreadPool
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QDialog,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QRadioButton,
    QSizePolicy,
    QSpacerItem,
    QTabWidget,
    QTextBrowser,
    QVBoxLayout,
    QWidget,
)

from .cavfind import CavFind, default_settings, param_types, tool_tips
from .mol_viewer import (
    create_object_from_pdbstr,
    display_cavities,
    get_coordinates,
    get_object_list,
    get_pdbstr,
    is_valid_object,
)
from .params import PluginParams
from .settings_gui import AdvancedSettings
from .structure_retrieval import retrieve_structures
from .thread import Worker
from .util import get_unique_name

PLUGIN_DIRECTORY = os.path.dirname(__file__)

logger = getLogger(__name__)

# threadpool for asynchronous execution
threadpool = QThreadPool()

# container for CavitOmiX objects
cavitomix_objs: Dict[str, CavFind] = {}

# global settings
global_settings = default_settings.copy()
settings_tooltips = tool_tips
settings_types = param_types


class MainWindow(QDialog):
    """
    class definition for the main window of the plugin
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.setWindowTitle("CavitOmiX")
        self.setLayout(QVBoxLayout())

        self.banner = Banner()  # banner with logo

        # widgets for the individual tabs
        self.struct_tab = StructWidget(self)
        self.struct_tab.setObjectName("struct_tab")
        self.run_tab = RunWidget(self)
        self.run_tab.setObjectName("run_tab")
        self.analyze_tab = AnalyzeWidget(self)
        self.analyze_tab.setObjectName("analyze_tab")
        self.settings_tab = SettingsWidget(self)
        self.settings_tab.setObjectName("settings_tab")
        self.log_tab = LogWidget(self)
        self.log_tab.setObjectName("log_tab")
        self.about_tab = AboutWidget(self)
        self.about_tab.setObjectName("about_tab")

        tabs = [
            (self.struct_tab, "&Retrieve Structures", "Retrieve protein structures."),
            (self.run_tab, "&Detect Cavities", "Setup and run cavity detection."),
            (self.analyze_tab, "Analyze &Cavities", "Analyze and display cavities."),
            (self.settings_tab, "&Settings", "Change global or object settings."),
            (self.log_tab, "&Log", "View logs."),
            (self.about_tab, "&About", "About CavitOmiX, credits,..."),
        ]

        self.tab_widget = QTabWidget(self)
        self.tab_widget.setTabPosition(QTabWidget.South)
        for i, (widget, tab_name, tab_tooltip) in enumerate(tabs):
            self.tab_widget.addTab(widget, tab_name)
            self.tab_widget.setTabToolTip(i, tab_tooltip)

        # assemble widgets
        self.layout().addLayout(self.banner)
        self.layout().addWidget(self.tab_widget)

        self.initialize_ui()
        self.setup_callbacks()

    def initialize_ui(self):
        """
        set initial state of the widget
        """
        self.resize(650, 420)  # set window size
        self.setMinimumSize(650, 420)

        # set the current tab to "Retrieve Structures"
        self.tab_widget.setCurrentWidget(self.struct_tab)

    def setup_callbacks(self):
        """
        register callback functions of the widget
        """
        self.tab_widget.currentChanged.connect(self.on_tab_change)

    def on_tab_change(self, index):
        if index == self.tab_widget.indexOf(self.run_tab):
            self.run_tab.refresh_object_list()
            self.run_tab.run_btn.setDefault(True)
        elif index == self.tab_widget.indexOf(self.analyze_tab):
            self.analyze_tab.show_btn.setDefault(True)


class Banner(QHBoxLayout):
    """
    definition of the banner
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # CavitOmiX logo
        logo = QLabel()
        logo_file = os.path.join(PLUGIN_DIRECTORY, "assets/cavitomix_logo.png")
        logo.setPixmap(QPixmap(logo_file))
        logo.setScaledContents(True)

        # add items to layout
        self.addSpacerItem(get_horizontal_spacer())
        self.addWidget(logo)


class StructWidget(QWidget):
    """
    definition of the structure retrieval widget
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setLayout(QHBoxLayout())

        left_column = QVBoxLayout()
        label_1 = QLabel(
            text="Enter PDB-codes, Uniprot accession numbers, or amino acid sequences (FASTA format)"
        )
        label_1.setWordWrap(True)
        self.seq_edit = QPlainTextEdit()
        left_column.addWidget(label_1)
        left_column.addWidget(self.seq_edit)

        right_column = QVBoxLayout()
        self.retrieve_btn = Button(text="Retrieve Structure(s)")
        self.align_chk = QCheckBox(text="Align Structures")
        self.align_chk.setToolTip(
            "Align all structures to the first structures by superimposing Calpha atoms."
        )
        label_2 = QLabel(text="Server to use for structure prediction:")
        label_2.setWordWrap(True)
        self.bionemo_btn = QRadioButton(text="BioNeMo\npowered by Nvidia")
        self.esmfold_btn = QRadioButton(text="ESMfold\npowered by Meta")
        self.pymolfold_btn = QRadioButton(text="PyMolFold")
        right_column.addSpacerItem(get_vertical_spacer())
        right_column.addWidget(self.retrieve_btn)
        right_column.addWidget(self.align_chk)
        right_column.addSpacerItem(get_vertical_spacer())
        right_column.addWidget(label_2)
        right_column.addWidget(self.bionemo_btn)
        right_column.addWidget(self.esmfold_btn)
        right_column.addWidget(self.pymolfold_btn)
        right_column.addSpacerItem(get_vertical_spacer())

        self.layout().addLayout(left_column)
        self.layout().addLayout(right_column)

        self.initialize_ui()
        self.setup_callbacks()

    def initialize_ui(self):
        """
        set initial state of the widget
        """
        self.retrieve_btn.setDefault(False)
        self.align_chk.setChecked(True)
        self.esmfold_btn.setChecked(True)

    def setup_callbacks(self):
        """
        register callback functions of the widget
        """
        self.retrieve_btn.clicked.connect(self.retrieve_models)

    def retrieve_models(self):
        """
        central dispatch of structure retrieval
        """

        query = self.seq_edit.toPlainText()

        if query:
            server = "BioNeMo"
            if self.esmfold_btn.isChecked():
                server = "ESMfold"
            elif self.pymolfold_btn.isChecked():
                server = "PyMolFold"

            if (server == "BioNeMo" and "bionemo" not in PluginParams.SERVER_LIST) or (
                server == "PyMolFold" and "pymolfold" not in PluginParams.SERVER_LIST
            ):
                logger.warning(
                    f"{server} server currently not available, using ESMfold instead!"
                )
                server = "ESMfold"

            align = self.align_chk.isChecked()

            self.retrieve_btn.setEnabled(False)
            self.seq_edit.setEnabled(False)
            # setup thread for structure retrieval
            _async_func = Worker(retrieve_structures, query, server, align)
            _async_func.signals.finished.connect(self._retrieve_finished)
            _async_func.signals.error.connect(_async_error)
            # _async_func.signals.result.connect(self._retrieve_result)
            # _async_func.signals.progress.connect(self._retrieve_progress)
            threadpool.start(_async_func)
        else:
            logger.warning("Query field contains no data!")
            _ = QMessageBox.warning(
                self,
                "Empty query!",
                "Query field contains no data!",
            )

    def _retrieve_finished(self):
        self.retrieve_btn.setEnabled(True)
        self.seq_edit.setEnabled(True)


class RunWidget(QWidget):
    """
    definition of the cavity detection widget
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setLayout(QHBoxLayout())

        left_column = QVBoxLayout()
        self.object_list = QListWidget(self)
        self.object_list.setToolTip("List of current PyMOL objects.")
        label = QLabel(text="Name of next CavitOmiX object:")
        self.name_edit = QLineEdit(self)
        left_column.addWidget(self.object_list)
        left_column.addWidget(label)
        left_column.addWidget(self.name_edit)

        right_column = QVBoxLayout()
        self.run_btn = Button(text="Run")
        self.run_btn.setToolTip("Run CavFind to detect cavities.")
        self.refresh_btn = Button(text="Refresh")
        self.refresh_btn.setToolTip("Update PyMOL object list.")
        self.settings_btn = Button(text="Settings")
        self.refresh_btn.setToolTip("Edit settings for the next run.")
        self.progress_label = QLabel()
        self.progress_label.setAlignment(Qt.AlignCenter)
        self.progress_label.setWordWrap(True)
        self.load_btn = Button(text="Load CavObj")
        self.load_btn.setToolTip("Load CavitOmiX object.")
        right_column.addWidget(self.run_btn)
        right_column.addWidget(self.refresh_btn)
        right_column.addWidget(self.settings_btn)
        right_column.addSpacerItem(get_vertical_spacer())
        right_column.addWidget(self.progress_label)
        right_column.addSpacerItem(get_vertical_spacer())
        right_column.addWidget(self.load_btn)

        self.layout().addLayout(left_column)
        self.layout().addLayout(right_column)

        self.initialize_ui()
        self.setup_callbacks()

    def initialize_ui(self):
        """
        set initial state of the widget
        """
        self.run_btn.setEnabled(False)
        self.name_edit.setEnabled(False)
        self.progress_label.setText("")

    def setup_callbacks(self):
        """
        register callback functions of the widget
        """
        self.object_list.currentItemChanged.connect(self.select_item)
        self.object_list.itemClicked.connect(self.select_item)
        self.run_btn.clicked.connect(self.run_cavity_detection)
        self.refresh_btn.clicked.connect(self.refresh_object_list)
        self.settings_btn.clicked.connect(self.edit_global_settings)
        self.load_btn.clicked.connect(self.load_cavitomix_object)

    def refresh_object_list(self):
        """
        fill 'object_list' with the current list of objects from PyMOL
        """
        self.object_list.clear()  # clear current list
        self.object_list.addItems(get_object_list())

        if self.object_list.count() > 0:  # set the current item to the first list item
            self.object_list.setCurrentItem(self.object_list.item(0))

    def select_item(self, item):
        """
        This function is called, when an item in the list is clicked/selected/made the 'current item'
        :param item: new current item or 'None', when the list was cleared
        """
        # necessary, because the signal 'currentItemChanged' is also emitted upon 'clear'
        if item is not None:
            self.run_btn.setEnabled(True)  # enable run button
            obj = self.object_list.currentItem().text()
            if obj.startswith("("):
                obj = obj[1:-1]
            name = get_unique_name(f"{obj}_cav", cavitomix_objs.keys())
            self.name_edit.setText(name)
            self.name_edit.setEnabled(True)
        else:
            self.run_btn.setEnabled(False)  # disable run button
            self.name_edit.clear()  # clear the name input field
            self.name_edit.setEnabled(False)  # disable the name input field

    def load_cavitomix_object(self):
        """
        This function is called, when 'load_btn' is clicked
        """
        filename = open_file_dialog(
            self,
            "Load CavitOmiX Object",
            "CavitOmiX Files (*.cavitomix);;All Files (*.*)",
        )

        if filename:
            with bz2.BZ2File(filename, "r") as infile:
                try:
                    cav_find = pickle.load(infile)
                except ModuleNotFoundError as e:
                    logger.error(f"Unable to load CavitOmiX object: {e}")
                    return
                except OSError as e:
                    logger.error(f"Unable to load CavitOmiX object: {e}")
                    return

            # display PDB-structure in PyMOL
            obj_name = cav_find.obj_name
            logger.info(
                f"Name of structure/selection used for cavity calculation: {obj_name}"
            )

            # if cavities were calculated from a selection
            if obj_name[0] == "(" and obj_name[-1] == ")":
                obj_name = obj_name[1:-1]  # remove parentheses
                if obj_name == "sele":  # reserved name for selection
                    obj_name = "STRUCT"

            obj_name = create_object_from_pdbstr(cav_find.pdb.get_pdbstr(), obj_name)
            cav_find.obj_name = obj_name
            logger.info(f"PyMOL object created: {obj_name}")
            self.refresh_object_list()

            name = get_unique_name(f"{obj_name}_cav", cavitomix_objs)
            logger.info(f"Name of CavitOmiX object: {name}")
            cavitomix_objs[name] = cav_find

            # update object list and settings box
            # references to other widgets
            tab_widget = self.parent().parent()
            analyze_tab = tab_widget.findChild(QWidget, "analyze_tab")
            settings_tab = tab_widget.findChild(QWidget, "settings_tab")

            analyze_tab.object_list.addItem(name)
            analyze_tab.object_list.setCurrentIndex(len(cavitomix_objs) - 1)
            settings_tab.settings_box.addItem(name)

            logger.info(f"CavitOmiX object '{name}' loaded from: {filename}")

    def edit_global_settings(self):
        """
        This function is called, when 'settings_btn' is clicked
        """
        tab_widget = self.parent().parent()
        settings_tab = tab_widget.findChild(QWidget, "settings_tab")
        settings_tab.settings_box.setCurrentIndex(0)
        tab_widget.setCurrentWidget(settings_tab)

    def run_cavity_detection(self):
        """
        Find cavities in selected object from 'object_list' and
        generate a new CavitOmiX object
        """
        obj = self.object_list.currentItem().text()
        if obj.startswith("("):
            obj = obj[1:-1]

        name = self.name_edit.text()

        if name in cavitomix_objs:
            reply = QMessageBox.question(
                self,
                "Object already exists",
                f"CavitOmiX object {name:s} already exists and will be overwritten.\n"
                "Do you want to continue?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No,
            )
            if reply == QMessageBox.No:
                return

        if not is_valid_object(obj):
            logger.warning(f"Object '{obj}' contains no atoms fulfilling all criteria!")
            _ = QMessageBox.warning(
                self,
                "No atoms in object",
                f"Object {obj} contains no atoms\nfulfilling all criteria!",
            )
            return

        logger.info(f"Object '{obj}' selected for LigSite calculation")

        # generate CavFind object
        cav_find = CavFind(name=name, obj_name=obj, settings=global_settings)
        # add structure from PyMOL 'pdbstr'
        cav_find.struct_from_pdb(get_pdbstr(obj))
        n_atoms = len(cav_find.pdb.atom)
        logger.info(f"Number of atoms fulfilling all criteria: {n_atoms}")

        cavitomix_objs[name] = cav_find

        # find unique name for the next object
        next_name = get_unique_name(f"{obj}_cav", cavitomix_objs.keys())
        self.name_edit.setText(next_name)

        # setup thread for LigSite calculation
        _async_func = Worker(cav_find.run, rerun=False)
        _async_func.signals.result.connect(self._ligsite_result)
        _async_func.signals.finished.connect(self._ligsite_finished)
        _async_func.signals.error.connect(_async_error)
        _async_func.signals.progress.connect(self._ligsite_progress)

        time = datetime.now()
        logger.info(
            f"LigSite calculation started: {time.strftime('%Y/%m/%d, %H:%M:%S')}"
        )
        logger.info("++ Parameters:")
        for key, value in cavitomix_objs[name].settings.items():
            if key not in ["split_files", "shape_dmax", "shape_object", "cp_limit"]:
                logger.info(f"++ {key}: {value}")
        logger.info("++")

        self.run_btn.setEnabled(False)
        threadpool.start(_async_func)

    def _ligsite_result(self, cav_find):
        name = cav_find.name
        self.progress_label.setText(f"{name}:\nFinished")
        logger.info(f"{name}: Finished")

        # references to other widgets
        tab_widget = self.parent().parent()
        analyze_tab = tab_widget.findChild(QWidget, "analyze_tab")
        settings_tab = tab_widget.findChild(QWidget, "settings_tab")

        # update object list and settings box
        index = analyze_tab.object_list.findText(name)
        if index == -1:  # new name
            analyze_tab.object_list.addItem(name)
            settings_tab.settings_box.addItem(name)
            analyze_tab.object_list.setCurrentIndex(len(cavitomix_objs) - 1)
        else:  # name already exists
            current_index = analyze_tab.object_list.currentIndex()
            analyze_tab.object_list.setCurrentIndex(index)
            settings_tab.settings_box.setCurrentIndex(index + 1)
            if index == current_index:  # no automatic refresh of cavity list
                analyze_tab.select_cavitomix_object(index)

    def _ligsite_finished(self):
        time = datetime.now()
        logger.info(
            f"LigSite calculation finished: {time.strftime('%Y/%m/%d, %H:%M:%S')}"
        )
        self.run_btn.setEnabled(True)

    def _ligsite_progress(self, value):
        self.progress_label.setText(str(value))


class AnalyzeWidget(QWidget):
    """
    definition of the cavity analysis widget
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setLayout(QVBoxLayout())

        top_row = QHBoxLayout()
        self.object_list = QComboBox(self)
        self.object_list.setSizePolicy(
            QSizePolicy(QSizePolicy.MinimumExpanding, QSizePolicy.Fixed)
        )
        self.object_list.setToolTip("Available CavitOmiX objects.")
        self.settings_btn = Button(text="Settings")
        self.settings_btn.setToolTip("Edit settings of selected CavitOmiX object.")
        self.saveobj_btn = Button(text="Save CavObj")
        self.saveobj_btn.setToolTip("Save selected CavitOmiX object.")
        top_row.addWidget(self.object_list)
        top_row.addWidget(self.settings_btn)
        top_row.addWidget(self.saveobj_btn)

        bottom_row = QHBoxLayout()
        left_column = QVBoxLayout()
        label_1 = QLabel(text="Cavities")
        label_1.setWordWrap(True)
        self.cavity_list = QListWidget(self)
        self.cavity_list.setSelectionMode(QAbstractItemView.MultiSelection)
        self.cavity_list.setToolTip("List of detected cavities.")
        left_column.addWidget(label_1)
        left_column.addWidget(self.cavity_list)

        middle_column = QVBoxLayout()
        label_2 = QLabel(text="Annotations")
        label_2.setWordWrap(True)
        self.annotation_list = QListWidget(self)
        self.annotation_list.setSelectionMode(QAbstractItemView.MultiSelection)
        self.annotation_list.setToolTip("Available annotations.")
        middle_column.addWidget(label_2)
        middle_column.addWidget(self.annotation_list)

        right_column = QVBoxLayout()
        self.savecav_btn = Button(text="Save Cavities")
        self.savecav_btn.setToolTip("Save cavities as PDB or CSV-files.")
        self.show_btn = Button(text="Display Cavities")
        self.show_btn.setToolTip("Display cavities in PyMOL window.")
        self.showenv_chk = QCheckBox(text="Show\nSurrounding\nResidues")
        self.showenv_chk.setToolTip("Display residues around cavities.")
        self.shape_btn = Button(text="Shape Cavities")
        self.shape_btn.setToolTip("Shape all cavities based on PyMOL object.")
        self.rerun_btn = Button(text="Rerun LigSite")
        self.rerun_btn.setToolTip("Rerun cavity detection in existing grid.")
        right_column.addWidget(self.savecav_btn)
        right_column.addWidget(self.show_btn)
        right_column.addWidget(self.showenv_chk)
        right_column.addSpacerItem(get_vertical_spacer())
        right_column.addWidget(self.shape_btn)
        right_column.addWidget(self.rerun_btn)

        bottom_row.addLayout(left_column)
        bottom_row.addLayout(middle_column)
        bottom_row.addLayout(right_column)

        self.layout().addLayout(top_row)
        self.layout().addLayout(bottom_row)

        self.initialize_ui()
        self.setup_callbacks()

    def initialize_ui(self):
        """
        set initial state of the widget
        """
        self.settings_btn.setEnabled(False)
        self.saveobj_btn.setEnabled(False)
        self.savecav_btn.setEnabled(False)
        self.show_btn.setEnabled(False)
        self.showenv_chk.setEnabled(False)
        self.shape_btn.setEnabled(False)
        self.rerun_btn.setEnabled(False)

    def setup_callbacks(self):
        """
        register callback functions of the widget
        """
        self.cavity_list.itemSelectionChanged.connect(self.select_cavity)
        self.object_list.currentIndexChanged.connect(self.select_cavitomix_object)
        self.settings_btn.clicked.connect(self.edit_local_settings)
        self.saveobj_btn.clicked.connect(self.save_cavitomix_object)
        self.savecav_btn.clicked.connect(self.save_cavities)
        self.show_btn.clicked.connect(self.show_cavities)
        self.shape_btn.clicked.connect(self.shape_cavities)
        self.rerun_btn.clicked.connect(self.rerun_cavity_detection)

    def select_cavity(self):
        """
        This function is called, when the selection in the cavity list is changed
        """
        sel_list = self.cavity_list.selectedItems()
        if len(sel_list) > 0:
            self.savecav_btn.setEnabled(True)
            self.show_btn.setEnabled(True)
            self.showenv_chk.setEnabled(True)
        else:
            self.savecav_btn.setEnabled(False)
            self.show_btn.setEnabled(False)
            self.showenv_chk.setEnabled(False)

    def select_cavitomix_object(self, index):
        """
        This function is called, when a CavitOmiX object is selected in the ComboBox
        :param index: new current index
        """
        if index is not None:
            obj = self.object_list.currentText()
            self.settings_btn.setEnabled(True)
            self.saveobj_btn.setEnabled(True)
            self.shape_btn.setEnabled(True)
            self.rerun_btn.setEnabled(True)

            # empty the lists
            self.cavity_list.clear()
            self.annotation_list.clear()

            # fill 'cavityList'
            if obj in cavitomix_objs:
                cav_list = cavitomix_objs[obj].cavities
                if cav_list:
                    item_list = [
                        f"#{i + 1}, {cav.size} pts., {cav.volume:.0f} Å^3"
                        for i, cav in enumerate(cav_list)
                    ]
                    self.cavity_list.addItems(item_list)

                    # fill 'annotationList'
                    # 'LIG' is not added to the list, because the LigSite score
                    # is mapped to B-factor column by default
                    # if no annotation is selected, only cavities annotated with the LigSite score
                    # are saved or displayed
                    item_list = [
                        PluginParams.annotation_names[key]
                        for key in cav_list[0].annotations.keys()
                        if key != "LIG"
                    ]
                    self.annotation_list.addItems(item_list)

    def edit_local_settings(self):
        """
        This function is called, when 'settings_btn' is clicked
        """
        tab_widget = self.parent().parent()
        settings_tab = tab_widget.findChild(QWidget, "settings_tab")
        settings_tab.settings_box.setCurrentIndex(self.object_list.currentIndex() + 1)
        tab_widget.setCurrentWidget(settings_tab)

    def save_cavitomix_object(self):
        """
        This function is called, when 'saveobj_btn' is clicked
        """
        name = self.object_list.currentText()
        filename = save_file_dialog(
            self, "Save CavitOmiX Object As...", "CavitOmiX File (*.cavitomix)"
        )

        if filename:
            obj = cavitomix_objs[name]
            with bz2.BZ2File(filename, "w") as outfile:
                pickle.dump(obj, outfile)
            logger.info(f"CavitOmiX object '{name}' written to: {filename}")

    def _get_annotation_list(self):
        return [
            PluginParams.annotation_keys[item.text()]
            for item in self.annotation_list.selectedItems()
        ]

    def save_cavities(self):
        """
        This function is called, when the 'savecav_btn' is clicked
        """
        cavity_indices = [
            self.cavity_list.row(item) for item in self.cavity_list.selectedItems()
        ]
        annotation_list = self._get_annotation_list()

        filename = save_file_dialog(
            self, "Save Cavities As...", "PDB File (*.pdb);;CSV File (*.csv)"
        )

        if filename:
            basename, extension = os.path.splitext(filename)
            file_format = extension[1:]  # get file_format from extension

            name = self.object_list.currentText()  # get the CavitOmiX object
            cav_find = cavitomix_objs[name]

            header = ""
            if file_format == "pdb" or file_format == "csv":
                tag = "REMARK" if file_format == "pdb" else "#"
                header = "\n".join(
                    [
                        f"{tag} CavitOmiX object: {cav_find.name}",
                        f"{tag} created: {cav_find.time_of_creation.strftime('%Y/%m/%d, %H:%M:%S')}",
                        "",
                    ]
                )

            cav_find.write_cavities(
                basename,
                file_format,
                cavity_indices,
                annotation=annotation_list,
                header=header,
            )

    def show_cavities(self):
        """
        This function is called, when 'show_btn' is clicked
        """
        obj = self.object_list.currentText()  # get the CavitOmiX object
        cavity_indices = [
            self.cavity_list.row(item) for item in self.cavity_list.selectedItems()
        ]
        annotation_list = self._get_annotation_list()

        display_cavities(
            cavitomix_objs[obj],
            cavity_indices,
            annotation_list,
            self.showenv_chk.isChecked(),
            4.5,
        )

    def shape_cavities(self):
        """
        shape all cavities of the active CavitOmiX object based on atoms in 'sele'
        cavity points are kept, if they are closer than 'shape_dmax' from any atom in 'sele'
        cavities are not removed from the list, if no points are left after shaping
        """
        obj = self.object_list.currentText()
        cav_find = cavitomix_objs[obj]

        shape_object = cav_find.settings["shape_object"]
        dmax = cav_find.settings["shape_dmax"]

        if not is_valid_object(shape_object):
            _ = QMessageBox.warning(
                self,
                "No or empty object",
                f"'{shape_object}' does not exist\nor is empty!",
            )
            return

        shape_coords = get_coordinates(shape_object)

        logger.info(
            f"Shaping cavities in '{obj}' using atoms in '{shape_object}' "
            f"with a max. distance of {dmax:.2f} Å."
        )

        shaped_cavities = []  # new list of shaped cavities
        for cav in cav_find.cavities:
            retained_points = cav.shape(shape_coords, dmax)
            if retained_points > 0:  # only keep non-empty cavities
                shaped_cavities.append(cav)

        cav_find.cavities = shaped_cavities

        self.select_cavitomix_object(self.object_list.currentIndex())

    def rerun_cavity_detection(self):
        """
        rerun cavity detection. e.g. with different 'ligsite_cutoff' or 'gap_parameter'
        old cavities are deleted
        """
        obj = self.object_list.currentText()

        logger.info(f"Object '{obj}' selected for re-calculation of cavities")

        # setup thread for LigSite calculation
        _async_func = Worker(cavitomix_objs[obj].run, rerun=True)
        _async_func.signals.result.connect(self._cav_detect_result)
        _async_func.signals.finished.connect(self._cav_detect_finished)
        _async_func.signals.error.connect(_async_error)
        # _async_func.signals.progress.connect(self._cav_detect_progress)

        time = datetime.now()
        logger.info(
            f"Cavity detection in existing grid started: {time.strftime('%Y/%m/%d, %H:%M:%S')}"
        )
        logger.info("++ Parameters:")
        for key, value in cavitomix_objs[obj].settings.items():
            if key in [
                "ligsite_cutoff",
                "gap",
                "min_size",
                "max_size",
                "original_ligsite",
                "max_dist",
                "annotate",
            ]:
                logger.info(f"++ {key}: {value}")
        logger.info("++")

        self.rerun_btn.setEnabled(False)
        self.shape_btn.setEnabled(False)
        self.settings_btn.setEnabled(False)
        self.saveobj_btn.setEnabled(False)
        self.cavity_list.clear()
        self.annotation_list.clear()
        threadpool.start(_async_func)

    def _cav_detect_result(self, cav_find):
        name = cav_find.name
        index = self.object_list.findText(name)
        if index >= 0:
            if index != self.object_list.currentIndex():
                self.object_list.setCurrentIndex(index)
            else:
                self.select_cavitomix_object(index)

    def _cav_detect_finished(self):
        time = datetime.now()
        logger.info(f"Cavity detection finished: {time.strftime('%Y/%m/%d, %H:%M:%S')}")
        self.rerun_btn.setEnabled(True)
        self.shape_btn.setEnabled(True)
        self.settings_btn.setEnabled(True)
        self.saveobj_btn.setEnabled(True)


class SettingsWidget(QWidget):
    """
    definition of the settings widget
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setLayout(QVBoxLayout())

        top_row = QHBoxLayout()
        self.settings_box = QComboBox()
        self.settings_box.setSizePolicy(
            QSizePolicy(QSizePolicy.MinimumExpanding, QSizePolicy.Fixed)
        )

        self.load_btn = Button(text="Load Settings")
        self.load_btn.setToolTip("Load settings from file.")
        self.save_btn = Button(text="Save Settings")
        self.save_btn.setToolTip("Save settings to file.")
        top_row.addWidget(self.settings_box)
        top_row.addWidget(self.load_btn)
        top_row.addWidget(self.save_btn)

        self.settings_table = AdvancedSettings(
            parent=None,
            settings=global_settings,
            tool_tips=settings_tooltips,
            param_types=settings_types,
            disabled_entries=[],
        )
        self.settings_table.setSizePolicy(
            QSizePolicy(QSizePolicy.MinimumExpanding, QSizePolicy.MinimumExpanding)
        )

        self.layout().addLayout(top_row)
        self.layout().addWidget(self.settings_table)

        self.initialize_ui()
        self.setup_callbacks()

    def initialize_ui(self):
        """
        set initial state of the widget
        """
        self.settings_box.addItem("global")

    def setup_callbacks(self):
        """
        register callback functions of the widget
        """
        self.settings_box.currentIndexChanged.connect(self.update_settings_table)
        self.save_btn.clicked.connect(self.save_settings)
        self.load_btn.clicked.connect(self.load_settings)

    def update_settings_table(self):
        name = self.settings_box.currentText()
        w = self.settings_table

        if name != "global":
            disabled_entries = [  # settings that cannot be changed at the object level
                "keep_hydrogens",  # PDB interpretation
                "keep_hetatms",
                "keep_waters",
                "resn_water",
                "alternate",
                "grid_spacing",  # Grid setup
                "probe_radius",
                "softness",
                "cushion",
                "radii",
                "radius_factor",
            ]
            w.settings = cavitomix_objs[name].settings
            w.disabled_entries = disabled_entries
        else:
            w.settings = global_settings
            w.disabled_entries = []
        w.model.clear()
        w.populateData()
        w.formatTable()

    def save_settings(self):
        """
        Save CavitOmiX settings: global or on object level
        """
        filename = save_file_dialog(
            self, "Save CavitOmiX Settings As...", "CavitOmiX Settings (*.cavset)"
        )

        if filename:
            obj = self.settings_box.currentText()  # what to save

            if obj == "global":  # global settings
                with bz2.BZ2File(filename, "w") as outfile:
                    pickle.dump(global_settings, outfile)
                logger.info(f"Global settings written to: {filename}")
            else:
                with bz2.BZ2File(filename, "w") as outfile:
                    pickle.dump(cavitomix_objs[obj].settings, outfile)
                logger.info(
                    f"Settings for CavitOmiX object '{obj}' written to: {filename}"
                )

    def load_settings(self):
        """
        Load CavitOmiX settings: global or on object level
        """
        filename = open_file_dialog(
            self,
            "Load CavitOmiX settings",
            "CavitOmiX Settings (*.cavset);;All Files (*.*)",
        )

        if filename:
            with bz2.BZ2File(filename, "r") as infile:
                new_settings = pickle.load(infile)

            obj = self.settings_box.currentText()
            if obj == "global":  # global settings
                global_settings.update(new_settings)
                logger.info(f"Global settings read from: {filename}")
            else:
                cavitomix_objs[obj].settings.update(new_settings)
                logger.info(
                    f"Settings for CavitOmiX object '{obj}' read from: {filename}"
                )

            self.update_settings_table()


class LogWidget(QWidget):
    """
    definition of the widget showing info-, warning and error-messages
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setLayout(QVBoxLayout())
        self.log_browser = QTextBrowser()
        self.layout().addWidget(self.log_browser)


class AboutWidget(QWidget):
    """
    definition of the widget showing the 'about'-message, credits, references,...
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setLayout(QVBoxLayout())

        text = (
            f"CavitOmiX v{PluginParams.__version__}\n"
            + f"{PluginParams.__copyright__}\n"
            + ", ".join(PluginParams.__developers__)
        )
        self.label = QLabel(text)
        self.label.setAlignment(Qt.AlignCenter)

        self.about_text = QTextBrowser()
        with open(os.path.join(PLUGIN_DIRECTORY, "assets/about.txt"), "r") as f:
            self.about_text.setText(
                f.read()
                .replace("__version__", f"{PluginParams.__version__}")
                .replace("__year__", f"{PluginParams.__year__}")
            )

        self.layout().addWidget(self.label)
        self.layout().addWidget(self.about_text)


class Button(QPushButton):
    """
    definition of a standard push button used in the CavitOmiX plugin
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setMinimumSize(QSize(150, 0))
        self.setAutoDefault(False)
        self.setSizePolicy(QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed))


def get_horizontal_spacer():
    return QSpacerItem(40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)


def get_vertical_spacer():
    return QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)


def _async_error(returned_data):
    exc_type, value, trace_back = returned_data
    logger.error(
        "\nAn exception was raised:\n  Type: {}\n  Text: {}\n".format(
            exc_type.__name__, value
        )
    )
    logger.error(trace_back)


def open_file_dialog(widget, dialog_title="Load File", format_filter="All files (*.*)"):
    # file dialog based on examples shown at https://pythonspot.com/pyqt5-file-dialog/
    options = QFileDialog.Options()
    options |= QFileDialog.DontUseNativeDialog
    file_name, _ = QFileDialog.getOpenFileName(
        widget, dialog_title, "", format_filter, options=options
    )

    return file_name


def save_file_dialog(
    widget, dialog_title="Save As...", format_filter="All files (*.*)"
):
    """
    Return a file name, append extension from filter if no extension provided.
    """
    # file dialog based on the function 'getSaveFileNameWithExt' from pymol.Qt.utils
    options = QFileDialog.Options()
    options |= QFileDialog.DontUseNativeDialog
    file_name, selected_filter = QFileDialog.getSaveFileName(
        widget, dialog_title, "", format_filter, options=options
    )

    if file_name:
        if "." not in os.path.split(file_name)[-1]:
            m = re.search(r"\*(\.[\w.]+)", selected_filter)
            if m:
                # append first extension from filter
                file_name += m.group(1)

        return file_name
