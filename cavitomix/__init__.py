# -*- coding: utf-8 -*-
"""
PyMOL Qt-Plugin

CavitOmiX App
"""
import logging
import os
import sys

logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)

dir_path = os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)
)
# append the directory of the CavitOmiX module
# this is necessary for pickled CavitOmiX objects to be readable in different PyMOL installations
# and different computers
# otherwise, the CavitOmiX class definitions are not found upon opening pickle-files
sys.path.append(dir_path)

# Avoid importing "expensive" modules here (e.g. scipy), since this code is
# executed on PyMOL's startup. Only import such modules inside functions.


def __init_plugin__(app=None):
    """
    Add an entry to the PyMOL 'Plugin' menu
    """
    from pymol.plugins import addmenuitemqt

    addmenuitemqt("CavitOmiX", command=run_plugin_gui, menuName="PluginQt")


# global reference to avoid garbage collection of our dialog
dialog = None


def run_plugin_gui():
    """
    Open our custom dialog
    """
    from .cavitomix import MainWindow
    from .structure_retrieval import contact_inno_server

    global dialog

    if dialog is None:  # if nonexistent, create window object
        contact_inno_server()
        dialog = MainWindow()
        handler = QTextEditLoggingHandler(dialog.log_tab.log_browser)
        handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)

    else:  # otherwise refresh the content
        dialog.run_tab.refresh_object_list()

    dialog.show()


class QTextEditLoggingHandler(logging.Handler):

    COLORS = {logging.WARNING: "orange", logging.INFO: "green", logging.DEBUG: "black"}

    def __init__(self, text_edit):
        logging.Handler.__init__(self)
        self.text_edit = text_edit

    def emit(self, record):
        color = self.COLORS.get(record.levelno, "red")
        text = self.format(record)
        self.text_edit.append("<font color={}>{}</font>".format(color, text))
