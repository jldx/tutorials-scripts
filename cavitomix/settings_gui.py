# -*- coding: utf-8 -*-
from logging import getLogger

from PyQt5.QtCore import QSortFilterProxyModel, Qt
from PyQt5.QtGui import QStandardItem, QStandardItemModel
from PyQt5.QtWidgets import QLineEdit, QTableView, QVBoxLayout, QWidget

logger = getLogger(__name__)


class AdvancedSettings(QWidget):
    def __init__(
        self,
        parent=None,
        settings=None,
        tool_tips=None,
        param_types=None,
        disabled_entries=None,
    ):
        super().__init__(parent)

        self.settings = settings or {}
        self.tool_tips = tool_tips or {}
        self.param_types = param_types or {}
        self.disabled_entries = disabled_entries or []

        self.model = QStandardItemModel(self)
        self.proxy_model = QSortFilterProxyModel(self)
        self.proxy_model.setSourceModel(self.model)

        layout = QVBoxLayout(self)
        self.setLayout(layout)
        self.filter_le = QLineEdit(self)
        layout.addWidget(self.filter_le)
        self.filter_le.setPlaceholderText("Filter")
        self.filter_le.textChanged.connect(self.proxy_model.setFilterRegExp)

        self.populateData()

        self.table = QTableView(self)
        self.table.setModel(self.proxy_model)
        layout.addWidget(self.table)

        self.formatTable()

        self.model.itemChanged.connect(self.itemChanged)

    def populateData(self):
        """
        Fill the model with data from dictionaries
        """
        for name, value in self.settings.items():
            value_item = QStandardItem()
            name_item = QStandardItem(name)

            name_item.setFlags(Qt.ItemIsEnabled)
            name_item.setToolTip(self.tool_tips.get(name, ""))

            param_type = self.param_types.get(name, None)
            if param_type in (0, 1, 2, 3, 4):  # if a type is available
                # types of setting parameters:
                # 0 ... boolean
                # 1 ... integer
                # 2 ... float
                # 3 ... str
                # 4 ... tuple
                if param_type == 0:  # boolean
                    value_item.setCheckable(True)
                    value_item.setEditable(False)  # Can't edit text (but toggles)
                    if value:
                        value_item.setCheckState(Qt.Checked)
                elif param_type in (1, 2, 3):  # int, float, or str
                    value_item.setText(str(value))
                elif param_type == 4:  # tuple
                    value_item.setText(", ".join(value))
            else:
                value_item.setText(str(value))  # unknown types are disabled
                name_item.setEnabled(False)
                value_item.setEnabled(False)

            if name in self.disabled_entries:
                name_item.setEnabled(False)
                value_item.setEnabled(False)

            self.model.appendRow([name_item, value_item])
            value_item.setData((name, param_type))

    def formatTable(self):
        """
        Set up the table to look appropriately
        """
        hh = self.table.horizontalHeader()
        hh.setStretchLastSection(True)
        hh.setVisible(False)
        self.table.verticalHeader().setVisible(False)
        self.table.setFocus()
        self.table.hide()
        self.table.resizeColumnsToContents()
        self.table.show()

    def itemChanged(self, item):
        """
        Called every time an item in the table is changed, only value items
        are changeable.
        @param item: The item which has changed
        @type  item: QStandardItem
        """

        name, param_type = item.data()

        # types of setting parameters:
        # 0 ... boolean
        # 1 ... integer
        # 2 ... float
        # 3 ... str
        # 4 ... tuple
        if item.isCheckable():
            checked = item.checkState() == Qt.Checked
            self.settings[name] = checked
        elif param_type in [1, 2]:
            old_value = self.settings[name]
            try:
                new_value = float(item.text())
                if param_type == 1:
                    new_value = int(new_value)
                self.settings[name] = new_value
                item.setText(str(new_value))
            except ValueError:
                logger.error(f"Invalid value '{item.text()}' for parameter '{name}'")
                logger.error(f"Resetting parameter to '{old_value}'.")
                item.setText(str(old_value))
        elif param_type == 3:
            new_value = item.text().strip()
            self.settings[name] = new_value
            item.setText(new_value)
        elif param_type == 4:
            items = [t.strip() for t in item.text().split(",")]
            self.settings[name] = tuple(items)
            item.setText(", ".join(items))
