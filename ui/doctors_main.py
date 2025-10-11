#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 18:46:01 2022

@author: Xin liu
"""

import sys

import os

import subprocess

from subprocess import Popen, PIPE

import webbrowser

from PyQt5 import QtCore

from PyQt5.QtWidgets import (

    QApplication, QDialog, QMainWindow, QMessageBox, QFileDialog,

)

from PyQt5.uic import loadUi

#import resource

from mainwindow_ui import Ui_MainWindow

class Window(QMainWindow, Ui_MainWindow):

    def __init__(self, parent=None):

        super().__init__(parent)

        self.setupUi(self)
        
        self.initializeUI()

        self.connectSignalsSlots()

    def initializeUI(self):
        # Update the setting
        self.tabWidget.tabBar().setTabTextColor(0, QtCore.Qt.red)
        self.tabWidget.tabBar().setTabTextColor(1, QtCore.Qt.red)
        self.tabWidget.tabBar().setTabTextColor(2, QtCore.Qt.red)
        self.tabWidget.tabBar().setTabTextColor(3, QtCore.Qt.red)
        self.tabWidget.tabBar().setTabTextColor(4, QtCore.Qt.red)
        self.tabWidget.tabBar().setTabTextColor(5, QtCore.Qt.red)
        self.tabWidget.tabBar().setTabTextColor(6, QtCore.Qt.red)
        
       
        # Add the tooltips
        self.geometryOpenPushButton.setToolTip("Open a dialog box to select CT volume data")



    def connectSignalsSlots(self):

       # self.action_Exit.triggered.connect(self.close)

       self.actionAbout.triggered.connect(self.about)
       
       self.actionUser_Manual.triggered.connect(self.usermanual)
       
       self.geometryOpenPushButton.clicked.connect(self.geometryOpenPushButton_clicked)
       
       self.xsOpenPushButton.clicked.connect(self.xsOpenPushButton_clicked)
       
       self.sourceEnergyPushButton.clicked.connect(self.sourceEnergyPushButton_clicked)
       
       self.detUpdatePushButton.clicked.connect(self.detUpdatePushButton_clicked)
       
       self.inputUpdatePushButton.clicked.connect(self.inputUpdatePushButton_clicked)
       
       self.inputSavePushButton.clicked.connect(self.inputSavePushButton_clicked)
       
       self.launchSolverPushButton.clicked.connect(self.launchSolverPushButton_clicked)
 
       

    def usermanual(self):
        # Open the PDF user manual
        pdfpath = 'DoctorsUserManual_v1.pdf'
        webbrowser.open_new(pdfpath)
        

    def geometryOpenPushButton_clicked(self):

        # Open the CT volume data file
        # 
        
        filename = QFileDialog.getOpenFileName(self, "Select a CT Data File", 
                                                     "./data/",
                                                     "Binary Files (*.bin);;All Files (*)")
        
        print(filename)
        if filename == '':
            print("Error, failed to select CT data file!")
        else:
            self.geometryFileLineEdit.setText(filename[0])
            
        #dialog = FindReplaceDialog(self)

        #dialog.exec()
        
    def xsOpenPushButton_clicked(self):

        # Open cross section data file
        # 
        
        filename = QFileDialog.getOpenFileName(self, "Select a cross section file", 
                                                     "./data/",
                                                     "Text Files (*.dtfr);;All Files (*)")
        
        print(filename)
        if filename == '':
            print("Error, failed to select cross section data file!")
        else:
            self.xsFileLineEdit.setText(filename[0])
            
    def sourceEnergyPushButton_clicked(self):

        # Open source energy spectrum data file
        # 
        
        filename = QFileDialog.getOpenFileName(self, "Select a source spectrum file", 
                                                     "./data/",
                                                     "Text Files (*.spec);;All Files (*)")
        
        print(filename)
        if filename == '':
            print("Error, failed to select source spectrum data file!")
        else:
            self.spectrumLineEdit.setText(filename[0])            
            
    def detUpdatePushButton_clicked(self):

        # Collect varialbes values and Update the input file
        # 
        
        sod = self.sodSpinBox.value()
        odd = self.oddSpinBox.value()
        if(sod != 0.0):
            geomag = (sod + odd) / sod
            self.magValueLabel.setText(str(geomag))
        else:
            msg = QMessageBox(QMessageBox.Critical, "Error", "SOD cannot be zero!")
            msg.exec_()
            
        self.xbinValueLabel.setText(str(self.xBinSpinBox.value()))
        self.ybinValueLabel.setText(str(self.yBinSpinBox.value()))
        self.zbinValueLabel.setText(str(self.zBinSpinBox.value()))
            
        self.detCenterxLabel.setText(str(self.isoXSpinBox.value()))
        self.detCenterzLabel.setText(str(self.isoZSpinBox.value()))
        detCenterY = geomag * (self.isoYSpinBox.value() - self.sourceYDoubleSpinBox.value()) \
                    + self.sourceYDoubleSpinBox.value()
        self.detCenteryLabel.setText(str(detCenterY))
        
        self.detxlenValueLabel.setText(str(self.xLendoubleSpinBox.value() * geomag))
        self.detzlenValueLabel.setText(str(self.zLendoubleSpinBox.value() * geomag))
        
        

    def inputUpdatePushButton_clicked(self):
        #
        # Collect varialbes values and Update the input file
        #
        
        if(self.mcnpCheckBox.isChecked()):
            mcnp_input_gen = 1
        else:
            mcnp_input_gen = 0
            
        xsection_filename = self.xsFileLineEdit.text()
        spec_filename = self.spectrumLineEdit.text()
        ctdata_filename = self.geometryFileLineEdit.text()
        
        source_x = self.sourceXDoubleSpinBox.text()
        source_y = self.sourceYDoubleSpinBox.text()
        source_z = self.sourceZDoubleSpinBox.text()
        source_phi = self.sourcePhiSpinBox.text()
        source_theta = self.sourceThetaSpinBox.text()
        source_number = self.sourceNSpinBox.text()
        source_type = self.sourceTypeComboBox.currentIndex()
        
        iso_x = self.isoXSpinBox.text()
        iso_y = self.isoYSpinBox.text()
        iso_z = self.isoZSpinBox.text()
        
        xbins = self.xBinSpinBox.text()
        ybins = self.yBinSpinBox.text()
        zbins = self.zBinSpinBox.text()
        
        xlen = self.xLendoubleSpinBox.text()
        ylen = self.yLendoubleSpinBox.text()
        zlen = self.zLendoubleSpinBox.text()
        
        material_type = self.matTypeComboBox.currentIndex()
        
        quadrature = self.quadData1ComboBox.currentText()
        
        geomag = self.magValueLabel.text()
        
        # 0: anisotropic; 1: isotropic
        solver_type = self.solverTypeComboBox.currentIndex() - 1
        
        if (self.solverGpuCheckBox.isChecked()):
            GPU_flag = 1
        else:
            GPU_flag = 0
            
        start_ang = self.startAngDoubleSpinBox.value()
        ang_range = self.angularRangeDoubleSpinBox.value()
        num_ang = self.numAngSpinBox.value()
        
        self.inputFileTextEdit.clear()

        self.inputFileTextEdit.append("xsection= " + xsection_filename)
        self.inputFileTextEdit.append("spectrum= " + spec_filename)
        self.inputFileTextEdit.append("source_x= " + source_x)
        self.inputFileTextEdit.append("source_y= " + source_y)
        self.inputFileTextEdit.append("source_z= " + source_z)
        self.inputFileTextEdit.append("source_phi= " + source_phi)
        self.inputFileTextEdit.append("source_theta= " + source_theta)
        self.inputFileTextEdit.append("source_number= " + source_number)
        self.inputFileTextEdit.append("source_type= " + str(source_type))
        self.inputFileTextEdit.append("iso_x= " + iso_x)
        self.inputFileTextEdit.append("iso_y= " + iso_y)
        self.inputFileTextEdit.append("iso_z= " + iso_z)
        self.inputFileTextEdit.append("ctdata= " + ctdata_filename)
        self.inputFileTextEdit.append("xbins= " + xbins)
        self.inputFileTextEdit.append("ybins= " + ybins)
        self.inputFileTextEdit.append("zbins= " + zbins)
        self.inputFileTextEdit.append("xlen= " + xlen)
        self.inputFileTextEdit.append("ylen= " + ylen)
        self.inputFileTextEdit.append("zlen= " + zlen)
        self.inputFileTextEdit.append("geomag= " + geomag)
        self.inputFileTextEdit.append("material_type= " + str(material_type))
        self.inputFileTextEdit.append("quadrature= " + quadrature)
        self.inputFileTextEdit.append("GPU_flag= " + str(GPU_flag))
        self.inputFileTextEdit.append("ISO_flag= " + str(solver_type))
        self.inputFileTextEdit.append("start_angle= " + str(start_ang))
        self.inputFileTextEdit.append("angular_range= " + str(ang_range))
        self.inputFileTextEdit.append("num_angle= " + str(num_ang))
        self.inputFileTextEdit.append("mcnp_input_gen= " + str(mcnp_input_gen))
        
    def inputSavePushButton_clicked(self):
        #
        # Save variales to a file
        #
        path, filter = QFileDialog.getSaveFileName(self, "Save fille", "", "Text files (*.txt)")
        
        if not path:
            msg = QMessageBox(QMessageBox.Warning, "Error!", "Cannot save file!")
            msg.exec_()
            
        self._save_to_path(path)
        
        
        
    def _save_to_path(self, path):
        #
        # save text to file
        #
        text = self.inputFileTextEdit.toPlainText()
        try:
            with open(path, 'w') as f:
                f.write(text)
        except Exception as e:
            msg = QMessageBox(QMessageBox.Warning, "Error!", str(e))
            msg.exec_()
         
    def launchSolverPushButton_clicked(self):
        #
        # Run the LBTE solver
        #
        path, filter = QFileDialog.getOpenFileName(self, "Select fille", "", "Text files (*.txt)")
        
        if not path:
            msg = QMessageBox(QMessageBox.Warning, "Error!", "Cannot open file!")
            msg.exec_()
        
        proc = Popen([r'.\DOCTORS', path], stdout=PIPE, stderr=PIPE, \
                      universal_newlines=True, shell=False)
        
        # Poll proc.stdout to show stdout live
        while True:
            output = proc.stdout.readline()
            if proc.poll() is not None:
                break
            if output:
                print (output)
        rc = proc.poll()
        proc.kill()
        return rc
        
        # stdout, stderr = proc.communicate()
        # print(stderr)
        # print(stdout)
        
        
            
    def about(self):

        QMessageBox.about(

            self,

            "About DOCTORS",

            "<p>Discrete Ordinates Computed TOmography and Radiography Simulator </p>"

            "<p>(i.e. DOCTORS) is a computer code package for producing a neutral </p>"
            
            "<p>particle (e.g. photon) fluence distribution in a voxilized object <p>"
            
            "<p>by solving linear Boltzmann transport equation (LBTE) using discrete </p>"

            "<p>ordinates method and GPU parallel computation techniques. This code </p>"
            
            "<p>package was orginally developed by Edward Norris and Xin Liu at  </p>"
            
            "<p>Missouri University of Science and Technology from 2016~2018. </p>",


        )


class FindReplaceDialog(QDialog):

    def __init__(self, parent=None):

        super().__init__(parent)

        loadUi("ui/find_replace.ui", self)


if __name__ == "__main__":

    app = QApplication(sys.argv)

    win = Window()

    win.show()

    sys.exit(app.exec())