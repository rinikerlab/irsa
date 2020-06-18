from PyQt5 import QtCore, QtGui, QtWidgets
import res_rc
import ui as gui
import global_ as g
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

import pickle
import sys
import random
import numpy as np
import argparse
import Reader
import Spectrum
import Algorithm
class Button:
    def __init__(self,button,entry):
        self.button = button
        self.entry = entry
        g.values[self.entry] = self.button.value()
        self.button.valueChanged.connect(self.valuechange)
    def valuechange(self):
        g.values[self.entry] = self.button.value()
        if(self.entry == "lb" or self.entry == "hb"):
            for i in range(len(g.canvas_list)):
                try:
                    g.canvas_list[i].change_x()
                except:
                    pass
        elif(self.entry == "w" or self.entry == "T"):
            try:
                g.Spectrum.Boltzmann_weight_IR()
                g.canvas_list[2].plot_IR_theo()
                g.Spectrum.Boltzmann_weight_VCD()
                g.canvas_list[3].plot_VCD_theo()
            except:
                pass
class Click_Button:
    def __init__(self,button,entry,args = None):
        self.button = button
        self.entry = entry
        self.button.clicked.connect(self.click)
        self.args = args
    def click(self):
        if(self.entry == "normalize_1"):
            try:
                tmp = (g.exp_ir[:,0] <= g.values["hb"]) & (g.exp_ir[:,0] >= g.values["lb"])
                g.exp_ir[:,1] = g.exp_ir[:,1]/np.max(g.exp_ir[tmp,1])
                g.canvas_list[0].plot_IR()
            except:
                pass
            try:
                tmp = (g.exp_vcd[:,0] <= g.values["hb"]) & (g.exp_vcd[:,0] >= g.values["lb"])
                g.exp_vcd[:,1] = g.exp_vcd[:,1]/np.max(np.abs(g.exp_vcd[tmp,1]))
                g.canvas_list[1].plot_VCD()            
            except:
                pass
        elif(self.entry == "normalize_2"):
            try:
                tmp = (g.theo_ir[:,0] <= g.values["hb"]) & (g.theo_ir[:,0] >= g.values["lb"])
                g.theo_ir[:,1] = g.theo_ir[:,1]/np.max(g.theo_ir[tmp,1])
                g.canvas_list[2].plot_IR()
            except:
                pass
            try:
                tmp = (g.theo_vcd[:,0] <= g.values["hb"]) & (g.theo_vcd[:,0] >= g.values["lb"])
                g.theo_vcd[:,1] = g.theo_vcd[:,1]/np.max(np.abs(g.theo_vcd[tmp,1]))
                g.canvas_list[3].plot_VCD()            
            except:
                pass
        elif(self.entry == "automatic"):
            try:
                tmp = (g.exp_ir[:,0] <= g.values["hb"]) & (g.exp_ir[:,0] >= g.values["lb"])
                tmp_ir = np.asarray(g.exp_ir[tmp])
                g.peak_list_x = []
                g.peak_list_y = []
                g.peak_list_VCD_y = []
                for i in range(1,len(tmp_ir)-1):
                    if(tmp_ir[i-1,1]<=tmp_ir[i,1]>=tmp_ir[i+1,1]):
                        g.peak_list_x.append(tmp_ir[i,0])
                        g.peak_list_y.append(tmp_ir[i,1])
                print(g.peak_list_x)
                g.canvas_list[0].plot_peaks()
                g.exp_peaks = np.zeros((len(g.peak_list_x),2))
                g.exp_peaks[:,0] = np.asarray(g.peak_list_x)
                g.exp_peaks[:,1] = np.asarray(g.peak_list_y)
                for peak in g.peak_list_x:
                    g.peak_list_VCD_y.append(g.exp_vcd[abs(g.exp_vcd[:,0]-peak)<10e-1,1][0])
                g.canvas_list[1].plot_peaks_VCD()
            except:
                pass
        elif(self.entry == "align"):
            try:
                del Algo
                print("del")
            except: 
                pass
            Algo = Algorithm.Algorithm()
            if(g.set_VCD==False):
                g.returnvalue, g.old_freq, g.freq_new, g.inten_new = Algo.Needleman_IR()
            else:
                g.returnvalue, g.old_freq, g.freq_new, g.inten_new,g.inten_VCD_new = Algo.Needleman_IR()
            g.canvas_list[4].plot_IR_assigned()
            g.Spectrum.IR_shifted()
            if(g.set_VCD==True):
                g.Spectrum.VCD_shifted()
                p_ir,p_vcd = g.Spectrum.integrate()

                self.args.setText("Score: " + str(g.returnvalue)[0:6]+"\np_ir: " + str(p_ir)[0:4]+"\np_vcd: " + str(p_vcd)[0:4]+"\n")
            else:
                p_ir = g.Spectrum.integrate()
                self.args.setText("Score: " + str(g.returnvalue)[0:6]+"\np_ir: " + str(p_ir)[0:4]+"\n")
            g.canvas_list[5].plot_IR_shifted()
class Load_Button:
    def __init__(self,button,entry):
        self.button = button
        self.entry = entry
        self.button.clicked.connect(self.click)
    def click(self):
        if(self.entry == "experimental IR"):
            try:
                fileName, _ = QtWidgets.QFileDialog.getOpenFileName(None,"Select "+self.entry +" Spectrum","","")
                g.exp_ir = np.loadtxt(fileName,usecols=(0,1))
                g.exp_ir = g.exp_ir[g.exp_ir[:,0].argsort()]
                g.canvas_list[0].plot_IR()
                g.set_IR = True
            except:
                pass
        elif(self.entry == "experimental VCD"):
            try:
                fileName, _ = QtWidgets.QFileDialog.getOpenFileName(None,"Select "+self.entry +" Spectrum","","")
                g.exp_vcd = np.loadtxt(fileName,usecols=(0,1,))
                g.exp_vcd = g.exp_vcd[g.exp_vcd[:,0].argsort()]
                g.canvas_list[1].plot_VCD()
                g.set_VCD = True
            except:
                pass
        elif(self.entry == "energies"):
            fileName_energy, _ = QtWidgets.QFileDialog.getOpenFileName(None,"Select "+self.entry +" Energies","","")
            g.E = pickle.load(open(fileName_energy,"rb"))
        elif(self.entry == "theoretical IR"):
            fileName_IR, _ = QtWidgets.QFileDialog.getOpenFileName(None,"Select "+self.entry +" Spectrum","","")
            g.theo_ir = pickle.load(open(fileName_IR,"rb"))
            g.Spectrum.Boltzmann_weight_IR()
            g.canvas_list[2].plot_IR_theo()
        elif(self.entry == "theoretical VCD"):
            try:
                fileName_VCD, _ = QtWidgets.QFileDialog.getOpenFileName(None,"Select "+self.entry +" Spectrum","","")
                g.theo_vcd = pickle.load(open(fileName_VCD,"rb"))
                g.Spectrum.Boltzmann_weight_VCD()
                g.canvas_list[3].plot_VCD_theo()
            except:
                pass
class Canvas(FigureCanvas):
    def __init__(self, single = None, parent = None, Button = None, dpi = 100):
        height = parent.height()/100.
        width = parent.width()/100.
        fig = Figure(figsize=(width, height), dpi=dpi)
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)
        self.ax = self.figure.add_subplot(111)
        #Button
    def change_x(self):
        self.ax.set_xlim(g.values["lb"],g.values["hb"])
        self.draw()
    def plot_IR(self):
        self.delete()
        self.ax.plot(g.exp_ir[:,0],g.exp_ir[:,1],color="black")
        self.ax.set_ylim(0,1.05)
        self.ax.set_xlim(g.values["lb"],g.values["hb"])
        self.draw()
    def plot_IR_shifted(self):
        self.delete()
        self.ax.plot(g.exp_ir[:,0],g.exp_ir[:,1],color="black")
        self.ax.plot(g.IR_shifted[:,0],g.IR_shifted[:,1],color="red")
        if(g.set_VCD==True):
            self.ax.plot(g.exp_vcd[:,0],g.exp_vcd[:,1],"--",color="black")
            self.ax.plot(g.VCD_shifted[:,0],g.VCD_shifted[:,1],"--",color="red")
        self.ax.set_ylim(-1.05,1.05)
        self.ax.set_xlim(g.values["lb"],g.values["hb"])
        self.draw()
    def plot_peaks(self):
        self.delete()
        self.ax.plot(g.exp_ir[:,0],g.exp_ir[:,1],color="black")
        self.ax.set_ylim(0,1.05)
        self.ax.set_xlim(g.values["lb"],g.values["hb"])
        self.ax.plot(g.peak_list_x,g.peak_list_y,"o",color="blue")
        self.draw()
    def plot_peaks_VCD(self):
        self.delete()
        self.ax.plot(g.exp_vcd[:,0],g.exp_vcd[:,1],"--",color="black")
        self.ax.plot(g.peak_list_x,g.peak_list_VCD_y,"o",color="blue")
        self.ax.set_ylim(-1.05,1.05)
        self.ax.set_xlim(g.values["lb"],g.values["hb"])
        self.draw()
    def delete(self):
        try:
            while(len(self.ax.lines)>0):
                self.ax.lines[-1].remove()
        except:
            pass
    def plot_VCD(self):
        self.delete()
        self.ax.plot(g.exp_vcd[:,0],g.exp_vcd[:,1],"--",color="black")
        self.ax.set_ylim(-1,1)
        self.ax.set_xlim(g.values["lb"],g.values["hb"])
        self.draw()
    def plot_IR_theo(self):
        self.delete()
        self.ax.plot(g.spectrum_boltzmann[:,0],g.spectrum_boltzmann[:,1],color="red") ##x_axis 0..2000
        self.ax.set_ylim(0,1)
        self.ax.set_xlim(g.values["lb"],g.values["hb"])
        tmp_ir = g.spectrum_boltzmann[g.values["lb"]:g.values["hb"]]
        g.theo_peaks_x = []
        g.theo_peaks_y = []
        for i in range(1,len(tmp_ir)-1):
            if(tmp_ir[i-1,1]<=tmp_ir[i,1]>=tmp_ir[i+1,1]):
                g.theo_peaks_x.append(tmp_ir[i,0])
                g.theo_peaks_y.append(tmp_ir[i,1])
        self.ax.plot(g.theo_peaks_x,g.theo_peaks_y,"o",color="blue")
        g.theo_peaks = np.zeros((len(g.theo_peaks_x),2))
        g.theo_peaks[:,0] = np.asarray(g.theo_peaks_x)
        g.theo_peaks[:,1] = np.asarray(g.theo_peaks_y)
        self.draw()
    def plot_IR_assigned(self):
        self.delete()
        #g.returnvalue, g.old_freq, g.freq, g.inten
        self.ax.plot(g.spectrum_boltzmann[:,0],g.spectrum_boltzmann[:,1],color="red") ##x_axis 0..2000
        self.ax.plot(g.exp_ir[:,0],g.exp_ir[:,1],color="black")
        if(g.set_VCD==True):
            self.ax.plot(g.exp_vcd[:,0],g.exp_vcd[:,1],"--",color="black")
            self.ax.plot(g.spectrum_boltzmann_vcd[:,0],g.spectrum_boltzmann_vcd[:,1],"--",color="red")
        for i in range(len(g.old_freq)):
            self.ax.plot([g.old_freq[i],g.freq_new[i]],[g.inten_new[i],g.inten_new[i]],color="blue")
        self.ax.set_ylim(-1,1)
        self.ax.set_xlim(g.values["lb"],g.values["hb"])
        tmp_ir = g.spectrum_boltzmann[g.values["lb"]:g.values["hb"]]

        self.draw()
    def plot_VCD_theo(self):
        self.delete()
        self.ax.plot(g.spectrum_boltzmann_vcd[:,0],g.spectrum_boltzmann_vcd[:,1],"--",color="red") ##x_axis 0..2000
        tmp_vcd = g.spectrum_boltzmann_vcd[g.values["lb"]:g.values["hb"]]
        self.ax.set_ylim(-1,1)
        g.peak_list_VCD_y_theo = []
        for peak in g.theo_peaks_x:
            g.peak_list_VCD_y_theo.append(tmp_vcd[np.abs(tmp_vcd[:,0]-peak)<10e-3,1][0])
        self.ax.plot(g.theo_peaks_x,g.peak_list_VCD_y_theo,"o",color="blue")
        self.ax.set_xlim(g.values["lb"],g.values["hb"])
        self.draw()
if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = gui.Ui_MainWindow()
    ui.setupUi(MainWindow)
    g.Spectrum = Spectrum.Spectrum()
    g.canvas_list = []
    g.canvas_list.append(Canvas(parent = ui.exp_ir_graph))
    g.canvas_list.append(Canvas(parent = ui.exp_vcd_graph))
    g.canvas_list.append(Canvas(parent = ui.theo_ir_graph))
    g.canvas_list.append(Canvas(parent = ui.theo_vcd_graph))
    g.canvas_list.append(Canvas(parent = ui.assignment_graph))
    g.canvas_list.append(Canvas(parent = ui.shifted_graph))
    g.list_buttons = []
    g.list_buttons.append(Button(ui.w,"w"))
    g.list_buttons.append(Button(ui.lb,"lb"))
    g.list_buttons.append(Button(ui.hb,"hb"))
    g.list_buttons.append(Button(ui.sigma_1,"s0"))
    g.list_buttons.append(Button(ui.sigma_2,"s1"))
    g.list_buttons.append(Button(ui.mu,"mu"))
    g.list_buttons.append(Button(ui.cutoff,"c"))
    g.list_buttons.append(Button(ui.temperature,"T"))

    g.list_buttons.append(Load_Button(ui.Load_EXP,"experimental IR"))
    g.list_buttons.append(Load_Button(ui.Load_EXP_VCD,"experimental VCD"))
    g.list_buttons.append(Load_Button(ui.load_theo_exp,"theoretical IR"))
    g.list_buttons.append(Load_Button(ui.load_theo_vcd,"theoretical VCD"))
    g.list_buttons.append(Load_Button(ui.load_theo_exp_2,"energies"))
    #g.list_buttons.append(Load_Button(ui.load_theo_vcd,"theoretical VCD"))


    g.list_buttons.append(Click_Button(ui.normalize_1,"normalize_1"))
    g.list_buttons.append(Click_Button(ui.normalize_2,"normalize_2"))
    g.list_buttons.append(Click_Button(ui.automatic,"automatic"))
    g.list_buttons.append(Click_Button(ui.align,"align",args = ui.results))

    MainWindow.show()
    sys.exit(app.exec_())
