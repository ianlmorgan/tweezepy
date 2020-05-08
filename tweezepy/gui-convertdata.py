# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 09:38:00 2020

@author: ianmo
"""
import tkinter as tk
from matplotlib.figure import Figure

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.backend_bases import key_press_handler
from pathlib import Path # easy path manipulation
from tkinter import Tk,filedialog


import sys
sys.path.append(r'C:/Users/ianmo/Google Drive/Tweezer Analysis/src/scriptsnew')

#from Plotting import zfplot
from Traceanalysis_v2 import Trace, Refidx, find_surface 

def plotzf(fig,notes):
    df = Trace(notes['dpath'], 
               notes['nrefs'], 
               notes['nexps']).process(notes['refbeads'], 
                                       notes['expbeads'],
                                       idxcorr=notes['idxcorr'], 
                                       pxn=notes['pxn']).df
    ax = fig.subplots(ncols=2,sharex=True)
    #zrkeys = [key for key in df.keys() if 'Zref' in key]
    for i in notes['refbeads']:
        ax[0].plot('Time','Zref%s'%i,data=df,alpha=0.6,c='C%s'%i)
    for i in notes['expbeads']:
        ax[1].plot('Time','Zexp%s'%i,data=df,alpha=0.6,c='C%s'%i)
    #df.plot(ax = ax[0],x='Time',y=zrkeys,alpha = 0.6)
    #zekeys = [key for key in df.keys() if 'Zexp' in key]
    #df.plot(ax = ax[1],x='Time',y=zekeys,alpha = 0.6)
    ax[0].set_ylabel('Z (nm)')
    for a in ax:
        a.set_xlabel('Time (s)')
        a.legend()
    surfs = find_surface(df)
    for surf in surfs:
        ax[1].axhline(surf)
    return ax,surfs    
class Application(tk.Frame):
    def __init__(self, master=None):
        super().__init__(master)
        self.master = master
        self.create_widgets(master)
        self.test()
    def create_widgets(self,master):
        tk.Button(master, text='Select zero force trace file', 
                  command=self.load_data_dir).grid(row=0, column=0,
                                                   columnspan=2)        
        text = ['Salt','Concentration (in mM)', 'Calibration #',
                '# of ref beads', '# of exp beads', 'Pixels per nm',
                'Frame acquisition frequency']
        for i,t in enumerate(text):
            tk.Label(master,text=t).grid(row=i+1, column=0,sticky = tk.W)

        # Create dictionary of entries
        entries = {}
        # Create a Tkinter variable
        tkvar = tk.StringVar(master)
        # Dictionary with options
        choices = ['NaCl','GuHCl','Urea']
        tkvar.set('NaCl')
        entries['salt'] = tkvar
        # Salt dropdown menu
        ddmenu = tk.OptionMenu(master,tkvar,*choices)
        ddmenu.grid(row=1, column=1)
        params = ['conc','cal','nrefs','nexps','pxn','freq']
        # Param entries
        for i,p in enumerate(params):
            entries[p] = tk.Entry(master)
            entries[p].grid(row=i+2,column=1)   
        self.entries = entries
        tk.Button(master, text='Quit', 
                  command=self._quit).grid(row=i+3, column=0, 
                                            sticky=tk.W, pady=4)
        tk.Button(master, text='Show', 
                  command=self.creategraphwindow).grid(row=i+3, column=1, 
                                            sticky=tk.W, pady=1)
                                                       
    def load_data_dir(self):
        initialdir = Path("../")
        dir_name = filedialog.askopenfilename(parent=self.master,
                                              initialdir=initialdir,
                                              title='Select zero force trace file now!')
        if not dir_name:
            print('No file selected.')
        else:
            dpath = Path(dir_name)
            self.dpath = dpath
            
    def creategraphwindow(self):
        entries = self.entries
        notes = {'salt':entries['salt'].get(),
                 'conc':float(entries['conc'].get()),
                 'cal':int(entries['cal'].get()),
                 'nrefs':int(entries['nrefs'].get()),
                 'nexps':int(entries['nexps'].get()),
                 'pxn':float(entries['pxn'].get()),
                 'freq':float(entries['freq'].get())}
        notes['idxcorr'] = Refidx().idxcorr(notes['salt'],notes['conc'])
        notes['refbeads'] = list(range(notes['nrefs']))
        notes['expbeads'] = list(range(notes['nexps']))
        notes['dpath'] = self.dpath
        self.notes = notes
        self.graphwindow = tk.Toplevel(self.master)
        self.app = Graphwindow(self.graphwindow,self.notes)

        # Make test plot
        #fig = Figure(figsize = (5,4),dpi=100)
        #ax = fig.add_subplot(111)
        #x = np.linspace(0,10)
        #y = np.sin(x)
        #ax.plot(x,y)
 
    def _quit(self):
        self.master.quit()     # stops mainloop
        self.master.destroy()  # this is necessary on Windows to prevent
                               # Fatal Python Error: PyEval_RestoreThread: NULL tstate
    def test(self):
        self.entries['salt'].set('NaCl')
        self.entries['conc'].insert(0,'10')
        self.entries['cal'].insert(0,'1')
        self.entries['nrefs'].insert(0,'4')
        self.entries['nexps'].insert(0,'3')
        self.entries['pxn'].insert(0,'118.2')
        self.entries['freq'].insert(0,'400')
        self.dpath = Path(r"C:\Users\ianmo\Google Drive\Tweezer Analysis\data\raw\dsDNA\3-15-20\1trk traces\000")
    
class Graphwindow(tk.Frame):
    def __init__(self, master = None,notes = None):
        self.master = master
        self.notes = notes
        # Reference beads        
        tk.Label(master, text="Reference beads:").grid(row=0,column = 0, sticky=tk.W)
        refs = Checkbar(master, picks = ['Ref %s'%bead for bead in notes['refbeads']], command = self.replot)
        refs.grid(row = 0,column = 1,columnspan=3,sticky = tk.W)
        self.refs = refs
        # Experimental beads
        tk.Label(master, text="Experimental beads:").grid(row=1,column = 0, sticky=tk.W)
        exps = Checkbar(master, picks = ['Exp %s'%bead for bead in notes['expbeads']],command = self.replot)
        exps.grid(row = 1,column = 1,columnspan=3,sticky = tk.W)
        self.exps = exps
        # Button to check and replot refbeads and expbeads
        #tk.Button(master, text='Replot', command=self.replot).grid(row=2,columnspan=10,sticky = tk.W)
        # Setup figure in graphwindow
        self.makefig(master,row = 5)
        # Convert other data
        tk.Button(master, text='Convert other data',command=self.newanalysiswindow).grid(row=8)
    def makefig(self,frame,row):
        fig = Figure(figsize=(12,4),dpi=100)
        self.fig = fig
        ax,surfs = plotzf(fig,self.notes)
        self.ax = ax; self.notes['surfs'] = surfs
        # Create the plot canvas in the graph window
        canvas = FigureCanvasTkAgg(fig, master=frame)  # A tk.DrawingArea.
        canvas.draw() # Draw the figure
        canvas.get_tk_widget().grid(row=row+1,column=0,columnspan=9) # Assign grid spacing
        self.canvas = canvas
        # Toolbar
        # First create toolbar frame
        # Prevents issues between tk.pack and tk.grid
        toolbarFrame = tk.Frame(master=frame)
        # Move/stick toolbar frame to left
        # Prevents it from jumping around
        toolbarFrame.grid(row=row,column=0,columnspan=9,sticky = tk.W+tk.E)
        # Put the toolbar in the toolbar frame
        toolbar = NavigationToolbar2Tk(canvas, toolbarFrame)
        toolbar.update()
        self.toolbar = toolbar
        #canvas.get_tk_widget().grid(row=0,column=0)
        # Is this even working right now?
        # Do I need this?
        canvas.mpl_connect("key_press_event", self.on_key_press)
        
    def on_key_press(self,event):
        print("you pressed {}".format(event.key))
        key_press_handler(event, self.canvas, self.toolbar)
        
    def replot(self):
        self.notes['refbeads'] = [i for i, e in enumerate(list(self.refs.state())) if e == 1]
        self.notes['expbeads'] = [i for i, e in enumerate(list(self.exps.state())) if e == 1]
        #print(self.notes['refbeads'],self.notes['expbeads'])
        for a in self.ax:
            self.fig.clf()
        ax,surfs = plotzf(self.fig,self.notes)
        self.ax = ax;self.surfs = surfs
        self.canvas.draw()
        
    def newanalysiswindow(self):
        self.analysiswindow = tk.Toplevel(self.master)
        self.app = Analysiswindow(self.analysiswindow,self.notes)
    
class Analysiswindow(tk.Frame):
    def __init__(self, master = None,notes = None):
        self.master = master
        self.notes = notes
        
        # Allow salt adjustment
        # Salt label
        tk.Label(master,text="Salt").grid(row=0, column=0,sticky = tk.W)
        # Create a Tkinter variable
        tkvar = tk.StringVar(master)
        # Dictionary with options
        choices = ['NaCl','GuHCl','Urea']
        tkvar.set(notes['salt'])
        self.tkvar = tkvar
        ddmenu = tk.OptionMenu(master,tkvar,*choices)
        ddmenu.grid(row=0, column=1)
        # Concentration label
        tk.Label(master,text="Concentration (in mM)").grid(row=1, column=0,sticky = tk.W)
        # Concentration entry
        conc = tk.Entry(master)
        conc.grid(row=1,column=1,sticky = tk.W)
        conc.insert(0,notes['conc'])
        self.conc = conc
        tk.Button(master, text='Analyze folder',command=self.analyze).grid(row=3,columnspan=2)
        # Work on progressbar for later
        #prog = tk.Progressbar(master,orient='horizontal',length=100.)
    def analyze(self):
        notes = self.notes
        notes['salt'] = self.tkvar.get()
        notes['conc'] = float(self.conc.get())
        initialdir = notes['dpath'].parents[1]
        dir_name = filedialog.askdirectory(parent=self.master,
                                              initialdir=initialdir,
                                              title='Select directory now!')
        if dir_name:
            dir_name2 = dir_name + ' converted'
            path = Path(dir_name)
            path2 = Path(dir_name2)
            if not path2.exists():
                path2.mkdir()
            for p in path.iterdir():
                trace = Trace(p, 
                           notes['nrefs'],
                           notes['nexps'])
                trace.process(notes['refbeads'],
                              notes['expbeads'],
                              idxcorr=notes['idxcorr'], 
                              pxn=notes['pxn'])
                trace.surface_correction(notes['surfs'])
                df = trace.df
                
                fname = p.name + '.csv'
                fpath = path2 / fname
                df.to_csv(fpath,index=False)
class Checkbar(tk.Frame):
   def __init__(self, parent=None, picks=[],command = None, side=tk.LEFT, anchor=tk.W):
      tk.Frame.__init__(self, parent)
      self.vars = []
      for pick in picks:
         var = tk.IntVar(value=1)
         chk = tk.Checkbutton(self, text=pick, command = command, variable=var)
         chk.pack(side=side, anchor=anchor, expand=tk.YES)
         self.vars.append(var)
   def state(self):
      return map((lambda var: var.get()), self.vars)
        
def datadir():
    root = Tk()
    root.withdraw()
    initialdir = Path("C:/Users/ianmo/Google Drive/Tweezer Analysis/data")
    dirname = filedialog.askdirectory(parent=root,initialdir=initialdir,title='Please select a directory')
    dpath = Path(dirname)
    print(dpath)
    return dpath,root


if __name__ == '__main__':
    root = tk.Tk()
    app = Application(master=root)
    app.master.title("My Do-Nothing Application")
    #app.master.maxsize(1000, 400)
    app.mainloop()
    # Select data directory
    #dpath,root = datadir()
    #notes = {'dpath':dpath}