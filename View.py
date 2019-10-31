from ExtraWidgets import *
from sys import platform
from TkinterDnD2 import *
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['axes.linewidth'] = 3
matplotlib.rcParams['xtick.major.width'] = 3
matplotlib.rcParams['ytick.major.width'] = 3
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
plt.ioff()
try:
    import Tkinter as Tk
    import ttk
    import tkFileDialog as FileDialog
    import Queue as queue
except ImportError:
    import tkinter as Tk
    import tkinter.ttk as ttk
    from tkinter import filedialog as FileDialog
    import queue

class Files():
    '''
    Main window for GUI
    Widgets:
        file_list - lists files

    '''
    def __init__(self,que):
        self.root = TkinterDnD.Tk()
        self.root.title('DECA: File List')
        self.root.geometry("600x200+0+0")
        self.menu()
        frame = Tk.Frame(self.root)
        frame.place(relx=0.025, rely=0.025, relheight=0.95, relwidth=0.95)
        self.file_list = ScrolledTreeView(frame)
        self.file_list.pack()
        columns = ('File', 'Proteins', 'States', 'Exposures', 'xCorrected', 'Recombined')
        self.file_list.configure(columns=columns, show='headings',height=8)
        for col in columns:
            self.file_list.column(col, width='100', stretch=0)
            self.file_list.heading(col, text=col)
        self.pb = Progress(frame,que)
        self.popup_menu = Tk.Menu(self.file_list, tearoff=0)
        renameMenu = Tk.Menu(self.file_list, tearoff=0)
        self.popup_menu.add_cascade(label="Rename", menu=renameMenu)
        deleteMenu = Tk.Menu(self.file_list, tearoff=0)
        self.popup_menu.add_cascade(label="Delete", menu=deleteMenu)

    # Window For piping console commands
    def console(self):
        self.console = Tk.Toplevel()
        self.console.title('Console')
        self.console.geometry('400x400+0+0')
        self.console_box = Tk.Text(self.console)
        self.console_box.pack(side=Tk.TOP)
        self.console_copy = ttk.Button(self.console, text='Copy to Clipboard',command= lambda: self.console_copy_action())
        self.console_copy.pack(side=Tk.TOP)

    # Menubar
    def menu(self):
        menubar = Tk.Menu(self.root, font="TkMenuFont")
        self.root.configure(menu=menubar)
        self.file_menu = Tk.Menu(menubar)
        menubar.add_cascade(label="File", menu=self.file_menu)
        self.edit_menu = Tk.Menu(menubar)
        menubar.add_cascade(label="Edit", menu=self.edit_menu)
        self.modify_menu = Tk.Menu(menubar)
        menubar.add_cascade(label="Modify", menu=self.modify_menu)
        self.data_menu = Tk.Menu(menubar)
        menubar.add_cascade(label="Data", menu=self.data_menu)
        self.analyze_menu = Tk.Menu(menubar)
        menubar.add_cascade(label="Analyze", menu=self.analyze_menu)
        self.window_menu = Tk.Menu(menubar)
        menubar.add_cascade(label="Window", menu=self.window_menu)
        self.help_menu = Tk.Menu(menubar)
        menubar.add_cascade(label="Help", menu=self.help_menu)

class Help():
    def __init__(self):
        top = Tk.Toplevel()
        top.geometry("750x500")
        top.title("Documentation")
        top.configure(highlightcolor="black")

        frame = Tk.Frame(top)
        frame.place(relx=0, rely=0, relheight=1, relwidth=1)
        frame.configure(relief=Tk.GROOVE)
        frame.configure(borderwidth="2")
        frame.configure(relief=Tk.GROOVE)
        frame.configure(width=100)

        listbox_files = ScrolledListBox(frame)
        listbox_files.place(relx=0, rely=0, relwidth=0.3, relheight=1)
        listbox_files.configure(selectmode=Tk.SINGLE)
        helplist = {
'0) Welcome':
'''Welcome to the documentation for DECA.\nChoose a topic to the left to get started.
''',
'1) Import New CSV File':
'''MacOS: Select the File window, then click OPEN under FILE in the top menu.
Windows: Click OPEN under FILE in the menu on the File window.
If your file is one of the following, it will be added to the file and data windows:
1) csv state data exported from DynamX
2) csv data exported from HDXWorkbench
3) any csv file to match the organization and headers of the sample file.''',
'2) Back-Exchange Correction':
'''The function addresses two issues of HDX analysis:
1) Sequence dependent loss of deuterium during chromatographic steps
2) Systematic variation between time-points due to automated sample handling\n
Issue 1 may be corrected using a timepoint with fully deuterated controls, where
the back exchange observed in these peptides will be used as correction factors 
for the data from shorter time-points.\n
Issue 2, where higher back exchange is observed at longer time-points when using LEAP
PALs, is readily corrected using correction factors generally for each time-point.\n
This method of correction may also be used to roughly correct data-sets without
fully deuterated controls so long as it has at least one fully exchanging peptide
to use as the basis for generating the correction values.\n
Select a timepoint to correct all peptides on, and/or add exposure correction factors.
''',
'3) Peptide Recombination':
'''This function utilizes peptide overlaps to increase the resolution of the data-set. 
Whenever two peptides share a single terminus, a new peptide will be generated at the
other end, with the uptake calculated from the difference between the two peptides'
uptake values.
To prevent error propagation this should only be performed once per data-set.
''',
'4) Uptake Plot':
'''This function displays the uptake and std data per each peptide.
Multiple states may be displayed on the same plot by SHIFT-clicking them.
To batch export plots, SHIFT-click to select all desired peptides from the list.
''',

'5) Coverage/Heat Map':
'''This function creates a coverage or heat map from the data-set.
The Heat Map is generated using the values from the smallest peptide at each position.
''',

'6) Butterfly Plot':
'''This function creates a butterfly style plot for
the two chosen plots with all exposures shown.
 ''',
'7) Pymol Scripts':
'''This function generates PML files for usage in PyMOL.
These scripts assign Uptake data to b-values for each residue, 
and use these values to assign a coloring scheme.
''',
'8) Retention Time Prediction':
'''Data-sets may be used to predict retention time values for new peptides.
''',
'9) CSV Export':
'''This function will export the data-set with any changes made as the DynamX format.
'''}

        topiclist = sorted(helplist.keys())
        for i in topiclist:
            listbox_files.insert(Tk.END, i)
        textbox = Tk.Text(frame, height=100, width=50, wrap=Tk.WORD)
        textbox.place(relx=0.3, rely=0, relwidth=0.7, relheight=1)
        textbox.configure(state='disabled')
        self.listbox_files = listbox_files
        self.topiclist = topiclist
        self.helplist = helplist
        self.textbox = textbox

class Spectra():
    def __init__(self,main):
        main.widgets.append('Spectra')
        top = Tk.Toplevel(main.root)
        top.title('Spectra Generation')
        combobox_rawfile = ttk.Combobox(top)
        combobox_rawfile.place(relx=0.5, rely=0.2, relwidth=0.9, anchor='center')
        combobox_rawfile['state'] = 'readonly'
        entry_intensity = Tk.Entry(top)
        entry_intensity.place(relx=0.5, rely=0.4, relwidth=0.9, anchor='center')
        combobox_view = ttk.Combobox(top)
        combobox_view.place(relx=0.5, rely=0.6, relwidth=0.9, anchor='center')
        combobox_view['state'] = 'readonly'
        button_save = ttk.Button(top)
        button_save.place(relx=0.5, rely=0.8, relwidth=0.9, anchor='center')
        self.combobox_rawfile = combobox_rawfile
        self.entry_intensity = entry_intensity
        self.combobox_view = combobox_view
        self.button_save = button_save
        self.top = top

class Data():
    def __init__(self, main):
        main.widgets.append('Data')
        self.top = Tk.Toplevel(main.root)
        self.top.title('Data Table')
        main.root.update()
        self.top.geometry("600x450+0+250")
        frame = Tk.Frame(self.top)
        frame.place(relx=0, rely=0, relheight=1, relwidth=1)
        self.tv = ScrolledTreeView(frame)
        self.tv.place(relx=0.025, rely=0.025, relheight=0.95, relwidth=0.95)
        self.tv.configure(selectmode='extended')
        columns = ('Start', 'Stop', 'MaxUptake','Uptake', 'Uptake SD','Uptake_corr','Uptake_SD_corr','RT','Fragment')
        self.tv.configure(columns=columns)
        self.tv.heading('#0', text='Name')
        self.tv.column('#0', width="300", stretch=0)
        self.tv["displaycolumns"] = columns
        for col in columns:
            self.tv.column(col, width="75")
            self.tv.heading(col, text=col)
        #self.tv.delete(*self.tv.get_children())

class Merge():
    def __init__(self, main):
        main.widgets.append('Merge')
        self.main = main
        self.top = Tk.Toplevel(main.root)
        self.top.title("Merge Files")

        frame = Tk.Frame(self.top)
        frame.pack()

        self.file_list = ScrolledListBox(frame)
        self.file_list.grid(row=0,column = 0, columnspan=2)
        self.file_list.configure(height=5,width = 50)

        label = Tk.Label(frame)
        label.grid(row=2,column=0,columnspan=2)
        label.configure(text='''Merge Dataset Name''')

        self.entry_name = Tk.Entry(frame)
        self.entry_name.grid(row=3,column=0,columnspan=2)

        button_add = ttk.Button(frame)
        button_add.grid(row=4,column=0, sticky='E')
        button_add.configure(text='''+''')

        button_remove = ttk.Button(frame)
        button_remove.grid(row=4,column=1, sticky='W')
        button_remove.configure(text='''-''')

        self.button_merge_auto = ttk.Button(frame)
        self.button_merge_auto.grid(row=5,column=0, sticky='E')
        self.button_merge_auto.configure(text='''Auto''')

        self.button_merge_man = ttk.Button(frame)
        self.button_merge_man.grid(row=5,column=1, sticky='W')
        self.button_merge_man.configure(text='''Manual''')

        self.button_help = ttk.Button(frame)
        self.button_help.grid(row=6,column=0, columnspan=2)
        self.button_help.configure(text='''?''')

        self.button_add = button_add
        self.button_remove = button_remove

    def auto(self):
        self.top_auto = Tk.Toplevel(self.top)
        style = ttk.Style()
        style.theme_use('default')
        frame_auto = Tk.Frame(self.top_auto)
        frame_auto.pack(fill=Tk.BOTH, expand=1)
        Tk.Label(frame_auto, text='If matching States, Proteins, or Exposures are found,').pack()
        Tk.Label(frame_auto, text='Merge the following:').pack()
        combobox_keep = ttk.Combobox(frame_auto)
        combobox_keep.pack()
        combobox_keep.configure(values=['None','States','States&Proteins','States,Proteins&Exposures','All, Keep Duplicates', 'All, Ave Duplicates'])
        combobox_keep.set('None')
        combobox_keep['state'] = 'readonly'
        Tk.Label(frame_auto, text='').pack()
        Tk.Label(frame_auto, text='Normalize to Shared State').pack()
        combobox_normalize = ttk.Combobox(frame_auto)
        combobox_normalize.pack()
        combobox_normalize['state'] = 'readonly'
        style.map('TCombobox', fieldbackground=[('readonly', 'white')])
        style.map('TCombobox', foreground=[('readonly', 'black')])
        button_merge_auto_exe = ttk.Button(frame_auto, text = "Merge")
        button_merge_auto_exe.pack()

        self.combobox_keep = combobox_keep
        self.combobox_normalize = combobox_normalize
        self.button_merge_auto_exe = button_merge_auto_exe

    def manual(self):
        self.top_man = Tk.Toplevel(self.top)
        frame_man = Tk.Frame(self.top_man)
        frame_man.pack(fill=Tk.BOTH, expand=1)

        frame1_man = Tk.Frame(frame_man)
        frame1_man.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        tv1 = ttk.Treeview(frame1_man)
        tv1.configure(selectmode='extended')
        tv1.heading('#0', text='Files')
        tv1.pack(side=Tk.LEFT, fill='both', expand='True')
        tv2 = EditableTreeview(frame1_man, columns=('Type', 'File'), show='tree')
        tv2.pack(side=Tk.RIGHT, fill='both', expand='True')
        tv2.configure(selectmode='extended')
        tv2.heading('File', text='File')
        tv2.column('File', width="20")
        tv2["displaycolumns"] = ('File')

        frame2_man = Tk.Frame(frame_man)
        frame2_man.pack(side=Tk.BOTTOM)
        button_add_man = ttk.Button(frame2_man, text='Add')
        button_add_man.pack(side=Tk.LEFT)
        button_delete_man = ttk.Button(frame2_man, text='Delete')
        button_delete_man.pack(side=Tk.LEFT)
        button_addprotein_man = ttk.Button(frame2_man, text='+Protein')
        button_addprotein_man.pack(side=Tk.LEFT)
        button_addstate_man = ttk.Button(frame2_man, text='+State')
        button_addstate_man.pack(side=Tk.LEFT)

        frame3_man = Tk.Frame(frame_man)
        frame3_man.pack(side=Tk.BOTTOM)
        button_merge_man_exe = ttk.Button(frame3_man, text='Merge')
        button_merge_man_exe.pack()

        self.tv1 = tv1
        self.tv2 = tv2
        self.button_add_man = button_add_man
        self.button_delete_man = button_delete_man
        self.button_addprotein_man = button_addprotein_man
        self.button_addstate_man = button_addstate_man
        self.button_merge_man_exe = button_merge_man_exe

class Correct():
    def __init__(self, main):
        main.widgets.append('xCorrect')
        top = Tk.Toplevel(main.root)
        top.title('Back Exchange Correction')
        style = ttk.Style()
        style.theme_use('default')
        frame = Tk.Frame(top)
        frame.place(relx=0, rely=0, relwidth=1, relheight=1)
        top.title("Back-Exchange Correction")
        top.geometry("250x250")
        label1 = Tk.Label(top)
        label1.place(relx=0.2, rely=0.05, relheight = 0.05, anchor='center')
        label1.configure(text="100% D")
        label2 = Tk.Label(top)
        label2.place(relx=0.2, rely=0.15, relheight = 0.05, anchor='center')
        label2.configure(text="Global BX")

        combobox_exposure = ttk.Combobox(frame)
        combobox_exposure.place(relx=0.75, rely=0.05, relwidth=0.4, anchor='center')
        combobox_exposure['state'] = 'readonly'

        entry_global = Tk.Entry(frame)
        entry_global.place(relx=0.75, rely=0.15, relwidth=0.4, anchor='center')

        label3 = Tk.Label(top)
        label3.place(relx=0.5, rely=0.275, relheight = 0.15, anchor='center')
        label3.configure(text="L(ong) E(xposure)\nA(djustment) C(orrection)")
        label4 = Tk.Label(top)
        label4.place(relx=0.5, rely=0.375, relheight = 0.05, anchor='center')
        label4.configure(text="Double-Click to Edit")
        listbox_name = EditableTreeview(frame)
        columns = ('Exposure', 'Correction Factor')
        listbox_name.configure(columns=columns, show='headings')
        for col in columns:
            listbox_name.heading(col, text=col)
            listbox_name.column(col, width="75")

        listbox_name.place(relx=0.1, rely=0.4, relheight=0.45, relwidth=0.8)
        button_close = ttk.Button(frame, text=" Close ", command=lambda: top.destroy)
        button_close.place(relx=0.35, rely=0.925, anchor=Tk.CENTER)
        button_apply = ttk.Button(frame, text=" Apply ")
        button_apply.place(relx=0.65, rely=0.925, anchor=Tk.CENTER)
        button_help = ttk.Button(frame)
        button_help.place(relx=0.9, rely=0.925, height=25, width=25, anchor=Tk.CENTER)
        button_help.configure(takefocus="")
        button_help.configure(text='''?''')

        style.map('TCombobox', fieldbackground=[('readonly', 'white')])
        style.map('TCombobox', foreground=[('readonly', 'black')])
        self.top = top
        self.combobox_exposure = combobox_exposure
        self.entry_global = entry_global
        self.listbox_name = listbox_name
        self.button_apply = button_apply
        self.button_help = button_help

class Replicates():
    def __init__(self, root):
        top = Tk.Toplevel(root)
        top.title('Get Replicates')
        style = ttk.Style()
        style.theme_use('default')
        frame = Tk.Frame(top)
        frame.place(relx=0, rely=0, relwidth=1, relheight=1)
        top.geometry("400x300")
        label = Tk.Label(top)
        label.place(relx=0.5, rely=0.1, relheight = 0.05, anchor='center')
        label.configure(text="Double-Click to Edit")
        listbox_name = EditableTreeview(frame)
        columns = ('State','Exposure','# Replicates')
        listbox_name.configure(columns=columns, show='headings')
        for col in columns:
            listbox_name.heading(col, text=col)
            listbox_name.column(col, width="75")
        listbox_name.place(relx=0.1, rely=0.2, relheight=0.7, relwidth=0.8)
        button_close = ttk.Button(frame, text=" Close ", command=lambda: top.destroy)
        button_close.place(relx=0.35, rely=0.925, anchor=Tk.CENTER)
        button_apply = ttk.Button(frame, text=" Apply ")
        button_apply.place(relx=0.65, rely=0.925, anchor=Tk.CENTER)

        style.map('TCombobox', fieldbackground=[('readonly', 'white')])
        style.map('TCombobox', foreground=[('readonly', 'black')])
        self.top = top
        self.listbox_name = listbox_name
        self.button_apply = button_apply

class Recombine():
    def __init__(self, main):
        main.widgets.append('Recombine')
        top = Tk.Toplevel(main.root)
        top.title('Peptide Overlap Recombination')
        frame = Tk.Frame(top)
        frame.pack()
        Tk.Label(frame, text="   Peptide Recombination").pack()
        Tk.Label(frame, text="").pack()
        Tk.Label(frame, text="  This Function Generates: ").pack()
        text = Tk.Text(frame)
        text.pack()
        button_apply = ttk.Button(frame, text=" Apply Recombination ")
        button_apply.pack()
        button_help = ttk.Button(frame)
        button_help.place(relx=0.95, rely=0.95, height=25, width=25)
        button_help.configure(takefocus="")
        button_help.configure(text='''?''')
        self.top = top
        self.button_apply = button_apply
        self.button_help = button_help
        self.text = text

class Spectrum():
    def __init__(self,root):
        style = ttk.Style()
        if platform == "win32":
            style.theme_use('winnative')
        top = Tk.Toplevel(root)
        top.title("Mass Spectrum")
        top.configure(highlightcolor="black")
        top.resizable(0, 0)
        # top.geometry('450x300+0+0')

        frame1 = Tk.Frame(top)
        frame2 = Tk.Frame(top)
        frame3 = Tk.Frame(top)
        frame1.grid(row=0, column=0, columnspan=5)
        frame2.grid(row=1, column=0, columnspan=4, rowspan=2)
        frame3.grid(row=1, column=5, rowspan=2)

        label_state = Tk.Label(frame1)
        label_state.configure(text='State')
        label_state.pack(side=Tk.LEFT)

        combo_state = ttk.Combobox(frame1)
        combo_state.configure(width=10)
        combo_state.pack(side=Tk.LEFT)

        label_raw = Tk.Label(frame1)
        label_raw.configure(text='Raw')
        label_raw.pack(side=Tk.LEFT)

        combo_raw = ttk.Combobox(frame1)
        combo_raw.configure(width=10)
        combo_raw.pack(side=Tk.LEFT)

        label_exp = Tk.Label(frame1)
        label_exp.configure(text='Exposure')
        label_exp.pack(side=Tk.LEFT)

        combo_exp = ttk.Combobox(frame1)
        combo_exp.configure(width=10)
        combo_exp.pack(side=Tk.LEFT)

        label_z = Tk.Label(frame1)
        label_z.configure(text='Charge')
        label_z.pack(side=Tk.LEFT)

        combo_z = ttk.Combobox(frame1)
        combo_z.configure(width=10)
        combo_z.pack(side=Tk.LEFT)

        canvas = Tk.Canvas(frame2)
        canvas.pack()

        v = Tk.StringVar()
        radio1 = Tk.Radiobutton(frame3)
        radio1.configure(text='All', variable=v, value='all')
        radio1.pack(side=Tk.TOP, anchor=Tk.W)

        radio2 = Tk.Radiobutton(frame3)
        radio2.configure(text='m/z', variable=v, value='mz')
        radio2.pack(side=Tk.TOP, anchor=Tk.W)

        radio3 = Tk.Radiobutton(frame3)
        radio3.configure(text='RT', variable=v, value='rt')
        radio3.pack(side=Tk.TOP, anchor=Tk.W)

        radio4 = Tk.Radiobutton(frame3)
        radio4.configure(text='Mob', variable=v, value='mob')
        radio4.pack(side=Tk.TOP, anchor=Tk.W)

        figure = Figure(figsize=(8, 2), dpi=72, tight_layout=True, facecolor="white")
        plot = FigureCanvasTkAgg(figure, canvas)
        plot.draw()
        plot._tkcanvas.pack(ipadx=20, ipady=20)
        plot._tkcanvas.update()
        ax = figure.add_subplot(111, aspect='auto')

        self.popupMenu = Tk.Menu(plot.get_tk_widget(), tearoff=0)
        if platform == 'darwin':
            plot.get_tk_widget().bind("<Button-2>", lambda event: self.popupMenu.post(event.x_root, event.y_root))
        else:
            plot.get_tk_widget().bind("<Button-3>", lambda event: self.popupMenu.post(event.x_root, event.y_root))

        self.figure = figure
        self.plot = plot
        self.ax = ax
        self.top = top
        self.combo_state = combo_state
        self.combo_raw = combo_raw
        self.combo_exp = combo_exp
        self.combo_z = combo_z

class Significance():
    def __init__(self,root):
        top = Tk.Toplevel(root)
        top.title("Statistical Significance")
        top.configure(highlightcolor="black")
        top.geometry("%dx%d+%d+%d" % (650,685, 1260, 0))
        style = ttk.Style()
        if platform == "win32":
            style.theme_use('winnative')
        else:
            style.theme_use('default')
        frame0 = Tk.Frame(top)
        frame0.pack()
        label0 = Tk.Label(frame0, text='Exposure:')
        label0.grid(row=0, column=0)
        combo_exp = ttk.Combobox(frame0)
        combo_exp.grid(row=0, column=1)

        label1 = Tk.Label(frame0,text='Levene:')
        label1.grid(row=1, column=0)
        levene = Tk.Label(frame0,text='     ')
        levene.grid(row=1, column=1)

        label2 = Tk.Label(frame0,text='ANOVA:')
        label2.grid(row=2, column=0)
        anova = Tk.Label(frame0,text='     ')
        anova.grid(row=2, column=1)

        sep1 = ttk.Separator(frame0)
        sep1.grid(row=3, column=0, columnspan=2)

        frame1 = Tk.Frame(top)
        frame1.pack(fill=Tk.BOTH, expand=1)

        sep1 = ttk.Separator(frame1)
        sep1.pack()
        label3 = Tk.Label(frame1)
        label3.configure(text='Tukey Table')
        label3.pack()
        tukeytable = ScrolledTreeView(frame1)
        tukeytable.pack()
        columns = ('group1', 'group2', 'meandiff', 'lower', 'upper', 'RejectNull', 'p-value')
        tukeytable.configure(columns=columns, show='headings', height=6)
        for col in columns:
            tukeytable.column(col, width='100', stretch=0)
            tukeytable.heading(col, text=col)
        sep2 = ttk.Separator(frame1)
        sep2.pack()
        frame2 = Tk.Canvas(top)
        frame2.pack(expand=1)

        self.figure = Figure(figsize=(8, 6), dpi=72, tight_layout=True, facecolor="white")
        self.ax = self.figure.add_subplot(111, aspect='auto')
        self.plot = FigureCanvasTkAgg(self.figure, master=frame2)
        self.plot.draw()
        self.plot.get_tk_widget().pack(ipadx=20)

        self.popupMenu = Tk.Menu(self.plot.get_tk_widget(), tearoff=0)
        if platform == 'darwin':
            self.plot.get_tk_widget().bind("<Button-2>", lambda event: self.popupMenu.post(event.x_root, event.y_root))
        else:
            self.plot.get_tk_widget().bind("<Button-3>", lambda event: self.popupMenu.post(event.x_root, event.y_root))

        self.top = top
        self.combo_exp = combo_exp
        self.levene = levene
        self.anova = anova
        self.tukeytable = tukeytable
        self.frame2 = frame2

class Stats():
    def __init__(self, root):
        top = Tk.Toplevel(root)
        top.geometry("800x300")
        top.title('Outlier Statistics')
        frame = Tk.Frame(top)
        frame.pack(fill=Tk.BOTH, expand=1)
        tv = ScrolledTreeView(frame, columns=('PPM Mean', 'PPM Std','>10ppm','RT Mean','RT Std','RT>5%','Mob Mean','Mob Std','Mob >5%'))
        tv.configure(selectmode='extended')
        tv.heading('#0', text='Name')
        tv.column('#0', width="150", stretch=0)
        tv.place(relx=0.01,rely=0.01,relheight=.98,relwidth=.98)
        columns = ('PPM Mean', 'PPM Std','>10ppm','RT Mean','RT Std','RT>5%','Mob Mean','Mob Std','Mob >5%')
        tv["displaycolumns"] = columns
        for i in columns:
            tv.column(i, width="75", stretch=0)
            tv.heading(i, text=i)
        self.top = top
        self.tv = tv

class Plot():
    def __init__(self, main):
        self.main = main
        self.top = Tk.Toplevel(main.root)
        self.top.geometry("%dx%d+%d+%d" % (650,685, 605, 0))
        self.top.title('Uptake Plot')
        frame = Tk.Frame(self.top)
        frame.place(relx=0,rely=0,relheight=1,relwidth=1)

        canvas = Tk.Canvas(frame)
        canvas.place(relx=0, rely=0, relheight=1,relwidth=1)
        canvas.configure(highlightthickness=0)

        self.figure = Figure(figsize=(8, 8), dpi=72, tight_layout=True, facecolor="white")
        self.figcan = FigureCanvasTkAgg(self.figure, canvas)
        self.figcan.get_tk_widget().pack(ipadx=20, ipady=20)
        self.canvas = canvas

    class context_menu():
        def __init__(self, frame, main, *args, **kwargs):
            self.main = main
            sizes = range(10,50)
            fonts = sorted(set([f.name for f in matplotlib.font_manager.fontManager.ttflist]))

            popupMenu = Tk.Menu(frame, tearoff=0)

            genMenu = Tk.Menu(popupMenu, tearoff=0)
            fontMenu = Tk.Menu(genMenu, tearoff=0)
            lsMenu = Tk.Menu(genMenu, tearoff=0)
            msMenu = Tk.Menu(genMenu, tearoff=0)
            asMenu = Tk.Menu(genMenu, tearoff = 0)
            arMenu = Tk.Menu(genMenu, tearoff=0)
            flMenu = Tk.Menu(genMenu, tearoff=0)
            genMenu.add_cascade(label="Font", menu=fontMenu)
            genMenu.add_cascade(label="Line Stroke", menu=lsMenu)
            genMenu.add_cascade(label="Marker Size", menu=msMenu)
            genMenu.add_cascade(label="Axis Stroke", menu=asMenu)
            genMenu.add_cascade(label="Aspect Height/Width", menu=arMenu)
            genMenu.add_cascade(label="Line Fit", menu=flMenu)
            for i in list(range(0,11,1)):
                if i == 5:
                    lsMenu.add_command(label=i, command=lambda i=i: self.update('ls', i), background = 'green')
                else:
                    lsMenu.add_command(label=i, command=lambda i=i: self.update('ls', i))
            for i in list(range(0, 11, 1)):
                if i == 3:
                    asMenu.add_command(label=i, command=lambda i=i: self.update('as', i), background = 'green')
                else:
                    asMenu.add_command(label=i, command=lambda i=i: self.update('as', i))
            for i in list(range(0, 21, 2)):
                if i == 10:
                    msMenu.add_command(label=i, command=lambda i=i: self.update('ms', i), background = 'green')
                else:
                    msMenu.add_command(label=i, command=lambda i=i: self.update('ms', i))
            for i in fonts:
                if i == 'Arial':
                    fontMenu.add_command(label=i, command=lambda i=i: self.update('font', i), background = 'green')
                else:
                    fontMenu.add_command(label = i,command=lambda i=i:self.update('font', i))
            for i in [0.7,0.8,0.9,1.0,1.1]:
                if i == 1.0:
                    arMenu.add_command(label = i,command=lambda i=i:self.update('ar', i), background = 'green')
                else:
                    arMenu.add_command(label=i, command=lambda i=i: self.update('ar', i))
            flMenu.add_command(label='Log', command=lambda: self.update('fit', 'log'))
            flMenu.add_command(label='Lin', command=lambda: self.update('fit', 'lin'))
            popupMenu.add_cascade(label="General", menu=genMenu)

            titleMenu = Tk.Menu(popupMenu, tearoff=0)
            fMenut = Tk.Menu(titleMenu, tearoff=0)
            sMenut = Tk.Menu(titleMenu, tearoff=0)

            titleMenu.add_checkbutton(label="Show", variable=main.title_val, onvalue=1, offvalue=0, command=lambda: self.update('title', main.title_val.get()))
            #titleMenu.add_cascade(label="Format", menu=fMenut)
            titleMenu.add_cascade(label="Size", menu=sMenut)
            for i in list(range(16, 65, 4)):
                sMenut.add_command(label=i, command=lambda i=i: self.update('title_size', i))
            popupMenu.add_cascade(label="Title", menu=titleMenu)

            xMenu = Tk.Menu(popupMenu, tearoff=0)
            fMenux = Tk.Menu(xMenu, tearoff=0)
            sMenux = Tk.Menu(xMenu, tearoff=0)
            lMenux = Tk.Menu(xMenu, tearoff=0)
            xMenu.add_checkbutton(label="Show", variable=main.xaxis_val, onvalue=1, offvalue=0, command=lambda: self.update('xaxis', main.xaxis_val.get()))
            #xMenu.add_cascade(label="Format", menu=fMenux)
            xMenu.add_cascade(label="Size", menu=sMenux)
            for i in list(range(16, 65, 4)):
                sMenux.add_command(label=i, command=lambda i=i: self.update('xaxis_size', i))
            #xMenu.add_cascade(label="Labels", menu=lMenux)
            popupMenu.add_cascade(label="X Axis", menu=xMenu)

            yMenu = Tk.Menu(popupMenu, tearoff=0)
            uMenuy = Tk.Menu(yMenu, tearoff=0)
            sMenuy = Tk.Menu(yMenu, tearoff=0)
            lMenuy = Tk.Menu(yMenu, tearoff=0)

            yMenu.add_checkbutton(label="Show", variable=main.yaxis_val, onvalue=1, offvalue=0, command=lambda: self.update('yaxis', main.yaxis_val.get()))
            yMenu.add_cascade(label="Units", menu=uMenuy)
            yMenu.add_cascade(label="Size", menu=sMenuy)
            #yMenu.add_cascade(label="Labels", menu=lMenuy)
            for i in list(range(16, 65, 4)):
                sMenuy.add_command(label=i, command=lambda i=i: self.update('yaxis_size', i))
            uMenuy.add_command(label='Da',command=lambda: self.update('yaxis_units','Da'))
            uMenuy.add_command(label='%Da', command=lambda: self.update('yaxis_units','percentDa'))
            popupMenu.add_cascade(label="Y Axis", menu=yMenu)

            annMenu = Tk.Menu(popupMenu, tearoff=0)
            fMenuann = Tk.Menu(annMenu, tearoff=0)
            sMenuann = Tk.Menu(annMenu, tearoff=0)

            annMenu.add_checkbutton(label="Show", variable=main.ann_val, onvalue=1, offvalue=0, command=lambda: self.update('ann', main.ann_val.get()))
            #annMenu.add_cascade(label="Format", menu=fMenuann)
            annMenu.add_cascade(label="Size", menu=sMenuann)
            for i in list(range(12, 36, 2)):
                sMenuann.add_command(label=i, command=lambda i=i: self.update('ann_size', i))
            popupMenu.add_cascade(label="Annotation", menu=annMenu)

            legMenu = Tk.Menu(popupMenu, tearoff=0)
            sMenuleg = Tk.Menu(legMenu, tearoff=0)

            legMenu.add_checkbutton(label="Show", variable = main.leg_val, onvalue=1, offvalue=0, command=lambda: self.update('leg', main.leg_val.get()))
            legMenu.add_cascade(label="Size", menu=sMenuleg)
            for i in list(range(8, 24, 2)):
                sMenuleg.add_command(label=i, command=lambda i=i: self.update('leg_size', i))
            popupMenu.add_cascade(label="Legend", menu=legMenu)

            statMenu = Tk.Menu(popupMenu, tearoff=0)
            popupMenu.add_cascade(label="State Settings", menu=statMenu)
            self.popupMenu = popupMenu
            if platform == 'darwin':
                frame.bind("<Button-2>", lambda event: popupMenu.post(event.x_root, event.y_root))
            else:
                frame.bind("<Button-3>", lambda event: popupMenu.post(event.x_root, event.y_root))
            self.statMenu = statMenu
            self.main=main

        def update(self,k,v):
            print(str(k)+' '+str(v))
            if k == 'colors' or k== 'markers':
                self.main.settings[k][v[0]] = v[1]
            else:
                self.main.settings[k] = v
            self.main.select_item_data()

    class options(object):
        def __init__(self, parent):
            style = ttk.Style()
            if platform == "win32":
                style.theme_use('winnative')
            top = Tk.Toplevel(parent)
            top.title("Optional Features")
            top.configure(highlightcolor="black")
            top.resizable(0, 0)
            top.geometry('450x300+0+0')

            frame = Tk.Frame(top, padx=5)
            # frame.place(relx=0.0, rely=0.0, relheight=1, relwidth=1)
            frame.grid(row=0, column=0, sticky=Tk.W)

            frame2 = Tk.Frame(top, padx=5)
            # frame.place(relx=0.0, rely=0.0, relheight=1, relwidth=1)
            frame2.grid(row=0, column=1, sticky=Tk.E)

            frame3 = Tk.Frame(top)
            # frame.place(relx=0.0, rely=0.0, relheight=1, relwidth=1)
            frame3.grid(row=1, column=0, columnspan=2, sticky=Tk.S)

            check_title = Tk.Checkbutton(frame)
            # check_title.place(relx=0.0, rely=0.05, relheight=0.05, relwidth=0.6)
            check_title.pack()
            check_title.configure(justify=Tk.LEFT)
            check_title.configure(text='''Title''')
            var_check_title = Tk.IntVar()
            check_title.configure(variable=var_check_title)

            check_legend = Tk.Checkbutton(frame)
            # check_legend.place(relx=0.0, rely=0.05, relheight=0.05, relwidth=0.6)
            check_legend.pack()
            check_legend.configure(justify=Tk.LEFT)
            check_legend.configure(text='''Legend''')
            var_check_legend = Tk.IntVar()
            check_legend.configure(variable=var_check_legend)

            check_info = Tk.Checkbutton(frame)
            # check_info.place(relx=0.0, rely=0.15, relheight=0.05, relwidth=0.6)
            check_info.pack()
            check_info.configure(justify=Tk.LEFT)
            check_info.configure(text='''Peptide Info''')
            var_check_info = Tk.IntVar()
            check_info.configure(variable=var_check_info)

            label_font = Tk.Label(frame)
            # label_font.place(relx=0.05, rely=0.3, relheight=0.05)
            label_font.pack()
            label_font.configure(text='''Font:''')

            combobox_font = ttk.Combobox(frame)
            # combobox_font.place(relx=0.5, rely=0.3, relheight=0.05, relwidth=0.5)
            combobox_font.pack()
            list_fonts = sorted(set([f.name for f in matplotlib.font_manager.fontManager.ttflist]))
            combobox_font.configure(values=list_fonts)
            var_combobox_font = Tk.StringVar()
            var_combobox_font.set('Arial')
            combobox_font.configure(textvariable=var_combobox_font)



            label_xaxis = Tk.Label(frame)
            # label_xaxis.place(relx=0.05, rely=0.5, height=18)
            label_xaxis.pack()
            label_xaxis.configure(text='''X Axis:''')

            entry_xaxis = Tk.Entry(frame)
            # entry_xaxis.place(relx=0.5, rely=0.5, height=20, relwidth=0.5)
            entry_xaxis.pack()
            entry_xaxis.insert(0, 'Time (min)')

            label_yaxis = Tk.Label(frame)
            # label_yaxis.place(relx=0.05, rely=0.55, height=18)
            label_yaxis.pack()
            label_yaxis.configure(text='''Y Axis:''')

            entry_yaxis = Tk.Entry(frame)
            # entry_yaxis.place(relx=0.5, rely=0.55, height=20, relwidth=0.5)
            entry_yaxis.pack()
            entry_yaxis.insert(0, 'Deuterium Uptake (Da)')

            label_xlabels = Tk.Label(frame)
            label_xlabels.pack()
            label_xlabels.configure(text='''X Labels:''')

            entry_xlabels = Tk.Entry(frame)
            # entry_xaxis.place(relx=0.5, rely=0.5, height=20, relwidth=0.5)
            entry_xlabels.pack()
            entry_xlabels.insert(0, 'Fitted')

            label_colors = Tk.Label(frame2)
            # label_colors.place(relx=0.5, rely=0.65, height=18, anchor='center')
            label_colors.pack()
            label_colors.configure(text='''Colors''')

            allcmapkeys = list(mcolors.get_named_colors_mapping().keys())
            newcmapkeys = [i for i in allcmapkeys if (i[0:4] != 'xkcd') and (i[0:3] != 'tab') and (len(i) > 1)]

            tree_colors = ComboboxTreeview(frame2, height = 5, combovals=newcmapkeys)
            # tree_colors.place(relx=0.05, rely=0.69, relheight=0.23, relwidth=0.92)
            tree_colors.pack(fill='x')
            columns1 = ('State', 'Color')
            tree_colors.configure(columns=columns1, show='headings')
            tree_colors.heading('State', text='State')
            tree_colors.column('State', width="75")
            tree_colors.heading('Color', text='Color')
            tree_colors.column('Color', width="125")

            tree_markers = ComboboxTreeview(frame2, height = 5,
                                            combovals=['circle', 'triangle', 'square', 'plus', 'star', 'diamond', 'x'])
            # tree_markers.place(relx=0.05, rely=0.69, relheight=0.23, relwidth=0.92)
            tree_markers.pack(fill='x')
            columns2 = ('State', 'Marker')
            tree_markers.configure(columns=columns2, show='headings')
            tree_markers.heading('State', text='State')
            tree_markers.column('State', width="75")
            tree_markers.heading('Marker', text='Marker')
            tree_markers.column('Marker', width="125")

            button_cancel = ttk.Button(frame3)
            # button_cancel.place(relx=0.11, rely=0.92, height=26, width=68)
            button_cancel.pack(side=Tk.LEFT)
            button_cancel.configure(text='''Cancel''')

            button_save = ttk.Button(frame3)
            # button_save.place(relx=0.54, rely=0.92, height=26, width=56)
            button_save.pack(side=Tk.RIGHT)
            button_save.configure(text='''Save''')

            self.top = top
            self.check_legend = check_legend
            self.var_check_legend = var_check_legend
            self.check_info = check_info
            self.var_check_info = var_check_info
            self.combobox_font = combobox_font
            self.var_combobox_font = var_combobox_font
            self.check_title = check_title
            self.var_check_title = var_check_title
            self.entry_xaxis = entry_xaxis
            self.entry_yaxis = entry_yaxis
            self.entry_xlabels = entry_xlabels
            self.tree_colors = tree_colors
            self.tree_markers = tree_markers
            self.button_cancel = button_cancel
            self.button_save = button_save

    class save(object):
        def __init__(self, parent, multi=False):
            style = ttk.Style()
            if platform == "win32":
                style.theme_use('winnative')
            top = Tk.Toplevel(parent)
            top.title("Save Plot Options")
            top.configure(highlightcolor="black")

            frame = Tk.Frame(top)
            frame.pack()
            frame.configure(relief=Tk.GROOVE)
            frame.configure(borderwidth="2")
            frame.configure(width=185)

            label_dpi = Tk.Label(frame)
            label_dpi.grid(row=0, column=0)
            label_dpi.configure(text='''DPI:''')

            combobox_dpi = ttk.Combobox(frame)
            combobox_dpi.grid(row=0, column=1)
            list_dpi = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600]
            combobox_dpi.configure(values=list_dpi)
            combobox_dpi.set(300)
            combobox_dpi.configure(takefocus="")
            rownum = 1

            label_file_type = Tk.Label(frame)
            label_file_type.grid(row=rownum, column=0)
            label_file_type.configure(text='''Type:''')

            file_type_str = Tk.StringVar()
            combobox_file_type = ttk.Combobox(frame)
            combobox_file_type.grid(row=rownum, column=1)
            list_file_type = ["EPS","PDF","PNG","PS","SVG","TIF"]
            combobox_file_type.configure(values=list_file_type,textvariable=file_type_str)
            combobox_file_type.set("TIF")
            combobox_file_type.configure(takefocus="")
            combobox_file_type.configure(cursor="xterm")
            rownum = rownum + 1

            button_cancel = ttk.Button(frame)
            button_cancel.grid(row=rownum, column=0)
            button_cancel.configure(text='''Cancel''')

            button_save = ttk.Button(frame)
            button_save.grid(row=rownum, column=1)
            button_save.configure(text='''Save''')

            self.top = top
            self.combobox_dpi = combobox_dpi
            self.combobox_file_type = combobox_file_type
            self.file_type_str = file_type_str
            self.button_cancel = button_cancel
            self.button_save = button_save

class Map():
    def __init__(self, main):
        main.widgets.append('Map')
        top = Tk.Toplevel(main.root)
        style = ttk.Style()
        style.theme_use('default')
        top.title("Coverage/Heat Map Generator")
        frame1 = Tk.Frame(top)
        frame1.pack(side=Tk.TOP)
        frame2 = Tk.Frame(top)
        frame2.pack(side=Tk.TOP)
        type_label = Tk.Label(frame1, text='''Type''')
        type_label.grid(row=0,column=0)

        type_combobox = ttk.Combobox(frame1)
        type_combobox.grid(row=0,column=1)
        type_combobox['state'] = 'readonly'

        protein_label = Tk.Label(frame1)
        protein_label.grid(row=1,column=0)
        protein_label.configure(text='''Protein''')

        protein_combobox = ttk.Combobox(frame1)
        protein_combobox.grid(row=1,column=1)
        protein_combobox['state'] = 'readonly'

        state1_label = Tk.Label(frame1)
        state1_label.grid(row=2,column=0)
        state1_label.configure(text='''State 1''')

        state1_combobox = ttk.Combobox(frame1)
        state1_combobox.grid(row=2,column=1)
        state1_combobox['state'] = 'readonly'

        state2_label = Tk.Label(frame1)
        state2_label.grid(row=3,column=0)
        state2_label.configure(text='''State 2''')

        state2_combobox = ttk.Combobox(frame1)
        state2_combobox.grid(row=3,column=1)
        state2_combobox['state'] = 'disabled'

        exposure_label = Tk.Label(frame1,text='''Exposure''')
        exposure_label.grid(row=4,column=0)

        exposure_combobox = ttk.Combobox(frame1)
        exposure_combobox.grid(row=4,column=1)
        exposure_combobox['state'] = 'readonly'

        separator = ttk.Separator(frame1)
        separator.grid(row=5,column=0,columnspan=2)

        color_label = Tk.Label(frame1,text='''Color Scheme:''')
        color_label.grid(row=6,column=0)

        color_combobox = ttk.Combobox(frame1)
        color_combobox.grid(row=6,column=1)
        color_combobox['state'] = 'readonly'

        reverse_label = Tk.Label(frame1, text='''Reverse colors:''')
        reverse_label.grid(row=7, column=0)

        reverse_button = ttk.Checkbutton(frame1)
        reverse_button.grid(row=7, column=1)

        range_label = Tk.Label(frame1,text='''Color Range:''')
        range_label.grid(row=8,column=0)

        range_combobox = ttk.Combobox(frame1)
        range_combobox.grid(row=8,column=1)
        range_combobox['state'] = 'readonly'

        xlim_label = Tk.Label(frame1,text='''Row Length:''')
        xlim_label.grid(row=9,column=0)

        xlim_entry = Tk.Entry(frame1)
        xlim_entry.grid(row=9,column=1)

        cov_button = ttk.Button(frame2,text='''Coverage Map''')
        cov_button.pack(side=Tk.LEFT)

        heat_button = ttk.Button(frame2,text='''Heat Map''')
        heat_button.pack(side=Tk.LEFT)

        style.map('TCombobox', fieldbackground=[('readonly', 'white')])
        style.map('TCombobox', foreground=[('readonly', 'black')])

        style.map('TCombobox', fieldbackground=[('disabled', 'grey')])
        style.map('TCombobox', foreground=[('disabled', 'grey')])
        style.map('TCombobox', foreground=[('disabled', 'grey')])

        self.top = top
        self.frame1 = frame1
        self.frame2 = frame2
        self.state2_label = state2_label
        self.type_combobox = type_combobox
        self.protein_combobox = protein_combobox
        self.state1_combobox = state1_combobox
        self.state2_combobox = state2_combobox
        self.exposure_combobox = exposure_combobox
        self.color_combobox = color_combobox
        self.range_combobox = range_combobox
        self.xlim_entry = xlim_entry
        self.cov_button = cov_button
        self.heat_button = heat_button
        self.reverse_button = reverse_button

class Butterfly():
    def __init__(self, main):
        main.widgets.append('Butterfly')
        top = Tk.Toplevel(main.root)
        style = ttk.Style()
        style.theme_use('default')
        top.title("Butterfly Plot Generator")
        frame1 = Tk.Frame(top)
        frame1.pack()
        frame2 = Tk.Frame(top)
        frame2.pack(side=Tk.TOP)

        label_protein = Tk.Label(frame1)
        label_protein.grid(row=0,column=0)
        label_protein.configure(text='''Protein''')

        combobox_protein = ttk.Combobox(frame1)
        combobox_protein.grid(row=0,column=1)
        combobox_protein['state'] = 'readonly'

        label_state1 = Tk.Label(frame1)
        label_state1.grid(row=1,column=0)
        label_state1.configure(text='''State 1''')

        combobox_state1 = ttk.Combobox(frame1)
        combobox_state1.grid(row=1,column=1)
        combobox_state1['state'] = 'readonly'

        label_state2 = Tk.Label(frame1)
        label_state2.grid(row=2,column=0)
        label_state2.configure(text='''State 2''')

        combobox_state2 = ttk.Combobox(frame1)
        combobox_state2.grid(row=2,column=1)
        combobox_state2['state'] = 'readonly'

        label_exposure = Tk.Label(frame1)
        label_exposure.grid(row=3,column=0)
        label_exposure.configure(text='''Exposure''')

        combobox_exposure = ttk.Combobox(frame1)
        combobox_exposure.grid(row=3,column=1)
        combobox_exposure['state'] = 'readonly'

        button_wings = ttk.Button(frame2)
        button_wings.pack(side=Tk.LEFT)
        button_wings.configure(text='''Wings''')

        button_diff = ttk.Button(frame2)
        button_diff.pack(side=Tk.LEFT)
        button_diff.configure(text='''Diff''')

        button_help = ttk.Button(frame2)
        button_help.pack(side=Tk.LEFT)
        button_help.configure(takefocus="",text='''?''',width=10)

        style.map('TCombobox', fieldbackground=[('readonly', 'white')])
        style.map('TCombobox', foreground=[('readonly', 'black')])

        style.map('TCombobox', fieldbackground=[('disabled', 'grey')])
        style.map('TCombobox', foreground=[('disabled', 'grey')])
        style.map('TCombobox', foreground=[('disabled', 'grey')])

        self.top = top
        self.label_state2 = label_state2
        self.combobox_protein = combobox_protein
        self.combobox_state1 = combobox_state1
        self.combobox_state2 = combobox_state2
        self.combobox_exposure = combobox_exposure
        self.button_wings = button_wings
        self.button_diff = button_diff
        self.button_help = button_help

class PML():
    def __init__(self, main):
        main.widgets.append('PML')
        top = Tk.Toplevel(main.root)
        style = ttk.Style()
        style.theme_use('default')
        top.geometry("443x367")
        top.minsize(450,350)
        top.title("PyMOL Script Export")
        frame = Tk.Frame(top)
        frame.place(relx=0.0, rely=0.0, relheight=1.0, relwidth=1.0)

        label_representation = Tk.Label(frame)
        label_representation.place(relx=0.5, rely=0.025, anchor='center')
        label_representation.configure(text='''Data Representation''')

        combobox_representation = ttk.Combobox(frame)
        combobox_representation.place(relx=0.5, rely=0.075, anchor='center')
        combobox_representation['state'] = 'readonly'

        label_protein = Tk.Label(frame)
        label_protein.place(relx=0.5, rely=0.15, anchor='center')
        label_protein.configure(text='''Protein''')

        combobox_protein = ttk.Combobox(frame)
        combobox_protein.place(relx=0.5, rely=0.2, anchor='center')
        combobox_protein['state'] = 'readonly'

        label_state1 = Tk.Label(frame)
        label_state1.place(relx=0.25, rely=0.275, anchor='center')
        label_state1.configure(text='''State 1''')

        combobox_state1 = ttk.Combobox(frame)
        combobox_state1.place(relx=0.25, rely=0.325, anchor='center')
        combobox_state1['state'] = 'readonly'

        label_state2 = Tk.Label(frame)
        label_state2.place(relx=0.75, rely=0.275, anchor='center')
        label_state2.configure(text='''State 2''')

        combobox_state2 = ttk.Combobox(frame)
        combobox_state2.place(relx=0.75, rely=0.325, anchor='center')
        combobox_state2['state'] = 'disabled'

        label_exposure = Tk.Label(frame)
        label_exposure.place(relx=0.5, rely=0.4, anchor='center')
        label_exposure.configure(text='''Exposure''')

        combobox_exposure = ttk.Combobox(frame)
        combobox_exposure.place(relx=0.5, rely=0.45, anchor='center')
        combobox_exposure['state'] = 'readonly'

        separator = ttk.Separator(frame)
        separator.place(relx=0.5, rely=0.5, relwidth=1, anchor='center')

        button_pdb = ttk.Button(frame)
        #button_pdb.place(relx=0.5, rely=0.55, anchor='center')
        button_pdb.configure(text='''Select PDB (optional)''')

        label_pdbprotein = Tk.Label(frame)
        label_pdbprotein.place(relx=0.25, rely=0.625, anchor='center')
        label_pdbprotein.configure(text='''PDB Model:''')

        entry_protein = Tk.Entry(frame)
        entry_protein.place(relx=0.65, rely=0.625, anchor='center')

        label_pdbchain = Tk.Label(frame)
        label_pdbchain.place(relx=0.25, rely=0.7, anchor='center')
        label_pdbchain.configure(text='''PDB Chain:''')

        entry_chain = Tk.Entry(frame)
        entry_chain.place(relx=0.65, rely=0.7, anchor='center')

        label_color = Tk.Label(frame)
        label_color.place(relx=0.25, rely=0.775, anchor='center')
        label_color.configure(text='''Color Scheme:''')

        combobox_color = ttk.Combobox(frame)
        combobox_color.place(relx=0.65, rely=0.775, anchor='center')

        label_range = Tk.Label(frame)
        label_range.place(relx=0.25, rely=0.85, anchor='center')
        label_range.configure(text='''Color Range:''')

        combobox_range = ttk.Combobox(frame)
        combobox_range.place(relx=0.65, rely=0.85, anchor='center')
        combobox_range['state'] = 'readonly'

        button_save = ttk.Button(frame)
        button_save.place(relx=0.5, rely=0.95, anchor='center')
        button_save.configure(text='''Save''')

        button_help = ttk.Button(frame)
        button_help.place(relx=0.9, rely=0.95, anchor='center')
        button_help.configure(takefocus="")
        button_help.configure(text='''?''')

        style.map('TCombobox', fieldbackground=[('readonly', 'white')])
        style.map('TCombobox', foreground=[('readonly', 'black')])

        style.map('TCombobox', fieldbackground=[('disabled', 'grey')])
        style.map('TCombobox', foreground=[('disabled', 'grey')])
        style.map('TCombobox', foreground=[('disabled', 'grey')])

        self.top = top
        self.label_state2 = label_state2
        self.combobox_representation = combobox_representation
        self.combobox_protein = combobox_protein
        self.combobox_state1 = combobox_state1
        self.combobox_state2 = combobox_state2
        self.combobox_exposure = combobox_exposure
        self.button_pdb = button_pdb
        self.button_save = button_save
        self.entry_protein = entry_protein
        self.entry_chain = entry_chain
        self.combobox_color = combobox_color
        self.combobox_range = combobox_range
        self.button_help = button_help

class RT():
    def __init__(self, root):
        top = Tk.Toplevel(root)
        label_sequence = Tk.Label(top)
        label_sequence.configure(text='Peptide Sequence:')
        label_sequence.pack()
        text = Tk.Entry(top)
        text.pack()
        label_RT = Tk.Label(top)
        label_RT.configure(text='RT:')
        label_RT.pack()
        button_calculate = ttk.Button(top)
        button_calculate.configure(text='Calculate')
        button_calculate.pack()
        self.top = top
        self.text = text
        self.label_RT = label_RT
        self.button_calculate = button_calculate

class Label():
    def __init__(self,root):
        self.root = root
        self.top = Tk.Toplevel(self.root)
        self.stv = Tk.StringVar()
        self.label = Tk.Label(self.top,textvariable = self.stv)
        self.stv.set('NA')
        self.label.pack()

class Progress():
    def __init__(self,top,que):
        self.que = que
        self.que.put(0)
        self.top = top
        self.frame = Tk.Frame(self.top)
        self.loading_label = Tk.Label(self.frame,text='Loading')
        self.loading_label.pack(side=Tk.LEFT,anchor=Tk.CENTER,expand=1,fill=Tk.BOTH)
        self.dvar = Tk.DoubleVar()
        self.progbar = ttk.Progressbar(self.frame,variable=self.dvar,maximum=100)
        self.progbar.pack(side=Tk.LEFT,anchor=Tk.CENTER,expand=1,fill=Tk.BOTH)
        self.num_label = Tk.Label(self.frame, textvariable=self.dvar)
        self.num_label.pack(side=Tk.LEFT,anchor=Tk.CENTER,expand=1,fill=Tk.BOTH)