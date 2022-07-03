# 2020-05-XX Added details to map names
# 2020-05-XX Added blue-white-magenta color scheme
# 2020-05-XX Added Custom ranges for maps

# 2020-06-04 Adjusted maps to starts at first covered residue
# 2020-06-04 Adjusted maps to fix colorbars on different linelengths

print("Loading, please wait")
import View
from Model import *
from operator import itemgetter
from itertools import cycle,islice
from os.path import split, splitext, expanduser, exists
from os import system, makedirs, environ, path
from sys import platform, exit, exc_info
import sys
from copy import copy, deepcopy
import threading
import numpy as np
from scipy.optimize import leastsq
import queue
import textwrap
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import rcParams
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['axes.linewidth'] = 3
matplotlib.rcParams['xtick.major.width'] = 3
matplotlib.rcParams['ytick.major.width'] = 3
from matplotlib.widgets import Slider
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from mpl_toolkits.mplot3d import Axes3D
plt.ioff()
if getattr(sys, 'frozen', False):
    application_path = sys._MEIPASS
else:
    application_path = path.dirname(path.abspath(__file__))
TK_DND_PATH = path.join(application_path,'tkdnd2.8')

from TkinterDnD2 import *
import atexit

try:
    import Tkinter as Tk
    import ttk
    import tkFileDialog as FileDialog
except ImportError:
    import tkinter as Tk
    import tkinter.ttk as ttk
    from tkinter import filedialog as FileDialog

class Main():
    '''
    Main Controller Class
    Windows:
        Files
        Data
        Back Exchange
        Peptide Segmentation
        Plot
        Maps
        Butterfly Plots
        Significance
        Ion Stats
    Start Files window and right click menu
    Set Drag and Drop bindings to Files widget
    Set Button click bindings to update all windows upon file selections

    '''
    # Initialize GUI
    def __init__(self):
        self.widgets = []
        self.widget_dict = {}
        self.state_obj = {}
        self.file_list = []
        self.dir = ''
        # Initialize Main Window
        self.view = View.Files(que)
        self.widgets.append('Main')
        self.view.window_menu.add_command(label='Main', command=lambda: self.focus(self.root))
        self.proteinsvar = Tk.Variable()
        self.statesvar = Tk.Variable()
        self.exposuresvar = Tk.Variable()
        self.peptidesvar = Tk.Variable()
        self.root = self.view.root
        self.menu()
        self.drag_and_drop()
        self.context_menu()

        # # Output console
        # stderr = StdRedirector(self.console_box)
        # stdout = StdRedirector(self.console_box)

        # Focus on new window
        if platform == 'darwin':  # How Mac OS X is identified by Python
            system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "Python" to true' ''')

        # Initialize Data Window
        self.data = Data(self)
        self.view.file_list.bind('<ButtonRelease-1>', self.select_item)

        # Fill Data Window if state_obj already filled
        if self.state_obj != {}:
            for key in self.state_obj.keys():
                fileid = self.view.file_list.insert('', 'end', text=key, values=(
                    key, str(', '.join(self.state_obj[key].proteins)),
                    str(', '.join(self.state_obj[key].states)),
                    str(', '.join([str(i) for i in self.state_obj[key].exposures])),
                    self.state_obj[key].corrected, "No"))
                self.corrected = self.state_obj[key].corrected
                self.view.file_list.selection_set(fileid)
                self.view.file_list.focus(fileid)
                self.select_item()
                print('Initialize ' + str(self.proteinsvar.get()))

    # Setup Drag and Drop to import files
    def drag_and_drop(self):
        self.view.file_list.drop_target_register(DND_FILES)
        self.view.file_list.dnd_bind('<<Drop>>', lambda event: self.open_file(event))

    # Pipes console messages
    def console_copy_action(self):
        #self.console.withdraw()
        self.console.clipboard_clear()
        self.console.clipboard_append(self.console_box.get(1.0, "end-1c"))
        #self.console.update()

    # Creates Right click menu for the file list
    def context_menu(self):
        self.view.popup_menu.add_command(label="Open", command=self.open_file)
        self.view.popup_menu.add_command(label="Close", command=self.close_file)

        # renameMenu.add_command(label="File", command=self.rename)
        # renameMenu.add_command(label="State", command=self.rename)
        # renameMenu.add_command(label="Protein", command=self.rename)
        #
        # deleteMenu.add_command(label="File", command=self.delete)
        # deleteMenu.add_command(label="State", command=self.delete)
        # deleteMenu.add_command(label="Protein", command=self.delete)

        self.view.file_list.bind("<Button-2>", lambda event: self.view.popup_menu.post(event.x_root, event.y_root))

    # Command to Focus on the specified window
    def focus(self, window):
        window.attributes('-topmost', 1)
        window.grab_set()
        window.grab_release()
        window.attributes('-topmost', 0)

    # Creates and fills Menubar
    def menu(self):
        self.view.file_menu.add_command(label="Open File", command=self.open_file)
        self.view.file_menu.add_command(label="Merge", command=self.merge)
        self.view.file_menu.add_command(label="Exit", command=lambda: self.root.destroy_window())

        self.view.edit_menu.add_command(label="Rename File/Protein/State", state="disabled")
        self.view.edit_menu.add_command(label="Delete Data", state="disabled")
        self.view.edit_menu.add_command(label="Find Incomplete Data", state="disabled")
        self.view.edit_menu.add_command(label="Filter Data", state="disabled")

        self.view.modify_menu.add_command(label="Correct Back Exchange", command=lambda: self.open_window('xCorrect'))
        self.view.modify_menu.add_command(label="Get Peptide Overlaps", command=lambda: self.open_window('Recombine'))

        self.view.data_menu.add_command(label="Mass Spectra", command=lambda: Spectra(self))
        self.view.data_menu.add_command(label="Uptake Plot", command=lambda: self.open_window('Plot'))
        self.view.data_menu.add_command(label="Coverage Map", command=lambda: self.open_window('Map'))
        self.view.data_menu.add_command(label="Butterfly Plot", command=lambda: self.open_window('Butterfly'))
        self.view.data_menu.add_command(label="PyMOL Script", command=lambda: self.open_window('PML'))
        self.view.data_menu.add_command(label='Export CSV', command=lambda: self.save_data())
        self.view.data_menu.add_command(label='Export Plot', command=lambda: self.plot.save())
        self.view.data_menu.entryconfig('Export Plot', state='disabled')

        self.view.analyze_menu.add_command(label="Retention Time", command=lambda: RT(self))
        self.view.analyze_menu.add_command(label="Spectrum", command=lambda: self.open_window('Spectrum'))
        self.view.analyze_menu.add_command(label='Significance', command=lambda: self.open_window('Significance'))
        self.view.analyze_menu.add_command(label='Ion Stats', command=lambda: self.open_window('Ion Stats'))

        self.view.help_menu.add_command(label="About", command=lambda:
        self.popup('About Deca', 'DECA: hydrogen Deuterium Exchange Correction and Analysis\nCreated by '
                                 'Ryan Lumpkin\nKomives Lab, UCSD\n2019'))
        self.view.help_menu.add_command(label="Documentation", command=lambda: self.help())
        self.view.help_menu.add_command(label='Log Errors',command=lambda: self.on_exit())

    # Creates or opens windows based on selection from menu
    def open_window(self,key):
        if key in self.widgets:
            self.focus(self.widget_dict[key])
        else:
            if key == 'xCorrect':
                self.correct = Correct(self)
            elif key == 'Recombine':
                self.recombine = Recombine(self)
            elif key == 'Plot':
                self.plot = Plot(self)
            elif key == 'Map':
                self.map = Map(self)
            elif key == 'Butterfly':
                self.map = Butterfly(self)
            elif key == 'PML':
                self.pml = PML(self)
            elif key == 'Spectrum':
                self.spectrum = Spectrum(self)
            elif key == 'Significance':
                self.significance = Significance(self)
            elif key == 'Ion Stats':
                self.stats = Stats(self)
            else:
                return

    # Prompts a file to open, imports file into DECA class, updates progess bar
    def open_file(self,item=None):
        print('Opening File')
        self.increment = 0
        if item == None:
            csv_file = FileDialog.askopenfilename(title="Choose CSV or DNX file",
                                                  filetypes=[('csv files', '*.csv'),('dnx files', '*.DnX'), ('all files', '*.*')], )
        else:
            if item.data[0] == '{':
                csv_file = str(item.data[1:-1])
            else:
                csv_file = str(item.data)
        if csv_file != '':
            print(csv_file)
            if ('Data' not in self.widgets):
                self.data = Data(self)
            csv_file_name = str(splitext(split(csv_file)[1])[0])
            if str(splitext(split(csv_file)[1])[1]).lower() in ['.csv','.dnx']:
                self.view.file_list.unbind('<ButtonRelease-1>')
                self.view.file_list.unbind('<<Drop>>')
                n = 1
                if csv_file_name in self.state_obj.keys():
                    while csv_file_name in self.state_obj.keys():
                        csv_file_name = csv_file_name + '(%s)' % n
                        n = n + 1
                    self.popup('Error!', 'Name renamed to:'+str(csv_file_name))
                # self.view.pb.frame.pack(side=Tk.TOP)
                # def thread_queue():
                try:
                    d = DECA(csv_file,que=que)
                except:
                    print('Exception')
                    self.view.file_list.bind('<ButtonRelease-1>', self.select_item)
                    self.view.file_list.dnd_bind('<<Drop>>', lambda event: self.open_file(event))
                    return
                self.state_obj[csv_file_name] = d
                if self.state_obj[csv_file_name].recombined:
                    recombined = 'Yes'
                else:
                    recombined = 'No'
                id = self.view.file_list.insert('', 'end', text=csv_file_name, values=(
                    csv_file_name, str(', '.join(self.state_obj[csv_file_name].proteins)),
                    str(', '.join(self.state_obj[csv_file_name].states)),
                    str(', '.join([str(i) for i in self.state_obj[csv_file_name].exposures])),
                    self.state_obj[csv_file_name].corrected, recombined))
                self.corrected = self.state_obj[csv_file_name].corrected
                self.view.file_list.selection_set(id)
                self.view.file_list.focus(id)
                self.select_item()
                self.file_list.append(csv_file)
                self.view.file_list.bind('<ButtonRelease-1>', self.select_item)
                self.view.file_list.dnd_bind('<<Drop>>', lambda event: self.open_file(event))
                #
                # def check_finished():
                #     while True:
                #         try:
                #             x = self.view.pb.que.get_nowait()
                #         except queue.Empty:
                #             if secondary_thread.is_alive():
                #                 self.root.after(50,check_finished)
                #                 break
                #             else:
                #                 self.view.pb.progbar.stop()
                #                 self.view.pb.frame.pack_forget()
                #                 print('Secondary Thread Dead')
                #                 break
                #         else:  # continue from the try suite
                #             if x == -1:
                #                 self.view.pb.progbar.stop()
                #                 self.view.pb.frame.pack_forget()
                #                 print('Secondary Thread Finished')
                #                 break
                #             else:
                #                 self.view.pb.dvar.set(x)
                #
                # secondary_thread = threading.Thread(target=thread_queue)
                # secondary_thread.start()
                #
                # # check the Queue in 50ms
                # self.root.after(50,check_finished)

            else:
                return
        else:
            self.popup("Open File Error", "No Item Chosen")

    # Removes a file from the file list
    def close_file(self):
        del self.state_obj[self.csvfile]
        self.view.file_list.delete(self.curitem)
        if 'Data' in self.widgets:
            self.data.delete()

    # Triggers the Merge feature
    def merge(self):
        self.merge = Merge(self)

    # Updates data list and plots upon the selection of an item in the file list
    def select_item(self, *args):
        print('Item Selected')
        if ('Data' not in self.widgets):
            print('Creating DataCtrl')
            self.data = Data(self)
        self.curitem = self.view.file_list.focus()
        self.csvfile = self.view.file_list.item(self.curitem)['text']
        self.dir = split(self.state_obj[self.csvfile].filename)[0]
        if self.csvfile != '':
            self.corrected = self.view.file_list.set(self.curitem, 'xCorrected')
            self.data.delete()
            self.data.insert_data(self.state_obj[self.csvfile])
            self.proteinsvar.set(tuple(sorted(self.state_obj[self.csvfile].proteins)))
            self.statesvar.set(tuple(sorted(self.state_obj[self.csvfile].states)))
            self.exposuresvar.set(tuple([str(i) for i in self.state_obj[self.csvfile].exposures]))
            self.set_colors = {}
            colors = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black']
            color_cycle = cycle(colors)
            for state in sorted(self.state_obj[self.csvfile].states):
                self.set_colors[state] = next(color_cycle)
        if str(splitext(split(self.state_obj[self.csvfile].filename)[1])[1]).lower() == '.dnx':
            self.view.data_menu.entryconfig('Mass Spectra',state='normal')
            self.view.analyze_menu.entryconfig('Spectrum', state='normal')
            self.view.analyze_menu.entryconfig('Ion Stats', state='normal')
        else:
            self.view.data_menu.entryconfig('Mass Spectra', state='disabled')
            self.view.analyze_menu.entryconfig('Spectrum', state='disabled')
            self.view.analyze_menu.entryconfig('Ion Stats', state='disabled')

    # Starts the Tkinter Mainloop
    def run(self):
        self.root.mainloop()

    # Popup Message
    def popup(self,titlestr,messagestr):
        top = Tk.Toplevel()
        top.title(titlestr)
        label = Tk.Label(top)
        label.pack(fill=Tk.BOTH, expand=True)
        label.configure(text=messagestr)

    # Help Documents
    def help(self, topic='0) Welcome'):
        def updatehelp():
            if len(h.listbox_files.curselection()) != 0:
                h.textbox.configure(state='normal')
                h.textbox.delete(1.0, Tk.END)
                h.textbox.insert(Tk.END,h.helplist[h.topiclist[h.listbox_files.curselection()[0]]])
                h.textbox.configure(state='disabled')
        h = View.Help()
        h.textbox.configure(state='normal')
        h.textbox.insert(Tk.END, h.helplist[topic])
        h.textbox.configure(state='disabled')
        h.listbox_files.bind('<<ListboxSelect>>', lambda e: updatehelp())

    # Custom save dialog for MacOS
    def save_figure(self, dpival, fig): # Thanks to Donal Fellows https://stackoverflow.com/questions/52282334/tkfiledialog-does-not-show-file-extension-options-on-osx10-12-6
        root = Tk.Toplevel()
        filetypes = "{{{eps files} *.eps} {{jpg files} *.jpg} {{pdf files} *.pdf} {png files} *.png} {{ps files} *.ps}" \
                    " {{svg files} *.svg} {{tif files} *.tif} {{all files} *}"
        filename = root.tk.eval('::tk::dialog::file:: save -filetypes {' + filetypes + '}')
        if filename == "":
            filename = None
        if split(filename)[0] != "":
            self.dir = split(filename)[0]
        fig.savefig(filename, dpi=dpival, transparent=True)
        root.destroy()

    # Properly close window and remove window from registry
    def on_closing(self,widget, name):
        self.widgets.remove(name)
        self.view.window_menu.delete(name)
        self.widget_dict[name] = None
        widget.destroy()

    # Export data to CSV
    def save_data(self):
        filename = FileDialog.asksaveasfilename(title="Choose save location", defaultextension=".csv",
                                                initialfile=self.csvfile + (self.corrected == 'Yes') * "_Corrected",
                                                initialdir=self.dir)
        if split(filename)[0] != "":
            self.dir = split(filename)[0]
        self.state_obj[self.csvfile].exportData(filename)
        # 2020-05-XX Create peptides for each individual plot
        # self.state_obj[self.csvfile].exportPlotData(filename)

    def on_exit(self):
        self.console = Console()
        sys.stderr = self.console
        sys.stdout = self.console
        atexit.register(self.console.exit)

class Spectra():
    '''
    Spectra Controller
    Widgets:
        Rawfile Combobox
        View Combobox
        Intensity Entry box
    Fill content
    Bind buttons
    '''

    # Initialize Spectra window, Fill widgets
    def __init__(self,main):
        self.main = main
        self.view = View.Spectra(main)
        main.view.window_menu.add_command(label='Spectra', command=lambda: main.focus(self.view.top))
        self.view.top.protocol("WM_DELETE_WINDOW", lambda: main.on_closing(self.view.top, 'Spectra'))
        self.view.button_save.configure(text='Generate', command = lambda: self.gen_spectra())
        rawfiles = main.state_obj[main.csvfile].rawfiles
        self.view.combobox_rawfile.configure(values=[(i['RawID'],i['State'],i['Exposure']) for i in rawfiles])
        self.view.combobox_rawfile.set('None')
        vcmd = (self.view.top.register(self.validate),
                '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')
        self.view.entry_intensity.configure(validate = 'key', validatecommand = vcmd)
        self.view.entry_intensity.insert(10000)
        self.view.combobox_view.configure(values=['3D','mz_rt','mz_mob','rt_mob'])
        self.view.combobox_view.set('None')

    # Ensures intensity values are numeric
    def validate(self, value_if_allowed, text):
        if text in '0123456789.':
            try:
                float(value_if_allowed)
                return True
            except ValueError:
                return False
        else:
            return False

    # Graphs spectra
    def gen_spectra(self):
        plt.ion()
        rawfile = self.view.combobox_rawfile.get()[0]
        intensity = self.view.entry_intensity.get()
        view = self.view.combobox_view.get()
        if rawfile == None:
            commands = ['Int > ' + str(intensity)]
        else:
            commands = ['Int > '+str(intensity), 'RawID = '+str(rawfile)]
        print(commands)
        ions = self.main.state_obj[self.main.csvfile].getDnxIons(commands)
        print(len(ions))
        nullions = self.main.state_obj[self.main.csvfile].getDnxIons(commands, null=True)
        print(len(nullions))
        x = [i['mz'] for i in ions]
        y = [i['RT'] for i in ions]
        z = [i['Drift'] for i in ions]
        xn = [i['mz'] for i in nullions]
        yn = [i['RT'] for i in nullions]
        zn = [i['Drift'] for i in nullions]
        fig = plt.figure()
        if (view == None) or (view == '3D'):
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(xn, yn, zn, c='b', s=1, marker='.')
            ax.scatter(x, y, z, c='r', s=1, marker='.')
            ax.set_xlabel('M/Z')
            ax.set_ylabel('Retention Time')
            ax.set_zlabel('Mobility')
            ax.view_init(elev=30, azim=225)
        elif view == 'mz_rt':
            ax = fig.add_subplot(111)
            ax.scatter(xn, yn, c='b', s=1, marker='.')
            ax.scatter(x, y, c='r', s=1, marker='.')
            ax.set_xlabel('M/Z')
            ax.set_ylabel('Retention Time')
        elif view == 'mz_mob':
            ax = fig.add_subplot(111)
            ax.scatter(xn, zn, c='b', s=1, marker='.')
            ax.scatter(x, z, c='r', s=1, marker='.')
            ax.set_xlabel('M/Z')
            ax.set_ylabel('Mobility')
        elif view == 'rt_mob':
            ax = fig.add_subplot(111)
            ax.scatter(yn, zn, c='b', s=1, marker='.')
            ax.scatter(y, z, c='r', s=1, marker='.')
            ax.set_xlabel('Retention Time')
            ax.set_ylabel('Mobility')
        plt.draw()
        plt.ioff()

class StdRedirector():
    '''
    Pipes Console messages to GUI widget
    From: #https://www.reddit.com/r/learnprogramming/comments/3vq0dm/python_how_can_i_print_text_out_in_the_gui_rather/
    '''

    # This class pipes console messages to a text widget
    def __init__(self, text_widget):
        self.text_space = text_widget

    # Writes a command to the text widget
    def write(self, string):
        self.text_space.config(state=Tk.NORMAL)
        self.text_space.insert("end", string)
        self.text_space.see("end")
        self.text_space.config(state=Tk.DISABLED)

class Data():
    '''
    Data Controller
    Widgets:
        Peptide tree
    Binds Mouse click and arrow key events
    Sets Right Click Menu
    '''
    # Initialize Data window, Create event bindings, start right click menu
    def __init__(self, main):
        self.main = main
        self.view = View.Data(main)
        main.widget_dict['Data'] = self.view.top
        main.view.window_menu.add_command(label='Data', command=lambda: main.focus(self.view.top))
        self.view.top.protocol("WM_DELETE_WINDOW", lambda: main.on_closing(self.view.top, 'Data'))
        self.view.tv.bind('<ButtonRelease-1>', lambda e: self.update())
        self.view.tv.bind_all('<Up>', lambda e: self.update())
        self.view.tv.bind_all('<Down>', lambda e: self.update())
        self.context()

    # Selection update triggers updates to other windows if open
    def update(self):
        if 'Significance' in self.main.widgets:
            self.main.significance.update_data()
        if 'Spectrum' in self.main.widgets:
            self.main.spectrum.update()
        if 'Plot' in self.main.widgets:
            self.main.plot.select_item_data()
        if 'Stats' in self.main.widgets:
            self.main.stats.getData()

    # Fill Right click menu, Create event binding
    def context(self):
        self.popup_menu = Tk.Menu(self.view.tv, tearoff=0)
        self.popup_menu.add_command(label="Sort")
        self.popup_menu.add_command(label="Filter")
        self.popup_menu.add_command(label="Delete")
        self.popup_menu.add_command(label="Show")
        self.popup_menu.add_command(label="Hide")
        self.view.tv.bind("<Button-2>", lambda event: self.popup_menu.post(event.x_root, event.y_root))

    # Adds content to window from item selected in main file list
    # Creates Xref to keep track of item position and the corresponding peptide
    def insert_data(self, obj):
        print('Inserting Data')
        self.xref = {}
        for protein in obj.proteins:
            print('Protein')
            protein_tvid = self.view.tv.insert('', 'end', text=protein, open=False)
            all_peptides = [i for i in obj.master_csv if i['Protein'] == protein]
            pepids = [i[0] for i in
                      sorted(
                          sorted(
                              list(
                                  set([(i['PepID'],i['Start'],i['End']) for i in all_peptides])),
                              key=itemgetter(2)),
                          key=itemgetter(1))]

            for pepid in pepids:
                print(pepid)
                all_states = [i for i in all_peptides if i['PepID'] == pepid]
                states = set([i['State'] for i in all_states])
                start = int(all_states[0]['Start'])
                stop = int(all_states[0]['End'])
                sequence = str(all_states[0]['Sequence'])
                maxuptake = str(round(float(all_states[0]['MaxUptake']), 2))
                modification = str(all_states[0]['Modification'])
                fragment = str(all_states[0]['Fragment'])
                pep_tvid = self.view.tv.insert(protein_tvid, 'end', text=sequence, open=False,
                                               values = (start, stop, maxuptake,'','','','','',fragment))
                self.xref[pep_tvid] = all_states
                for state in states:
                    all_exposures = [i for i in all_states if i['State'] == state]
                    exposures = [i['Exposure'] for i in all_exposures]
                    state_tvid = self.view.tv.insert(pep_tvid, 'end', text=state, open=False,
                                               values = (start, stop, maxuptake,'','','','','',fragment))
                    for exposure in exposures:
                        for row in [i for i in all_exposures if i['Exposure']==exposure]:
                            if row['Uptake'] is None:
                                uptake = 'None'
                                uptakesd = 'None'
                                rt = 'None'
                            else:
                                uptake = str(round(float(row['Uptake']), 2))
                                uptakesd = str(round(float(row['Uptake SD']), 2))
                                rt = str(row['RT'])
                            frag = str(row['Fragment'])
                            if obj.corrected == 'Yes':
                                try:
                                    uptake_corr = str(round(float(row['Uptake_corr']), 2))
                                except Exception as ex:
                                    print(row)
                                uptake_sd_corr = str(round(float(row['Uptake_SD_corr']), 2))
                            else:
                                uptake_corr = ''
                                uptake_sd_corr = ''
                            exp_tvid = self.view.tv.insert(state_tvid, 'end', text=exposure, open=False,
                                                           values=(start, stop, maxuptake, uptake, uptakesd, uptake_corr,
                                                                                  uptake_sd_corr, rt, frag))

    # Removes content from window
    def delete(self):
        self.view.tv.delete(*self.view.tv.get_children())

class Merge():
    '''
    Merge Window Controller
    Initial Window Prompts files to merge
    Subsequent windows facilitate manual a=or automatic merging
    '''
    # Initialize merge window, Sets button bindings
    def __init__(self, main):
        self.main = main
        self.files = []
        self.view = View.Merge(main)
        main.widget_dict['Merge'] = self.view.top
        self.main.view.window_menu.add_command(label='Merge', command=lambda: self.main.focus(self.view.top))
        self.view.button_merge_auto.configure(command = lambda: self.auto())
        self.view.button_merge_man.configure(command = lambda: self.manual())
        self.view.button_help.configure(command=lambda e: self.main.help('1) Import'))
        self.view.button_add.configure(command=lambda: self.addfile())
        self.view.button_remove.configure(command=lambda:self.delfile())
        self.view.top.protocol("WM_DELETE_WINDOW", lambda: main.on_closing(self.view.top, 'Merge'))

    # Adds an item to the list of files to be merged
    def addfile(self):
        file = FileDialog.askopenfilename(title="Choose file", filetypes=(('csv files', '*.csv'),('all files', '*.*')))
        if file == '':
            return
        elif file in self.files:
            self.main.popup('Error!', 'File already added')
            return
        else:
            self.view.file_list.insert(Tk.END, file)
            self.files.append(file)

    # Removes an item from the list of files to be merged
    def delfile(self):
        del self.files[self.view.file_list.curselection()[0]]
        self.view.file_list.delete(self.view.file_list.curselection()[0])

    # Triggers the manual merge feature
    def manual(self):
        def add_state():
            id = self.view.tv2.insert('', 'end', text='NewState', values=('state'))
        def add_protein():
            if len(self.view.tv2.selection()) != 0 and self.view.tv2.item(self.view.tv2.selection())['values'][
                0] == 'state':
                self.view.tv2.item(self.view.tv2.selection(), open='True')
                id = self.view.tv2.insert(self.view.tv2.selection(), 'end', text='NewProtein', values=('protein'))
            else:
                return
        def add_selection():
            if len(self.view.tv1.selection()) == 0:
                return
            for i in self.view.tv1.selection():
                if self.view.tv1.item(i)['values'][0] == 'file' and len(self.view.tv2.selection()) == 0:
                    for s in self.view.tv1.get_children(item=i):
                        sid = self.view.tv2.insert('', 'end', text=self.view.tv1.item(s)['text'],
                                                  values=self.view.tv1.item(s)['values'])
                        for p in self.view.tv1.get_children(item=s):
                            pid = self.view.tv2.insert(sid, 'end', text=self.view.tv1.item(p)['text'],
                                                      values=self.view.tv1.item(p)['values'])
                            for e in self.view.tv1.get_children(item=p):
                                self.view.tv2.insert(pid, 'end', text=self.view.tv1.item(e)['text'],
                                                    values=self.view.tv1.item(e)['values'])
                    continue
                elif self.view.tv1.item(i)['values'][0] == 'state' and len(self.view.tv2.selection()) == 0:
                    sid = self.view.tv2.insert('', 'end', text=self.view.tv1.item(i)['text'],
                                              values=self.view.tv1.item(i)['values'])
                    for p in self.view.tv1.get_children(item=i):
                        pid = self.view.tv2.insert(sid, 'end', text=self.view.tv1.item(p)['text'],
                                                  values=self.view.tv1.item(p)['values'])
                        for e in self.view.tv1.get_children(item=p):
                            self.view.tv2.insert(pid, 'end', text=self.view.tv1.item(e)['text'],
                                                values=self.view.tv1.item(e)['values'])
                    continue
                elif self.view.tv1.item(i)['values'][0] == 'state' and len(self.view.tv2.selection()) != 0:
                    if self.view.tv2.item(self.view.tv2.selection())['values'][0] == 'state':
                        for p in self.view.tv1.get_children(item=i):
                            pid = self.view.tv2.insert(self.view.tv2.selection(), 'end',
                                                      text=self.view.tv1.item(p)['text'],
                                                      values=self.view.tv1.item(p)['values'])
                            for e in self.view.tv1.get_children(item=p):
                                self.view.tv2.insert(pid, 'end', text=self.view.tv1.item(e)['text'],
                                                    values=self.view.tv1.item(e)['values'])
                        continue
                    else:
                        continue
                elif self.view.tv1.item(i)['values'][0] == 'protein' and self.view.tv2.selection() != 0:
                    if self.view.tv2.item(self.view.tv2.selection())['values'][0] == 'state':
                        pid = self.view.tv2.insert(self.view.tv2.selection(), 'end', text=self.view.tv1.item(i)['text'],
                                                  values=self.view.tv1.item(i)['values'])
                        for e in self.view.tv1.get_children(item=i):
                            self.view.tv2.insert(pid, 'end', text=self.view.tv1.item(e)['text'],
                                                values=self.view.tv1.item(e)['values'])
                        continue
                    elif self.view.tv2.item(self.view.tv2.selection())['values'][0] == 'protein':
                        for e in self.view.tv1.get_children(item=i):
                            self.view.tv2.insert(self.view.tv2.selection(), 'end', text=self.view.tv1.item(e)['text'],
                                                values=self.view.tv1.item(e)['values'])
                        continue
                    else:
                        continue
                elif self.view.tv1.item(i)['values'][0] == 'exposure' and self.view.tv2.selection() != 0:
                    if self.view.tv2.item(self.view.tv2.selection())['values'][0] == 'protein':
                        self.view.tv2.insert(self.view.tv2.selection(), 'end', text=self.view.tv1.item(i)['text'],
                                            values=self.view.tv1.item(i)['values'])
                        continue
                    else:
                        continue
                else:
                    continue
        def del_selection():
            for s in self.view.tv2.selection():
                if self.view.tv2.exists(s):
                    self.view.tv2.delete(s)
        def merge():
            dictlist = {}
            states = self.view.tv2.get_children()
            for state in states:
                state_name = self.view.tv2.item(state)['text']
                dictlist[state_name] = {}
                proteins = self.view.tv2.get_children(item=state)
                for protein in proteins:
                    protein_name = self.view.tv2.item(protein)['text']
                    dictlist[state_name][protein_name] = {}
                    exposures = self.view.tv2.get_children(item=protein)
                    for exposure in exposures:
                        exposure_name = self.view.tv2.item(exposure)['text']
                        if exposure_name not in dictlist[state_name][protein_name].keys():
                            dictlist[state_name][protein_name][exposure_name] = []
                        origin = self.view.tv2.item(exposure)['values']
                        dictlist[state_name][protein_name][exposure_name].append(origin[1])
                        for row in merge_each[origin[1]].data_nest[origin[2]][origin[3]][exposure_name]:
                            row['Protein'] = protein_name
                            row['State'] = state_name
                            row['File'] = origin[1]
            for file in merge_each.keys():
                merge_each[file].organizeData()
            self.main.state_obj[self.name].mergeFiles(dictlist, merge_each)
            id = self.main.view.file_list.insert('', 'end', text=self.name, values=(
                self.name, str(', '.join(self.main.state_obj[self.name].proteins)),
                str(', '.join(self.main.state_obj[self.name].states)),
                str(', '.join([str(i) for i in self.main.state_obj[self.name].exposures])),
                self.main.state_obj[self.name].corrected, "No"))
            self.view.top.destroy()
            self.main.view.file_list.selection_set(id)
            self.main.view.file_list.focus(id)
            self.main.select_item()
        files = self.view.file_list.get(0, Tk.END)
        self.name = str(self.view.entry_name.get())
        merge_each = {}
        n = 1
        if self.name in self.main.state_obj.keys():
            while self.name in self.main.state_obj.keys():
                self.name = self.name + '(%s)' % n
                n = n + 1
            self.main.popup('Error!', 'Name renamed to:' + str(self.name))
        self.main.state_obj[self.name] = DECA(empty=True)
        for i in range(0, len(files)):
            filename = str(split(files[i])[1])
            merge_each[filename] = DECA(files[i])
        self.view.manual()
        for file in merge_each.keys():
            fileid = self.view.tv1.insert('', 'end', text=file, values=('file'))
            states = {state:{} for state in merge_each[file].states}
            for protein in sorted(merge_each[file].data_nest.keys()):
                for state in merge_each[file].data_nest[protein].keys():
                    states[state][protein] = merge_each[file].data_nest[protein][state].keys()
            for state in states.keys():
                stateid = self.view.tv1.insert(fileid, 'end', text=state, values=('state'))
                for protein in sorted(states[state].keys()):
                    proteinid = self.view.tv1.insert(stateid, 'end', text=protein, values=('protein'))
                    for exposure in states[state][protein]:
                        self.view.tv1.insert(proteinid, 'end', text=exposure, values=('exposure', file, protein, state))
        self.view.button_add_man.configure(command=add_selection)
        self.view.button_delete_man.configure(command=del_selection)
        self.view.button_addprotein_man.configure(command=add_protein)
        self.view.button_addstate_man.configure(command=add_state)
        self.view.button_merge_man_exe.configure(command=merge)

    # Triggers the automatic merge feature
    def auto(self):
        def merge():
            keep_val = self.view.combobox_keep.get()
            if keep_val == 'None':
                s_bool = False
                p_bool = False
                e_bool = False
                d_bool = False
            elif keep_val == 'States':
                s_bool = True
                p_bool = False
                e_bool = False
                d_bool = False
            elif keep_val == 'States&Proteins':
                s_bool = True
                p_bool = True
                e_bool = False
                d_bool = False
            elif keep_val == 'States,Proteins&Exposures':
                s_bool = True
                p_bool = True
                e_bool = True
                d_bool = False
            elif keep_val == 'All, Keep Duplicates':
                s_bool = True
                p_bool = True
                e_bool = True
                d_bool = True
            elif keep_val == 'All, Ave Duplicates':
                s_bool = True
                p_bool = True
                e_bool = True
                d_bool = True
            else:
                print('Error')
            norm = self.view.combobox_normalize.get()
            states_sum = set()
            proteins_sum = set()
            exposures_sum = set()

            for i in range(0, len(files)):
                filename = str(split(files[i])[1])
                merge_each[filename] = DECA(files[i])
                states_sum.update(merge_each[filename].states)
                proteins_sum.update(merge_each[filename].proteins)
                exposures_sum.update(merge_each[filename].exposures)

            p_dict = dict.fromkeys(proteins_sum)
            ps_dict = dict.fromkeys(proteins_sum)

            for p in p_dict.keys():
                p_dict[p] = []
                ps_dict[p] = {}
                for i in merge_each.keys():
                    if p in merge_each[i].data_nest.keys():
                        p_dict[p].append(i)

            s_dict = dict.fromkeys(states_sum)
            sp_dict = dict.fromkeys(states_sum)
            spe_dict = dict.fromkeys(states_sum)

            for s in s_dict.keys():
                s_dict[s] = set()
                sp_dict[s] = {}
                spe_dict[s] = {}
                proteins = list()
                for i in merge_each:
                    for p in merge_each[i].proteins:
                        if s in merge_each[i].data_nest[p].keys():
                            s_dict[s].add(i)
                            proteins.append((p,i))
                for protein in list(set([p[0] for p in proteins])):
                    sp_dict[s][protein] = []
                    spe_dict[s][protein] = {}
                for p in proteins:
                    sp_dict[s][p[0]].append(p[1])
                for p in spe_dict[s].keys():
                    exposures = list()
                    for i in sp_dict[s][p]:
                        for e in merge_each[i].data_nest[p][s].keys():
                            exposures.append((e,i))
                    for exposure in list(set([e[0] for e in exposures])):
                        spe_dict[s][p][exposure] = []
                    for e in exposures:
                        spe_dict[s][p][e[0]].append(e[1])

            new_s_dict = {}
            new_sp_dict = {}
            new_spe_dict = {}

            for s in s_dict.keys():
                if len(s_dict[s]) > 1:
                    if not s_bool: # If not merging matches, then rename
                        n = 0
                        for i in s_dict[s]:
                            state_name = s + '-' + str(n)
                            new_s_dict[state_name] = i
                            new_sp_dict[state_name] = {}
                            new_spe_dict[state_name] = {}
                            for p in [j for j in merge_each[i].data_nest.keys() if s in merge_each[i].data_nest[j].keys()]:
                                ps_dict[p][state_name] = [i]
                                new_sp_dict[state_name][p] = [i]
                                new_spe_dict[state_name][p] = {}
                                for e in merge_each[i].data_nest[p][s].keys():
                                    new_spe_dict[state_name][p][e] = [i]
                                    for row in merge_each[i].data_nest[p][s][e]:
                                        row['State'] = state_name
                            n = n + 1
                    else:   # Merging matches, keep state name
                        new_s_dict[s] = s_dict[s]
                        new_sp_dict[s] = sp_dict[s]
                        new_spe_dict[s] = spe_dict[s]
                        for p in sp_dict[s].keys():
                            ps_dict[p][s] = sp_dict[s][p]

                else:   # No need to change, only one instance
                    new_s_dict[s] = s_dict[s]
                    new_sp_dict[s] = sp_dict[s]
                    new_spe_dict[s] = spe_dict[s]
                    for p in sp_dict[s].keys():
                        ps_dict[p][s] = sp_dict[s][p]

            s_dict = new_s_dict
            sp_dict = new_sp_dict
            spe_dict = new_spe_dict

            new_sp_dict = dict.fromkeys(sp_dict.keys())
            new_spe_dict = dict.fromkeys(spe_dict.keys())
            for s in sp_dict.keys():
                new_sp_dict[s] = {}
                new_spe_dict[s] = {}

            for p in ps_dict.keys():
                if any([len(ps_dict[p][s])>1 for s in ps_dict[p].keys()]):  # Duplicates found
                    if not p_bool:  # If not merging matches, then rename
                        n=0
                        for i in p_dict[p]:
                            protein_name = p + '-'+str(n)
                            for s in merge_each[i].data_nest[p].keys():
                                new_sp_dict[s][protein_name] = i
                                new_spe_dict[s][protein_name] = {}
                                for e in merge_each[i].data_nest[p][s].keys():
                                    new_spe_dict[s][protein_name][e] = [i]
                                    for row in merge_each[i].data_nest[p][s][e]:
                                        row['Protein'] = protein_name
                            n = n + 1
                    else:   # Merging matches, keep protein name
                        for i in p_dict[p]:
                            for state in [s for s in ps_dict[p].keys() if i in ps_dict[p][s]]:
                                new_sp_dict[state][p] = sp_dict[state][p]
                                new_spe_dict[state][p] = spe_dict[state][p]
                else:   # No need to change, only one instance
                    for i in p_dict[p]:
                        for state in [s for s in ps_dict[p].keys() if i in ps_dict[p][s]]:
                            new_sp_dict[state][p] = sp_dict[state][p]
                            new_spe_dict[state][p] = spe_dict[state][p]

            sp_dict = new_sp_dict
            spe_dict = new_spe_dict

            merge_dict = dict.fromkeys(s_dict.keys())
            for s in merge_dict.keys():
                merge_dict[s] = {}

            for s in s_dict.keys():
                for p in sp_dict[s].keys():
                    merge_dict[s][p] = {}
                    for e in [i for i in spe_dict[s][p].keys() if len(spe_dict[s][p][i]) > 0]:
                        merge_dict[s][p][e] = []
                        if len(spe_dict[s][p][e]) > 1:
                            if e_bool:  # Merge deeper
                                merge_dict[s][p][e] = list(spe_dict[s][p][e])
                                print('Exposures being merged, may result in duplicates')

                            else:  # POPUP NOTICE -> Exposures can't be renamed, data will be discarded
                                merge_dict[s][p][e] = [list(spe_dict[s][p][e])[0]]
                                print('Exposures cannot be renamed, Extra data will be discarded')
                        else:
                            merge_dict[s][p][e] = [list(spe_dict[s][p][e])[0]]
            for file in merge_each.keys():
                merge_each[file].organizeData()

            self.main.state_obj[self.name].mergeFiles(merge_dict, merge_each)
            id = self.main.fv.file_list.insert('', 'end', text=self.name, values=(
                self.name, str(', '.join(self.main.state_obj[self.name].proteins)),
                str(', '.join(self.main.state_obj[self.name].states)),
                str(', '.join([str(i) for i in self.main.state_obj[self.name].exposures])),
                self.main.state_obj[self.name].corrected, "No"))
            self.view.top.destroy()
            self.main.view.file_list.selection_set(id)
            self.main.view.file_list.focus(id)
            self.main.select_item()
        files = self.view.file_list.get(0, Tk.END)
        merge_each = {}
        n = 1
        self.name = str(self.view.entry_name.get())
        if self.name in self.main.state_obj.keys():
            while self.name in self.main.state_obj.keys():
                self.name = self.name + '(%s)' % n
                n = n + 1
            self.main.popup('Error!', 'Name renamed to:' + str(self.name))
        self.main.state_obj[self.name] = DECA(empty=True)
        self.view.auto()
        self.view.button_merge_auto_exe.configure(command=lambda:merge())

class Replicates():
    '''
    Replicate Window Controller
    Prompts entry of the numper of replicates per each Exposure-State
    Widgets
        Table of Replicates
        Commit button
    '''
    # Initializes the replicates prompt, adds content to widgets, sets button bindings
    def __init__(self,main):
        self.root = main.root
        self.main = main
        self.view = View.Replicates(self.root)
        main.widgets.append('Replicates')
        main.widget_dict['Replicates'] = self.view.top
        main.view.window_menu.add_command(label='Replicates', command=lambda: main.focus(self.view.top))
        dict_reps = {}
        for state in main.state_obj[main.csvfile].states:
            dict_reps[state] = {}
            for exposure in main.state_obj[main.csvfile].exposures:
                dict_reps[state][exposure] = 3
                self.view.listbox_name.insert('', 'end', values=(state, exposure,dict_reps[state][exposure]))
        self.view.button_apply.configure(command=lambda: self.commit(self.view.listbox_name.list_dict()))
        self.view.top.protocol("WM_DELETE_WINDOW", lambda: main.on_closing(self.view.top, 'Replicates'))
        self.dict_reps = dict_reps

    # Saves the content entered into the replicates window
    def commit(self, cor_dict):
        for i in list(self.dict_reps.keys()):
            for j in list(self.dict_reps[i].keys()):
                self.dict_reps[i][j] = int(cor_dict[str(i)+'-'+str(j)])
        self.main.state_obj[self.main.csvfile].replicates = self.dict_reps
        self.main.on_closing(self.view.top, 'Replicates')

class Correct():
    '''
    Correct Window Controller
    Widgets
        Dropdown for 100%D Control
        Global Back Exchange Entrybox
        Table of exposures and correction factors
    '''
    # Initializes Back Exchange Correction window, fills widgets, binds listbox selection
    def __init__(self, main):
        self.main = main
        self.view = View.Correct(main)
        main.widget_dict['xCorrect'] = self.view.top
        self.main.view.window_menu.add_command(label='xCorrect', command=lambda: self.main.focus(self.view.top))
        dict_exp = {}
        exposures = ['None']
        exposures.extend([str(i) for i in self.main.state_obj[self.main.csvfile].exposures])
        self.view.combobox_exposure.configure(values=exposures)
        self.view.combobox_exposure.set('None')
        self.view.entry_global.insert(0, "None")
        if main.state_obj[main.csvfile].merge_and_corrected == True:
            for i in sorted(main.state_obj[main.csvfile].exposures):
                dict_exp[float(i)] = 'multiple'
            self.view.button_apply['state'] = 'disabled'
            main.popup('Note', 'Back exchange corrected files were merged.\nFurther correction cannot be performed')
        else:
            if self.main.corrected == 'Yes':
                dict_exp = self.main.state_obj[self.main.csvfile].cor_factor
            else:
                for i in sorted(self.main.state_obj[self.main.csvfile].exposures):
                    dict_exp[float(i)] = 1
        for i in sorted(dict_exp.keys()):
            self.view.listbox_name.insert('', 'end', values=(i, dict_exp[i]))
        self.view.listbox_name.bind('<<ListboxSelect>>',lambda e: self.select())
        self.view.button_apply.configure(command = lambda: self.commit(self.view.listbox_name.list_dict()))
        self.view.button_help.configure(command = lambda: self.main.help('2) Back-Exchange Correction'))
        self.view.top.protocol("WM_DELETE_WINDOW", lambda: main.on_closing(self.view.top, 'xCorrect'))

    # Prints selected item in listbox
    def select(self):
        print(self.view.listbox_name.list_dict())

    # Saves Correction factors entered into listbox
    def commit(self, cor_dict):
        new_dict = {}
        for item in list(cor_dict.keys()):
            new_dict[float(item)] = cor_dict[item]
        global_bx = self.view.entry_global.get()
        timepoint = self.view.combobox_exposure.get()
        if timepoint == 'None':
            timepoint = None
        if global_bx == 'None':
            global_bx = None
        self.main.state_obj[self.main.csvfile].backCorrect(new_dict,global_bx,timepoint)
        self.main.view.file_list.set(self.main.curitem, 'xCorrected', 'Yes')
        self.main.on_closing(self.view.top, 'xCorrect')
        self.main.select_item()

class Recombine():
    '''
    Recombine Window Controller
    Widgets
        Text describing peptide increase
        Commit button
    '''
    # Initialize Recombination window, perform initial segmentation survey, apply button bindings
    def __init__(self, main):
        self.main = main
        if main.state_obj[self.main.csvfile].recombined:
            main.popup('Error!', 'Recombination can only be applied once.')
            return
        self.view = View.Recombine(main)
        main.widget_dict['Recombine'] = self.view.top
        self.main.view.window_menu.add_command(label='Recombine', command=lambda: self.main.focus(self.view.top))
        self.main.state_obj[self.main.csvfile].segmentPeptides()
        for state in self.main.state_obj[self.main.csvfile].states:
            for protein in self.main.state_obj[self.main.csvfile].proteins:
                if state in self.main.state_obj[self.main.csvfile].data_nest[protein].keys():
                    n = 0
                    tot_size = 0
                    for pep in self.main.state_obj[self.main.csvfile].recombine_dict[protein][state][0]:
                        n = n + 1
                        tot_size = tot_size + len(pep.get('Sequence', 0))
                    if n != 0 and tot_size != 0:
                        ave_size = tot_size / n
                    else:
                        ave_size = 0
                    self.view.text.insert(Tk.END, str(n) + " new peptides at " + str(ave_size) +
                                                  " ave. size for " + str(protein) + " in " + str(state) + '\n')
        self.view.button_apply.configure(command = lambda: self.commit())
        self.view.button_help.configure(command = lambda: self.main.help('3) Peptide Recombination'))
        self.view.text.config(state=Tk.DISABLED)
        self.view.top.protocol("WM_DELETE_WINDOW", lambda: main.on_closing(self.view.top, 'Recombine'))

    # Commit newly segmented peptides to the dataset
    def commit(self):
        self.main.state_obj[self.main.csvfile].applySegmentation()
        self.view.top.destroy()
        self.main.view.file_list.set(self.main.curitem, 'Recombined', 'Yes')
        self.main.select_item()
        self.main.on_closing(self.view.top, 'Recombine')

class Significance():
    '''
    Significance Window Controller
    Widgets
        Table of Statistics
        Exposure dropdown
        Plot of confidence intervals
    '''

    # Initialize Significance Window, Bind Exposure selection, Start right click menu
    def __init__(self,main):
        print('initiating sig')
        if main.state_obj[main.csvfile].filetype == 'CSV':
            main.rc = Replicates(main)
        # Initialize Viewer
        main.widgets.append('Significance')
        self.main = main
        self.view = View.Significance(main.root)

        main.widget_dict['Significance'] = self.view.top
        self.main.view.window_menu.add_command(label='Significance', command=lambda: self.main.focus(self.view.top))
        self.view.top.protocol("WM_DELETE_WINDOW", lambda: main.on_closing(self.view.top, 'Significance'))

        # Bind radiobutton change to Canvas update
        self.view.combo_exp.bind('<<ComboboxSelected>>', lambda event: self.update_data(True))
        self.stats = None
        self.ax = None
        self.exp = None

        self.view.popupMenu.add_command(label='Save', command=lambda: self.export())

        if 'Data' in main.widgets:
            if len(main.data.view.tv.selection())>0:
                self.update_data()

    # Save plot to file
    def export(self):
        savefile = FileDialog.asksaveasfilename(
            title="Choose save location",
            initialfile=self.protein + "_" + str(self.peptide['Start']) + "-" + str(self.peptide['End']) + "_Sig.png",
            initialdir=expanduser('~'),
            filetypes=(("PNG", '*.png'), ("all files", "*.*")))
        self.view.figure.savefig(savefile, dpi=300, transparent=True)

    # Update content in widgets
    def update_data(self, selection=False):
        self.view.ax.clear()
        self.view.plot.draw()
        #self.view.plot.get_tk_widget().update()
        self.view.tukeytable.delete(*self.view.tukeytable.get_children())
        if self.main.data.view.tv.focus() in self.main.data.xref.keys():
            curitem = self.main.data.xref[self.main.data.view.tv.focus()]
        else:
            print('not in keys')
            return

        # get states

        self.states = list(set([i['State'] for i in curitem]))

        # get protein
        self.protein = curitem[0]['Protein']

        # get exposures
        exposures = sorted(list(set([i['Exposure'] for i in curitem])))
        exposures.remove(0)

        if selection:
            exposure = float(self.view.combo_exp.get())
            self.exp = exposure
        elif self.exp is not None:
            exposure = self.exp
            self.view.combo_exp.set(exposure)
        else:
            exposure=float(exposures[-1])
            self.exp = exposure
            self.view.combo_exp.set(exposure)

        peptide = curitem[0]

        # Fill combobox with exposures

        self.view.combo_exp.configure(values=sorted(exposures))
        self.view.combo_exp.set(exposure)
        if not selection:
            self.main.state_obj[self.main.csvfile].getSignificance(peptide['PepID'])
            self.stats = self.main.state_obj[self.main.csvfile].stats_nest[peptide['PepID']]
        if self.stats[exposure]['Stats'] != {}:
            # Fill Levene Text
            if 'Levene' not in self.stats[exposure]['Stats'].keys():
                self.view.levene.configure(text='N/A')
            elif float(self.stats[exposure]['Stats']['Levene']['p']) < 0.05:
                self.view.levene.configure(text=self.stats[exposure]['Stats']['Levene']['p'],fg='red')
            else:
                self.view.levene.configure(text=self.stats[exposure]['Stats']['Levene']['p'],fg='green')

            # Fill Anova Text
            if 'Anova' not in self.stats[exposure]['Stats'].keys():
                self.view.levene.configure(text='N/A')
            elif float(self.stats[exposure]['Stats']['Anova']['p']) < 0.05:
                self.view.anova.configure(text=self.stats[exposure]['Stats']['Anova']['p'],fg='green')
            else:
                self.view.anova.configure(text=self.stats[exposure]['Stats']['Anova']['p'],fg='red')


            # Fill Tukey Table
            if 'Tukey' in self.stats[exposure]['Stats'].keys():
                columns = ('group1', 'group2', 'meandiff', 'lower', 'upper', 'RejectNull', 'p-value')
                self.view.tukeytable.configure(columns=columns, show='headings', height=6)
                for col in columns:
                    self.view.tukeytable.column(col, width='100', stretch=0)
                    self.view.tukeytable.heading(col, text=col)

                for i in list(self.stats[exposure]['Stats']['Tukey']):
                    data = [i[key] for key in ['group1', 'group2', 'meandiff', 'lower', 'upper', 'reject', 'p-value']]
                    self.view.tukeytable.insert('', 'end', values=data)

                # Add TukeyHSD Plot to Canvas
                if self.main.state_obj[self.main.csvfile].filetype == 'DNX':
                    self.stats[exposure]['Stats']['Tukey_obj'].plot_simultaneous(ax=self.view.ax)
                    self.view.ax.clear()
                    groups = self.stats[exposure]['Stats']['Tukey_obj'].groupsunique
                    halfwidths = self.stats[exposure]['Stats']['Tukey_obj'].halfwidths
                    means = self.stats[exposure]['Stats']['Tukey_obj']._multicomp.groupstats.groupmean
                    self.plot_simultaneous(groups, means, halfwidths, self.view.ax)
                else:
                    self.view.ax.clear()
                    groups = [i for i in self.stats[exposure].keys() if i != 'Stats']
                    means = [self.stats[exposure][i][0] for i in groups]
                    halfwidths = [self.stats[exposure]['Stats']['CI'][i] for i in groups]
                    self.plot_simultaneous(groups, means, halfwidths, self.view.ax)
                self.peptide = peptide
                self.view.figure.set_size_inches((8, 6))
                self.view.plot.draw()
                self.view.plot.get_tk_widget().update()
            elif 'Ttest' in self.stats[exposure]['Stats'].keys():
                columns = ('group1', 'group2', 'meandiff', 'RejectNull', 'p-value')
                self.view.tukeytable.configure(columns=columns, show='headings', height=6)
                for col in columns:
                    self.view.tukeytable.column(col, width='100', stretch=0)
                    self.view.tukeytable.heading(col, text=col)
                data = [self.stats[exposure]['Stats']['Ttest'][key] for key in ['group1', 'group2', 'meandiff', 'reject', 'p-value']]
                self.view.tukeytable.insert('', 'end', values=data)
                self.view.ax.clear()
                groups = [i for i in self.stats[exposure].keys() if i != 'Stats']
                means = [self.stats[exposure][i][0] for i in groups]
                halfwidths = [self.stats[exposure]['Stats']['CI'][i] for i in groups]
                self.plot_simultaneous(groups, means, halfwidths, self.view.ax)
                self.view.figure.set_size_inches((8, 6))
                self.view.plot.draw()
                self.view.plot.get_tk_widget().update()
            else:
                self.view.ax.clear()
        self.peptide = peptide

    # Unused
    def tukey_heat(self, axis):
        states1 = self.main.state_obj[self.main.csvfile].stats_nest[32][2.0]['Stats']['Tukey_obj'].groupsunique
        states2 = self.main.state_obj[self.main.csvfile].stats_nest[32][2.0]['Stats']['Tukey_obj'].groupsunique

        harvest = np.array([[0.9, 0.118, 0.005, 0.001],
                            [0.118, 0.9, 0.177, 0.011],
                            [0.005, 0.177, 0.9, 0.248],
                            [0.001, 0.011, 0.248, 0.9]])

        axis.imshow(harvest, cmap='bwr')

        axis.set_xticks(np.arange(len(states1)))
        axis.set_yticks(np.arange(len(states2)))
        axis.set_xticklabels(states1)
        axis.set_yticklabels(states2)

        # Rotate the tick labels and set their alignment.
        plt.setp(axis.get_xticklabels(), rotation=45, ha="right",
                 rotation_mode="anchor")

        # Loop over data dimensions and create text annotations.
        for i in range(len(states1)):
            for j in range(len(states2)):
                axis.text(j, i, harvest[i, j], ha="center", va="center", color="w")

        axis.set_title("ASB9 78-84 2min pvalues")

    # Unused
    def tukey_hsd(self, axis):
        pass

    # Function from Statsmodels, plot confidence intervals
    def plot_simultaneous(self, groups, means, halfwidths, axis):
        axis.cla()
        minrange = [means[i] - halfwidths[i] for i in range(len(means))]
        maxrange = [means[i] + halfwidths[i] for i in range(len(means))]
        for i in range(0,len(means)):
            axis.errorbar(means[i], list(range(0, len(means)))[i], xerr=[i / 2 for i in halfwidths][i],
                         marker='o', linestyle='None', markeredgecolor=self.main.set_colors[groups[i]],
                          ecolor=self.main.set_colors[groups[i]])
        axis.set_title('Confidence Intervals', fontname='Arial', fontsize = 24)
        axis.set_xlabel('Deuterium Uptake (Da)', fontname='Arial', fontsize=20)
        r = max(maxrange) - min(minrange)
        axis.set_ylim([-1, len(groups)])
        axis.set_xlim([min(minrange) - r / 10., max(maxrange) + r / 10.])
        axis.set_yticklabels(np.insert(groups, 0, ''))
        axis.set_yticks(np.arange(-1, len(means) + 1))
        for tick in axis.get_xticklabels():
            tick.set_fontname('Arial')
            tick.set_fontsize(16)
        for tick in axis.get_yticklabels():
            tick.set_fontname('Arial')
            tick.set_fontsize(16)

class Spectrum():
    '''
    Spectrum Window Controller
    Widgets:
        Dropdowns for Exposure, State, Rawfile, Charge
        Spectrum plot
    '''
    # Initializes Spectrum Window, Right click menu
    def __init__(self,main):
        # Initialize Viewer

        self.main = main
        self.view = View.Spectrum(main.root)
        main.widgets.append('Spectrum')
        main.widget_dict['Spectrum'] = self.view.top
        self.main.view.window_menu.add_command(label='Spectrum', command=lambda: self.main.focus(self.view.top))
        self.view.top.protocol("WM_DELETE_WINDOW", lambda: main.on_closing(self.view.top, 'Spectrum'))

        #self.main.state_obj[self.main.csvfile].findOutliers()
        if 'Data' in main.widgets:
            if len(main.data.view.tv.selection())>0:
                self.update()
        self.view.popupMenu.add_command(label='Save', command=lambda: self.export())

    # Save plot to file
    def export(self):
        savefile = FileDialog.asksaveasfilename(
            title="Choose save location",
            initialfile=self.protein + "_S:" + str(self.state) + "_E:" + str(self.exposure) + "_F:" +
                        str(self.raw) + "_Z:" + str(self.z) + "_" + str(self.peptide['Start']) + "-" +
                        str(self.peptide['End']) + "_Ions.png",
            initialdir=expanduser('~'),
            filetypes=(("PNG", '*.png'), ("all files", "*.*")))
        self.view.figure.savefig(savefile, dpi=300, transparent=True)

    # Update widget content
    def update(self):
        if self.main.data.view.tv.focus() in self.main.data.xref.keys():
            curitem = self.main.data.xref[self.main.data.view.tv.focus()]
        else:
            print('not in keys')
            return

        # get proteins
        self.protein = curitem[0]['Protein']

        # get states
        self.states = list(set([i['State'] for i in curitem]))

        # get exposures
        exposures = sorted(list(set([i['Exposure'] for i in curitem])))

        peptide = curitem[0]

        # get charges
        charges = list(set([i['z'] for i in self.main.state_obj[self.main.csvfile].getDnxCluster(['PepID = '+str(peptide['PepID'])])]))

        #get ions
        ion_nest = {}
        for state in self.states:
            ion_nest[state] = {}
            for exposure in exposures:
                ion_nest[state][exposure]={}
                rawfiles = set([i['RawID'] for i in self.main.state_obj[self.main.csvfile].rawfiles if
                                (i['State'] == state) and
                                (i['Exposure'] == exposure)])
                for raw in rawfiles:
                    ion_nest[state][exposure][raw] = charges
        self.ion_nest=ion_nest

        # Fill ScrolledListboxes and set defaults

        self.view.combo_state.configure(values=self.states)
        self.view.combo_state.set(self.states[0])

        self.view.combo_exp.configure(values=exposures)
        self.view.combo_exp.set(exposures[0])

        self.view.combo_raw.configure(values=list(ion_nest[self.states[0]][exposures[0]].keys()))
        self.view.combo_raw.set(list(ion_nest[self.states[0]][exposures[0]].keys())[0])

        self.view.combo_z.configure(values=charges)
        self.view.combo_z.set(charges[0])

        # Bind Scrollbox changes to plot update

        self.view.combo_state.bind('<<ComboboxSelected>>', lambda event: self.plot('state'))
        self.view.combo_raw.bind('<<ComboboxSelected>>', lambda event: self.plot('raw'))
        self.view.combo_exp.bind('<<ComboboxSelected>>', lambda event: self.plot('exp'))
        self.view.combo_z.bind('<<ComboboxSelected>>', lambda event: self.plot('z'))
        self.pepid = peptide['PepID']
        self.peptide = peptide
        self.plot('z')

    # Plot spectra
    def plot(self,trigger):
        if trigger == 'state':
            state = self.view.combo_state.get()
            exposures = list(self.ion_nest[state].keys())
            if float(self.view.combo_exp.get()) in exposures:
                exposure = float(self.view.combo_exp.get())
            else:
                exposure = exposures[0]
            raws = list(self.ion_nest[state][exposure].keys())
            raw = raws[0]
            charges = self.ion_nest[state][exposure][raw]
            z = charges[0]
            self.view.combo_exp.configure(values=exposures)
            self.view.combo_exp.set(exposure)
            self.view.combo_raw.configure(values=raws)
            self.view.combo_raw.set(raw)
            self.view.combo_z.configure(values=charges)
            self.view.combo_z.set(z)

        if trigger == 'exp':
            state = self.view.combo_state.get()
            exposure = float(self.view.combo_exp.get())
            raws = list(self.ion_nest[state][exposure].keys())
            raw = raws[0]
            charges = self.ion_nest[state][exposure][raw]
            z = charges[0]
            self.view.combo_raw.configure(values=raws)
            self.view.combo_raw.set(raw)
            self.view.combo_z.configure(values=charges)
            self.view.combo_z.set(z)

        if trigger == 'raw':
            state = self.view.combo_state.get()
            exposure = float(self.view.combo_exp.get())
            raw = float(self.view.combo_raw.get())
            charges = self.ion_nest[state][exposure][raw]
            z = charges[0]
            self.view.combo_z.configure(values=charges)
            self.view.combo_z.set(z)

        if trigger == 'z':
            state = self.view.combo_state.get()
            exposure = float(self.view.combo_exp.get())
            raw = float(self.view.combo_raw.get())
            z = float(self.view.combo_z.get())

        self.state = state
        self.exposure = exposure
        self.raw = raw
        self.z = z

        limits = [(i['mzCenterStart'], i['mzCenterEnd']) for i in self.main.state_obj[self.main.csvfile].getDnxCluster(['PepID = '+str(self.pepid),'a.chargeState = '+str(z)])][0]
        commands = ['PepID = '+str(self.pepid),'b.chargeState = '+str(z),'RawID = '+str(raw)]
        ions = [(i['mz'], i['Int'], i['RT'], i['Drift']) for i in self.main.state_obj[self.main.csvfile].getDnxIons(commands)]
        self.main.state_obj[self.main.csvfile].findOutliers(pep=self.pepid)
        obj = self.main.state_obj[self.main.csvfile].file_nest[self.pepid]['States'][state][exposure][raw][z]
        outliers = [obj['data'][obj['ppm']['data'].index(i)] for i in obj['ppm']['data'] if i>10]
        outions = [i for i in ions if i[0] in outliers]
        ions = [i for i in ions if i not in outions]
        nullions_full = self.main.state_obj[self.main.csvfile].getDnxIons(
            commands=['mz BETWEEN {} AND {}'.format(limits[0],limits[1]),
                      'rawFileID == {}'.format(raw),
                      'Int > 10000',
                      'RT BETWEEN {} AND {}'.format(self.peptide['RT']-0.05*self.peptide['RT'],self.peptide['RT']+0.05*self.peptide['RT'])],null=True)
        nullions = [(i['mz'],i['Int']) for i in nullions_full]

        self.view.ax.clear()
        self.view.ax.bar([i[0] for i in nullions], [i[1] for i in nullions], width=(limits[1] - limits[0]) * 0.002, color='black')
        maxint = 0
        if len(ions)>0:
            maxint = max([float(i[1]) for i in ions])
            self.view.ax.bar([i[0] for i in ions], [i[1] for i in ions], width=(limits[1] - limits[0]) * 0.004,
                             color='blue')
            self.view.ax.bar([i[0] for i in outions], [i[1] for i in outions], width=(limits[1] - limits[0]) * 0.004,
                             color='red')
        elif len(nullions)>0 and max([float(i[1]) for i in nullions])>maxint:
            maxint = max([float(i[1]) for i in nullions])
        else:
            maxint=0

        self.view.ax.set_xlim(limits[0],limits[1])
        self.view.ax.set_ylim(0, maxint*1.1)
        self.view.plot.draw()
        return

class Stats():
    '''
    Outlier Statistics Window Controller
    Widgets:
        Ion Statistics Table
    '''

    # Initialize Stats Window, Gets data if peptide selected
    def __init__(self, main):
        self.main = main
        self.main.widgets.append('Stats')
        self.view = View.Stats(self.main.root)
        self.main.widget_dict['Stats'] = self.view.top
        self.main.view.window_menu.add_command(label='Stats', command=lambda: self.main.focus(self.view.top))
        self.view.top.protocol("WM_DELETE_WINDOW", lambda: self.main.on_closing(self.view.top, 'Stats'))
        self.peplist = self.main.state_obj[self.main.csvfile].peplist
        if 'Data' in self.main.widgets:
            self.getData()

    # Retrieves outlier data
    def getData(self):
        if self.main.data.view.tv.focus() in self.main.data.xref.keys():
            curitem = self.main.data.xref[self.main.data.view.tv.focus()]
            self.data = self.main.state_obj[self.main.csvfile].findOutliers(curitem[0]['PepID'])
            self.insertData()

    # Calculates summary statistics for ions
    def calculateStats(self,data):
        rt_mean = mean(data['rt']['data'])
        mob_5pc = []
        mob_nums = []
        mob_means = []
        mob_stds = []
        for z in [key for key in data['drift']['data'].keys() if key != '_STATS_']:
            mob_mean = mean(data['drift']['data'][z])
            mob_means.append(mob_mean)
            mob_stds.append(std(data['drift']['data'][z]))
            mob_nums.append(len(data['drift']['data'][z]))
            mob_5pc.append(len([i for i in data['drift']['data'][z] if abs(i - mob_mean) > (mob_mean * 0.05)]))
        denom = sum([n for n in mob_nums])
        num = sum([n * s ** 2 for n, s in zip(mob_nums, mob_stds)]) + sum(
            [n * (mean(mob_means) - m) ** 2 for n, m in zip(mob_nums, mob_means)])
        std_all = sqrt(num / denom)
        pc_all = sum(mob_5pc) / sum(mob_nums)
        d = [mean(data['ppm']['data']),
             std(data['ppm']['data']),
             len([i for i in data['ppm']['data'] if i > 10]) / len(data['ppm']['data']),
             mean(data['rt']['data']),
             std(data['rt']['data']),
             len([i for i in data['rt']['data'] if abs(i - rt_mean) > (rt_mean * 0.05)]) / len(data['rt']['data']),
             tuple(mob_means),
             tuple(mob_stds),
             pc_all]
        for i in [0, 1, 2, 3, 4, 5, 8]:
            d[i] = round(d[i], 3)
            if d[i] == 0:
                d[i] = ''
        d[6] = ','.join(map(str, [round(i,3) for i in d[6]]))
        d[7] = ','.join(map(str, [round(i, 3) for i in d[7]]))
        return d

    # Adds content to table
    def insertData(self):
        self.view.tv.delete(*self.view.tv.get_children())
        for pep in self.data.keys():
            data = self.data[pep]['_STATS_']
            d = self.calculateStats(data)
            pid = self.view.tv.insert('', 'end', text=self.peplist[pep]['Sequence'], values=tuple(d),open=False)
            for state in [key for key in self.data[pep]['States'].keys() if key != '_STATS_']:
                data = self.data[pep]['States'][state]['_STATS_']
                d = self.calculateStats(data)
                sid = self.view.tv.insert(pid, 'end', text=state,
                                          values=tuple(d),
                                          open=True)
                for exp in [key for key in self.data[pep]['States'][state].keys() if key != '_STATS_']:
                    data = self.data[pep]['States'][state][exp]['_STATS_']
                    d = self.calculateStats(data)
                    eid = self.view.tv.insert(sid, 'end', text=exp,
                                              values=tuple(d),
                                              open=True)
                    for raw in [key for key in self.data[pep]['States'][state][exp].keys() if key != '_STATS_']:
                        data = self.data[pep]['States'][state][exp][raw]['_STATS_']
                        d = self.calculateStats(data)
                        rid = self.view.tv.insert(eid, 'end', text=raw,
                                                  values=tuple(d),
                                                  open=True)
                        for charge in [key for key in self.data[pep]['States'][state][exp][raw].keys() if key != '_STATS_']:
                            data = self.data[pep]['States'][state][exp][raw][charge]
                            rt_mean = mean(data['rt']['data'])
                            mob_mean = mean(data['drift']['data'])
                            mob_5pc = (len([i for i in data['drift']['data'] if
                                                abs(i - mob_mean) > (mob_mean * 0.05)]))/len(data['drift']['data'])
                            d = [mean(data['ppm']['data']),
                                 std(data['ppm']['data']),
                                 len([i for i in data['ppm']['data'] if i > 10]) / len(data['ppm']['data']),
                                 mean(data['rt']['data']),
                                 std(data['rt']['data']),
                                 len([i for i in data['rt']['data'] if abs(i - rt_mean) > (rt_mean * 0.05)]) / len(
                                     data['rt']['data']),
                                 mean(data['drift']['data']),
                                 std(data['drift']['data']),
                                 mob_5pc]
                            for i in [0, 1, 2, 3, 4, 5, 6, 7, 8]:
                                d[i] = round(d[i], 3)
                                if d[i] == 0:
                                    d[i] = ''
                            self.view.tv.insert(rid, 'end', text=charge,
                                                      values=tuple(d),
                                                      open=True)

class Plot():
    '''
    Plot Window Controller
    Widgets:
        Plot Canvas
        Reset Settings, Save Figure, Plot Options Buttons
    '''
    # Initialize window, define initial settings, define Tk variables, draw axis
    def __init__(self, main):
        self.main = main
        self.view = View.Plot(main)  # Initialize Plot GUI
        main.widget_dict['Plot'] = self.view.top
        main.widgets.append('Plot')
        self.main.view.data_menu.entryconfig('Export Plot', state='normal')
        self.main.view.window_menu.add_command(label='Plot', command=lambda: self.main.focus(self.view.top))

        self.title_val = Tk.IntVar()
        self.xaxis_val = Tk.IntVar()
        self.yaxis_val = Tk.IntVar()
        self.ann_val = Tk.IntVar()
        self.leg_val = Tk.IntVar()
        self.title_val.set(1)
        self.xaxis_val.set(1)
        self.yaxis_val.set(1)
        self.ann_val.set(1)
        self.leg_val.set(1)
        self.position = 0
        self.linestyles = ['solid','dotted','dashed','dashdot']
        self.markers = {'circle':'o','triangle':'^','square':'s','plus':'p','star':'*','diamond':'D','x':'X'}
        self.state_markers = {}
        self.state_linestyles = {}
        self.listbox_dict = {}
        self.last_states = 0
        self.view.top.protocol("WM_DELETE_WINDOW", lambda: self.on_close())
        self.default_options = {'font':'Arial', 'ls':'5', 'ms':'10', 'as':'3', 'ar':1,
                                'leg': 1, 'ann': 1, 'title': 1, 'xaxis': 1, 'yaxis': 1,
                                'title_size': 40, 'xaxis_size': 32, 'yaxis_size': 32, 'ann_size': 28, 'leg_size': 16,
                                'xaxis_labels': 'Fitted', 'yaxis_labels': '',
                                'xaxis_title':'Time (min)','yaxis_title':'Deuterium Uptake (Da)',
                                'colors': {}, 'markers': {}, 'linestyles':{}, 'fit':'lin', 'yaxis_units':'Da'}
        self.settings = copy(self.default_options)
        self.xref = {}
        self.initialize_axis()

    def on_close(self):
        self.main.on_closing(self.view.top, 'Plot')
        self.main.view.data_menu.entryconfig('Export Plot', state='disabled')
    # Draws Axis
    def initialize_axis(self):
        self.view.cm = self.view.context_menu(self.view.figcan.get_tk_widget(), self)
        self.view.cm.popupMenu.add_command(label='Combine', command=lambda: self.select_item_data(combine=True))
        self.view.cm.popupMenu.add_command(label='Advanced', command=lambda: self.options())
        self.view.cm.popupMenu.add_command(label='Reset', command=lambda: self.reset())
        self.view.cm.popupMenu.add_command(label='Save', command=lambda: self.save())
        self.view.cm.popupMenu.entryconfig('Combine', state='disabled')
        axis = self.view.figure.add_subplot(111, aspect='auto')
        axis.set_xscale('linear')
        axis.set_xbound((0, 10))
        axis.set_xlim((0, 10))
        axis.set_ybound((0, 10))
        axis.set_ylim((0, 10))
        axis.set_xlabel('Exposure (min)', fontname=self.settings['font'],
                        fontsize=self.settings['xaxis_size'])
        axis.set_ylabel('Deuterium Incorporated (Da)',
                        fontname=self.settings['font'],
                        fontsize=self.settings['yaxis_size'])
        axis.set_title('Select Peptide',
                       fontname=self.settings['font'],
                       fontsize=self.settings['title_size'])
        self.view.figcan.draw()
        self.axis = axis

    # Updates Canvas, maintains aspect ratio
    def update_canvas(self):
        self.view.canvas.update()
        h = self.view.canvas.winfo_height()
        w = self.view.canvas.winfo_width()
        if h / w > float(self.settings['ar']):
            self.view.figcan.get_tk_widget().config(height=w * self.settings['ar'], width=w)
        elif h/w < float(self.settings['ar']):
            self.view.figcan.get_tk_widget().config(height=h, width=h/self.settings['ar'])
        else:
            return
        self.view.figcan.draw()
        self.view.figcan.get_tk_widget().pack(ipadx=20,ipady=20)
        self.view.figcan.get_tk_widget().update()

    # Draws new peptide when peptide selected
    def select_item_data(self, combine=False):
        self.view.figure.clear()
        self.axis.clear()
        self.view.figcan.draw()
        self.view.figcan.get_tk_widget().update()

        if combine:
            self.settings['yaxis_units'] = 'percentDa'
            merged = {}
            for item in self.selected_peptides:
                for key in [i for i in item.keys() if i != 'options']:
                    if key in list(merged.keys()):
                        print('Duplicate State Skipped')
                    else:
                        merged[key] = item[key]
            merged['options'] = {'leg_bbox': '', 'ann_bbox': ''}
            self.selected_peptides = [merged]
            self.xref['Double'] = merged
        else:
            self.selected_peptides = []
            for item in self.main.data.view.tv.selection():
                print(item)
                if (item in self.main.data.xref.keys()) and (item not in self.xref.keys()):
                    print('New Pep')
                    peptides = self.main.data.xref[item]
                    pep = {}
                    pep['options'] = {}
                    pep['options']['leg_bbox'] = ''
                    pep['options']['ann_bbox'] = ''
                    n = 0
                    for state in list(set([i['State'] for i in peptides])):
                        pep[state] = [i for i in peptides if i['State'] == state]
                    self.selected_peptides.append(pep)
                    self.xref[item] = pep
                elif (item in self.main.data.xref.keys()) and (item in self.xref.keys()):
                    print('Old Pep')
                    self.selected_peptides.append(self.xref[item])

        if len(self.selected_peptides)>0:
            self.states = [i for i in self.selected_peptides[0].keys() if i != 'options']
            self.protein = self.selected_peptides[0][self.states[0]][0]['Protein']
            n = 0
            self.state_dict = {}
            for state in self.states:
                self.state_dict[str("State" + str(n))] = state
                n = n + 1

            self.update_canvas()
            self.axis = self.view.figure.add_subplot(111, aspect='auto')
            print(self.selected_peptides[0]['options'])
            self.uptake(self.protein,self.selected_peptides[0],self.state_dict,self.axis)
            self.view.figcan.draw()
            pep_range = []
            for pep in self.selected_peptides:
                pep_state = list([i for i in pep.keys() if i != 'options'])[0]
                pep_range.append((pep[pep_state][0]['Start'],pep[pep_state][0]['End']))
            if len(self.selected_peptides)>=2 and len(set(pep_range)) == 1:
                self.view.cm.popupMenu.entryconfig('Combine', state='normal')
        if combine:
            self.settings['yaxis_units'] = 'Da'
            self.view.cm.popupMenu.entryconfig('Combine', state='disabled')

    # Resets settings to defaults
    def reset(self):
        self.settings = copy(self.default_options)
        self.main.set_colors = {}
        colors = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black']
        color_cycle = cycle(colors)
        for state in sorted(self.main.state_obj[self.main.csvfile].states):
            self.main.set_colors[state] = next(color_cycle)
        for i in self.xref.keys():
            self.xref[i]['options'] = {'leg_bbox': '', 'ann_bbox': ''}
        self.title_val.set(1)
        self.xaxis_val.set(1)
        self.yaxis_val.set(1)
        self.ann_val.set(1)
        self.leg_val.set(1)
        self.select_item_data()

    # Plots data onto canvas widget
    def uptake(self, protein, pep, state_args, axis, export = 0):
        plt.ioff()
        ha_val = 'center'
        va_val = 'center'
        def lin_fit(x, t, y):
            return x[0] * (1 - np.exp(-x[1] * t)) + x[2] * (1 - np.exp(-0.01 * t)) - y

        def lin_fit_x1(x, t, y, x0):
            return x0[0] * (1 - np.exp(-x * t)) + x0[2] * (1 - np.exp(-0.01 * t)) - y

        def lin_fit_x0x2(x, t, y, x0):
            return x[0] * (1 - np.exp(-x0 * t)) + x[2] * (1 - np.exp(-0.01 * t)) - y

        def gen_lin_fit(x, t):
            return x[0] * (1 - np.exp(-x[1] * t)) + x[2] * (1 - np.exp(-0.01 * t))

        def log_fit(x, t, y):
            return x[0] * t + x[1] - y

        def gen_log_fit(x, t):
            return x[0] * t + x[1]

        def onpick(event):
            item = event.artist
            if item in legline_dict.keys():
                origline = legline_dict[item]
                if item.get_color() in colors:
                    new_color = next(islice(cycle(colors), colors.index(item.get_color()) + 1, None))
                else:
                    new_color = next(islice(cycle(colors), 0, None))
                origline[0].set_color(new_color)
                origline[1][0].set_color(new_color)
                for i in origline[1][1]:
                    i.set_color(new_color)
                for i in origline[1][2]:
                    i.set_color(new_color)
                item.set_color(new_color)
                item._legmarker.set_color(new_color)
                self.main.set_colors[origline[0].get_label()] = new_color
                self.settings['colors'] = {}
                for key in legline_dict.keys():
                    self.settings['colors'][legline_dict[key][0].get_label()] = key.get_color()
                self.view.figcan.draw()

        def on_release(event, ann, leg):
            'on button press we will see if the mouse is over us and store some data'
            legbox = leg.get_window_extent()
            legbox_data = axis.transAxes.inverted().transform((legbox.x0,legbox.y0))
            annbox = ann.get_window_extent()
            e = (event.x,event.y)
            if (legbox.x0<e[0]<legbox.x1) and (legbox.y0<e[1]<legbox.y1):
                pep['options']['leg_bbox'] = (legbox_data[0],legbox_data[1])
            elif (annbox.x0<e[0]<annbox.x1) and (annbox.y0<e[1]<annbox.y1):
                pep['options']['ann_bbox'] = (ann.xyann, ha_val, va_val)

        def update_annot(ind,i):
            x, y = i[0].get_data()
            annot.xy = (x[ind["ind"][0]], y[ind["ind"][0]])
            name = names[annot.xy]
            text = "{}, {}".format(annot.xy,name)
            annot.set_text(text)
            annot.get_bbox_patch().set_alpha(0.4)

        def hover(event):
            vis = annot.get_visible()

            if event.inaxes == axis:
                for i in axis.containers:
                    cont, ind = i[0].contains(event)
                    if cont:
                        update_annot(ind,i)
                        annot.set_visible(True)
                        self.view.figcan.draw_idle()
                    else:
                        if vis:
                            annot.set_visible(False)
                            self.view.figcan.draw_idle()

        self.resizing = 0

        def on_resize(event):
            print('resize detected')
            self.resizing = 1
            self.view.figcan.get_tk_widget().pack_forget()

        def done_resizing(event):
            if self.resizing:
                self.resizing = 0
                print('resize finished')
                self.update_canvas()

        rcParams['pdf.fonttype'] = 42
        rcParams['ps.fonttype'] = 42
        colors = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black']
        color_cycle = cycle(colors)
        # Prep Data
        default_pep = pep[state_args['State0']][0]
        pep_mass = default_pep['MHP']
        if default_pep['Fragment'] != '':
            pep_mass = default_pep['Fragment']+'-terminal Fragment'
        pep_sequence = default_pep['Sequence']
        protein_name = protein
        protein_range = str(default_pep['Start']) + "-" + str(default_pep['End'])
        max_uptake = float(default_pep['MaxUptake'])
        lines = []
        errs = []
        if self.settings['colors'] != {}:
            for item in self.settings["colors"].keys():
                self.main.set_colors[item] = str(self.settings["colors"][item])
        if self.settings['markers'] != '':
            for item in self.settings["markers"].keys():
                self.state_markers[item] = str(self.settings["markers"][item])
        if self.settings['linestyles'] != '':
            for item in self.settings["linestyles"].keys():
                self.state_linestyles[item] = str(self.settings["linestyles"][item])
        fontsetting = matplotlib.font_manager.FontProperties(family=self.settings['font'],
                                                             weight='normal',
                                                             style='normal', size=self.settings['leg_size'])
        names = {}
        self.view.cm.statMenu.delete(0,self.last_states)
        self.last_states = len(state_args.values())
        for state in sorted(state_args.values()):
            menu = Tk.Menu(self.view.cm.statMenu, tearoff=0)
            self.view.cm.statMenu.add_cascade(label=state, menu=menu)
            cmenu = Tk.Menu(menu, tearoff=0)
            menu.add_cascade(label='colors', menu=cmenu)
            for color in colors:
                cmenu.add_command(label=color, command=lambda state=state, color=color: self.view.cm.update('colors', (state,color)))
            mmenu = Tk.Menu(menu, tearoff=0)
            menu.add_cascade(label='markers', menu=mmenu)
            for marker in ['circle', 'triangle', 'square', 'plus', 'star', 'diamond', 'x']:
                mmenu.add_command(label=marker, command=lambda state=state, marker=marker: self.view.cm.update('markers', (state,marker)))
            lmenu= Tk.Menu(menu, tearoff=0)
            menu.add_cascade(label='linestyles', menu=lmenu)
            for linestyle in ['solid', 'dashed', 'dotted', 'dashdot']:
                lmenu.add_command(label=linestyle,
                                  command=lambda state=state, linestyle=linestyle: self.view.cm.update('linestyles',
                                                                                                 (state, linestyle)))
            if state in self.main.set_colors.keys():
                state_color = self.main.set_colors[state]
            else:
                state_color = next(color_cycle)
                self.main.set_colors[state] = state_color
            if self.settings["markers"] != {}:
                if state in self.settings["markers"].keys():
                    state_marker = self.markers[str(self.settings["markers"][state])]
                else:
                    state_marker = 'o'
            else:
                state_marker = 'o'
            if self.settings["linestyles"] != {}:
                if state in self.settings["linestyles"].keys():
                    state_linestyle = str(self.settings["linestyles"][state])
                else:
                    state_linestyle = 'solid'
            else:
                state_linestyle = 'solid'
            pep_exposures = list(set([round(float(i['Exposure']),2) for i in pep[state]]))
            if len(pep_exposures) >= 3:
                x_values = []
                y_values = []
                y_error = []
                x0 = np.array([6, 4, 1])

                # Data Points
                fit_exposures =sorted(pep_exposures)
                for x in fit_exposures:
                    peptides = pep[state]
                    if self.main.corrected == 'Yes':
                        for i in peptides:
                            x_values.append(float(i['Exposure']))
                            if self.settings["yaxis_units"] == 'Da':
                                y_values.append(float(i['Uptake_corr']))
                                y_error.append(float(i['Uptake_SD_corr']))
                            elif self.settings["yaxis_units"] == 'percentDa':
                                y_values.append(float(i['Uptake_corr'])/max_uptake)
                                y_error.append(float(i['Uptake_SD_corr'])/max_uptake)
                            else:
                                return
                    else:
                        for i in peptides:
                            x_values.append(float(i['Exposure']))
                            if self.settings["yaxis_units"] == 'Da':
                                y_values.append(float(i['Uptake']))
                                y_error.append(float(i['Uptake SD']))
                            elif self.settings["yaxis_units"] == 'percentDa':
                                y_values.append(float(i['Uptake'])/max_uptake)
                                y_error.append(float(i['Uptake SD'])/max_uptake)
                            else:
                                return
                x_np = np.array(x_values)
                y_np = np.array(y_values)
                if any([i<0 for i in y_np]):
                    y_np = np.zeros(len(y_np))
                y_error_np = np.array(y_error)

                # Line Calculation

                res_lsq = leastsq(lin_fit, x0, args=(x_np, y_np))
                lowest_non_zero = sorted(list(set(x_values)))[1]
                lowest_non_zero_index = x_values.index(lowest_non_zero)
                lowest_non_zero_y = y_values[lowest_non_zero_index]
                new_y = gen_lin_fit(res_lsq[0], float(lowest_non_zero / 4))
                if new_y > float(lowest_non_zero_y / 2):
                    res_lsq2 = leastsq(lin_fit_x1, 4, args=(
                    lowest_non_zero / 4, (y_np[x_values.index(lowest_non_zero)] * 2 / 3), res_lsq[0]))
                    res_lsq[0][1] = res_lsq2[0][0]
                    res_lsq = leastsq(lin_fit_x0x2, res_lsq[0], args=(x_np, y_np, res_lsq[0][1]))
                if self.settings['fit'] == 'log':
                    t_space = np.linspace(min(sorted(x_values)[1:]), max(x_values), len(x_values) * 50)
                    y_lsq = gen_lin_fit(res_lsq[0], t_space)
                    x_np = np.array(x_values[1:])
                    y_np = np.array(y_values[1:])
                    y_error_np = np.array(y_error[1:])
                else:
                    t_space = np.linspace(min(x_values), max(x_values), len(x_values) * 50)
                    y_lsq = gen_lin_fit(res_lsq[0], t_space)

                # # First exponential fitting solution
                # res_lsq2 = leastsq(fun2e_solvex1, 4, args = (0.1, 1, res_lsq[0]))
                # maxrate = res_lsq2[0][0]
                # if res_lsq[0][1] > maxrate:
                #     res_lsq[0][1] = maxrate
                #     res_lsq = leastsq(fun2e_solvex0x2, res_lsq[0], args=(x_np, y_np, maxrate))

                # # First exponential fitting solution
                # x_index = [index for index in range(0, len(x_np)) if x_np[index] > 0][0]  # Added to prevent recomb from zeroing out the plot slope 20180909
                # if res_lsq[0][1] >= (3 * y_np[x_index]):
                #    res_lsq[0][1] = (3 * y_np[x_index])

                # Figure Plotting
                err = axis.errorbar(x_np, y_np, yerr=y_error_np, fmt=state_marker,
                                    elinewidth=float(self.settings['ls']) / 2, capthick=2,
                                    capsize=5, ms=self.settings['ms'], color=state_color)
                for i in zip(x_np,y_np):
                    names[i] = state
                errs.append(err)
                line, = axis.plot(t_space, y_lsq, label = state, linestyle=state_linestyle, color=state_color,
                                  linewidth=self.settings['ls'], marker=state_marker, markevery=[0])
                lines.append(line)
                print("State=" + state + " y=" + str(round(float(res_lsq[0][0]),2)) + "*(1-e^(-" + str(round(float(res_lsq[0][1]),2)) + "*t))+" + str(
                    round(float(res_lsq[0][2]),2)) + "*(1-e^-0.01*t)")
                # print("State=" + state + " y=" + str(round(float(res_lsq[0][0]),3)) + "*t-" + str(round(float(res_lsq[0][1]),3)))
            else:
                print("Cannot Plot with <=2 points")

        if self.settings['fit'] == 'log':
            axis.set_xscale('log')
            axis.set_xbound((sorted(x_np)[0]-0.1, max(x_np) + .1))
            axis.set_xlim((sorted(x_np)[0]-0.1, max(x_np) + .1))
        else:
            axis.set_xscale('linear')
            axis.set_xbound((0, max(x_np) + .1))
            axis.set_xlim((0, max(x_np) + .1))

        if self.settings["yaxis_units"] == 'Da':
            axis.set_ybound((0, max_uptake))
            axis.set_ylim((0, max_uptake+.1))
        elif self.settings["yaxis_units"] == 'percentDa':
            axis.set_ybound((0, 1))
            axis.set_ylim((0, 1.1))


        # Make the Legend
        legline_dict = {}
        if self.settings['leg']:
            if pep['options']['leg_bbox'] != '':
                leg = axis.legend(loc=pep['options']['leg_bbox'], prop=fontsetting, fancybox=True,
                                    shadow=True, numpoints = 1)
                if not export:
                    leg.set_draggable(state=True)
            else:
                leg = axis.legend(loc=0, prop=fontsetting, fancybox=True,
                                    shadow=True, numpoints = 1)
                if not export:
                    leg.set_draggable(state=True)

            for legline, origline, err in zip(leg.get_lines(), lines, errs):
                legline.set_picker(10)
                legline_dict[legline] = (origline, err)
            leg.set_picker(10)
        else:
            leg = axis.legend(numpoints=1)
            leg.set_picker(10)
            for legline, origline, err in zip(leg.get_lines(), lines, errs):
                legline.set_picker(10)
                legline_dict[legline] = (origline, err)
            leg.remove()

        # Set Axis Bounds

        axis.spines['bottom'].set_linewidth(self.settings['as'])
        axis.spines['left'].set_linewidth(self.settings['as'])


        # Place Labels
        for tick in axis.get_xticklabels():
            tick.set_fontname(self.settings['font'])
            tick.set_fontsize(self.settings['xaxis_size'])
        for tick in axis.get_yticklabels():
            tick.set_fontname(self.settings['font'])
            tick.set_fontsize(self.settings['yaxis_size'])
        if self.settings['xaxis']:
            axis.set_xlabel("\n".join(textwrap.wrap(self.settings["xaxis_title"],40)), fontname=self.settings['font'],
                            fontsize=self.settings['xaxis_size'])
        if self.settings['yaxis']:
            if self.settings["yaxis_title"] not in ['Deuterium Uptake (Da)','Deuterium Uptake (%)']:
                axis.set_ylabel("\n".join(textwrap.wrap(self.settings["yaxis_title"], 40)),
                                fontname=self.settings['font'],
                                fontsize=self.settings['yaxis_size'])
            elif self.settings["yaxis_units"] == 'Da':
                axis.set_ylabel("\n".join(textwrap.wrap('Deuterium Uptake (Da)', 40)),
                                fontname=self.settings['font'],
                                fontsize=self.settings['yaxis_size'])
            elif self.settings["yaxis_units"] == 'percentDa':
                axis.set_ylabel("\n".join(textwrap.wrap('Deuterium Uptake (%)', 40)),
                                fontname=self.settings['font'],
                                fontsize=self.settings['yaxis_size'])
            else:
                return
        axis.minorticks_off()

        # Title
        if self.settings['title']:
            axis.set_title("\n".join(textwrap.wrap(protein + ' ' + str(default_pep['Start']) + '-' + str(default_pep['End']),30)), fontname=self.settings['font'],
                             fontsize=self.settings['title_size'])

        # X Label
        if self.settings['xaxis_labels'] == 'Fitted':
            axis.set_xticks(x_np)
            axis.set_xticklabels([round(i, 2) for i in x_np])
        else:
            xlabels = self.settings['xaxis_labels']
            xlabels = xlabels.split(',')
            print(xlabels)
            axis.set_xticks([round(float(i), 2) for i in xlabels])
            axis.set_xticklabels([round(float(i), 2) for i in xlabels])

        if float(default_pep['MaxUptake']) <= 10:
            step = 1
        elif 10 < float(default_pep['MaxUptake']) <= 15:
            if float(default_pep['MaxUptake']) % 2 == 0:
                step = 2
            else:
                step = int(float(default_pep['MaxUptake']) / 5)
        elif 15 < float(default_pep['MaxUptake']) <= 20:
            if float(default_pep['MaxUptake']) % 3 == 0:
                step = 3
            elif float(default_pep['MaxUptake']) % 2 == 0:
                step = 2
            else:
                step = int(float(default_pep['MaxUptake']) / 5)
        elif 20 < float(default_pep['MaxUptake']) <= 25:
            if float(default_pep['MaxUptake']) % 4 == 0:
                step = 4
            elif float(default_pep['MaxUptake']) % 3 == 0:
                step = 3
            else:
                step = int(float(default_pep['MaxUptake']) / 5)
        elif 25 < float(default_pep['MaxUptake']) <= 30:
            if float(default_pep['MaxUptake']) % 5 == 0:
                step = 5
            elif float(default_pep['MaxUptake']) % 4 == 0:
                step = 4
            else:
                step = int(float(default_pep['MaxUptake']) / 5)
        else:
            if float(default_pep['MaxUptake']) % 7 == 0:
                step = 7
            elif float(default_pep['MaxUptake']) % 6 == 0:
                step = 6
            elif float(default_pep['MaxUptake']) % 5 == 0:
                step = 5
            else:
                step = int(float(default_pep['MaxUptake']) / 5)
        yaxis = list(range(0, int(float(default_pep['MaxUptake'])), step))
        yaxis.append(int(float(default_pep['MaxUptake'])))
        if self.settings["yaxis_units"] == 'Da':
            axis.set_yticks(list(set(yaxis)))
        elif self.settings["yaxis_units"] == 'percentDa':
            axis.set_yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
        # axis.set_yticklabels()
        chunks = [pep_sequence[i:i + 15] for i in range(0, len(pep_sequence), 15)]
        ann_text = str(pep_mass) + "\n" + protein_name + " " + protein_range
        if default_pep['Fragment'] == 'C':
            chunks[0] = '-'+chunks[0]
        for i in chunks:
            ann_text = ann_text + "\n" + i
        if self.settings['ann']:
            bbox_props = dict(boxstyle="square,pad=0.3", alpha=0, ec="b", lw=0)
            if pep['options']['ann_bbox'] != '':
                ann = axis.annotate(ann_text, xy = pep['options']['ann_bbox'][0], fontname=self.settings['font'],
                                      xycoords='axes fraction',
                                      fontsize=self.settings['ann_size'], ha=pep['options']['ann_bbox'][1], va=pep['options']['ann_bbox'][2],
                                      bbox = bbox_props, wrap=True)
                if not export:
                    ann.draggable()
            else:
                if max(y_np) < (max_uptake * 0.5):
                    ha_val = 'left'
                    va_val = 'top'
                    ann = axis.annotate(ann_text, xy=(0.01, 0.99), fontname=self.settings['font'],
                                             xycoords='axes fraction',
                                          fontsize=self.settings['ann_size'], ha=ha_val, va=va_val, bbox = bbox_props, wrap=True)
                    if not export:
                        ann.draggable()
                else:
                    ha_val = 'right'
                    va_val = 'bottom'
                    ann = axis.annotate(ann_text, xy=(0.99, 0.01), fontname=self.settings['font'],
                                             xycoords='axes fraction',
                                          fontsize=self.settings['ann_size'],
                                          ha=ha_val, va=va_val, bbox = bbox_props, wrap=True)
                    if not export:
                        ann.draggable()
        else:
            return

        # Create Slider
        # class FakeSlider(object):
        #     def __init__(self, val):
        #         self.val = val
        # def update_x(val):
        #     pos = spos.val
        #     axis.axis([0, pos, 0, max_uptake])
        #     self.pv.plot.draw_idle()
        # if any([0.0<(abs(e-other_e)/max(pep_exposures))<0.1 for other_e in pep_exposures for e in pep_exposures]):
        #     #self.pv.figure.subplots_adjust(bottom=0.05)
        #     axis.axis([0,max(pep_exposures),0,max_uptake])
        #     #cax = self.pv.figure.add_subplot(212)
        #     #self.pv.figure.subplots_adjust(hspace=0.1)
        #     divider = make_axes_locatable(axis)
        #     cax = divider.new_vertical(size="3%", pad=1, pack_start=True)
        #     self.pv.figure.add_axes(cax)
        #     min_val = float(min(list(sorted(pep_exposures))[1:]))
        #     spos = Slider(cax, 'Scroll X Axis', 0, max(pep_exposures), valinit=max(pep_exposures), slidermin=FakeSlider(min_val))
        #     spos.on_changed(update_x)
        #


        # Hovering Annotation
        annot = axis.annotate("", xy=(0, 0), xytext=(0, 20), textcoords="offset points",fontsize=10,fontname='Arial',
                            bbox=dict(boxstyle="round", fc="w"),
                            arrowprops=dict(arrowstyle="->"))
        annot.set_visible(False)
        self.view.figcan.mpl_connect("motion_notify_event", hover)
        #


        # Annotation Selection
        self.view.figcan.mpl_connect('pick_event', onpick)
        self.view.figcan.mpl_connect('button_release_event', lambda e: on_release(e,ann,leg))
        #

        # Resize Management
        self.view.figcan.mpl_connect('resize_event', on_resize)
        self.view.top.bind('<ButtonRelease>', done_resizing)
        self.view.top.bind('<Enter>', done_resizing)
        self.view.top.bind('<Leave>', done_resizing)

    # Save Dialog
    def save(self):
        def single():
            filetype = str(self.psv.file_type_str.get())
            print(filetype)
            print("lower "+filetype.lower())
            dpival = float(self.psv.combobox_dpi.get())
            default_pep = self.selected_peptides[0][self.states[0]][0]
            if str(default_pep['Fragment']) == "" and str(default_pep['Modification'])=="":
                initpath = self.protein + "_" + str(default_pep['Start']) + "-" + str(default_pep['End']) + "." + \
                           filetype.lower()
            elif str(default_pep['Fragment']) == "" and str(default_pep['Modification']) != "":
                initpath = self.protein + "_" + str(default_pep['Start']) + "-" + str(default_pep['End']) + "_" + \
                           str(default_pep['Modification']) + "." + filetype.lower()
            elif str(default_pep['Fragment']) != "" and str(default_pep['Modification']) == "":
                initpath = self.protein + "_" + str(default_pep['Start']) + "-" + str(default_pep['End']) + "_" + \
                           str(default_pep['Fragment']) + "." + filetype.lower()
            else:
                initpath = self.protein + "_" + str(default_pep['Start']) + "-" + str(default_pep['End']) + "_" + \
                           str(default_pep['Fragment']) + "_" + str(default_pep['Modification']) + "." + filetype.lower()
            print("Path "+initpath)
            savefile = FileDialog.asksaveasfilename(
                title="Choose save location",
                initialfile=initpath,
                initialdir=expanduser('~'))
            self.view.figure.savefig(savefile, dpi=dpival, transparent=True)
            self.psv.top.destroy()
        def multi():
            # Prompting save directory
            savedir = FileDialog.askdirectory(title="Choose save location", initialdir=expanduser('~'))
            filetype = str(self.psv.file_type_str.get())
            dpival = float(self.psv.combobox_dpi.get())
            n = 0
            for i in self.selected_peptides:
                self.view.figure.clear()
                #self.view.figure.set_size_inches(8 / self.settings['ar'], 8)
                self.view.figure.set_dpi(dpival)
                axis = self.view.figure.add_subplot(111)
                self.uptake(self.protein, i, self.state_dict, axis, export=1)
                self.view.figcan.draw()
                default_pep = i[self.states[0]][0]
                filename = self.protein + "_" + str(default_pep['Start']) + "-" + str(default_pep['End']) + "_" + \
                           str(default_pep['Fragment']) + "_" + str(
                    default_pep['Modification']) + "." + filetype.lower()
                self.view.figure.savefig(savedir + "/" + filename, dpi=dpival, transparent=True)
                n = n + 1
            self.psv.top.destroy()
            self.view.figure.set_dpi(72)
            self.view.figcan.draw()
        if len(self.selected_peptides) > 1:
            print('Save Multi')
            self.psv = self.view.save(self.view.top, multi=True)
            self.psv.button_save.configure(command=multi)
        else:
            self.psv = self.view.save(self.view.top)
            self.psv.button_save.configure(command = single)

    # Plot options window
    def options(self):
        def savedefaults(new_options):
            self.settings.update(new_options)
            self.select_item_data()
            pov.top.destroy()
        pov = self.view.options(self.view.top)

        if self.settings['leg']:
            pov.check_legend.select()
        else:
            pov.check_legend.deselect()

        if self.settings['ann']:
            pov.check_info.select()
        else:
            pov.check_info.deselect()

        if self.settings['title']:
            pov.check_title.select()
        else:
            pov.check_title.deselect()

        for i in sorted(self.main.set_colors.keys()):
            pov.tree_colors.insert('', 'end', values=(str(i), str(self.main.set_colors[i])))
        for i in self.states:
            if i in self.state_markers.keys():
                pov.tree_markers.insert('', 'end', values=(str(i), self.state_markers[i]))
            else:
                pov.tree_markers.insert('', 'end', values=(str(i), 'circle'))
        for i in self.states:
            if i in self.state_linestyles.keys():
                pov.tree_linestyles.insert('', 'end', values=(str(i), self.state_linestyles[i]))
            else:
                pov.tree_linestyles.insert('', 'end', values=(str(i), 'solid'))
        pov.button_save.configure(command=lambda: savedefaults({'leg': pov.var_check_legend.get(),
                                                            'ann': pov.var_check_info.get(),
                                                            'font': str(pov.var_combobox_font.get()),
                                                            'title': pov.var_check_title.get(),
                                                            'xaxis': pov.entry_xaxis.get(),
                                                            'yaxis': pov.entry_yaxis.get(),
                                                            'colors': pov.tree_colors.treeview_dict(),
                                                            'markers': pov.tree_markers.treeview_dict(),
                                                            'linestyles': pov.tree_linestyles.treeview_dict(),
                                                            'xaxis_labels': pov.entry_xlabels.get()}))
        pov.button_cancel.configure(command=lambda: pov.top.destroy())

class Map():
    '''
    Map Generator Window Controller
    Widgets:
        Protein Dropdown
        State Dropdown
        Exposure Dropdown
        Representation Dropdown
        Colormap Dropdown
        Range Dropdown
    Determine Desired settings before plotting map in new window
    '''

    # Initialize Map Generator Window, Fill widgets, Bind Buttons, Create Tk Variables
    def __init__(self, main):
        print('initiating map')
        self.style = ttk.Style()
        self.main = main

        # Initialize Tk window
        self.view = View.Map(main)
        main.widget_dict['Map'] = self.view.top
        self.main.view.window_menu.add_command(label='Map', command=lambda: self.main.focus(self.view.top))

        # Get data from root
        self.states = main.statesvar.get()
        self.proteins = main.proteinsvar.get()
        self.exposures = ['Ave']
        self.exposures.extend(main.exposuresvar.get())

        # Initialize Tk variables
        self.representation = Tk.StringVar()
        self.protein = Tk.StringVar()
        self.state1 = Tk.StringVar()
        self.state2 = Tk.StringVar()
        self.exposure = Tk.StringVar()
        self.color = Tk.StringVar()
        self.range = Tk.StringVar()

        # Load Data and Bindings into Tk
        self.view.type_combobox.configure(values=
                                          ['Uptake', 'SD', 'Uptake Diff',
                                           'Fractional Uptake', 'Fractional SD (RSD)',
                                           'Fractional Uptake Diff'])
        self.view.type_combobox.bind("<<ComboboxSelected>>", lambda e: self.rep_selected())
        self.view.type_combobox.configure(textvariable=self.representation)
        self.view.type_combobox.set('Fractional Uptake')
        self.view.protein_combobox.configure(values=self.proteins)
        self.view.protein_combobox.configure(textvariable=self.protein)
        self.view.protein_combobox.set(self.proteins[0])
        self.view.state1_combobox.configure(values=self.states)
        self.view.state1_combobox.configure(textvariable=self.state1)
        self.view.state1_combobox.set(self.states[0])
        self.view.state2_combobox.configure(values=[])
        self.view.state2_combobox.configure(textvariable=self.state2)
        self.view.state2_combobox.set(['None'])
        self.view.exposure_combobox.configure(values=self.exposures)
        self.view.exposure_combobox.configure(textvariable=self.exposure)
        self.view.exposure_combobox.set(str(max([float(i) for i in main.exposuresvar.get()])))
        cmap1 = mcolors.LinearSegmentedColormap.from_list("mycmap", ["blue", "white", "magenta"])
        self.cmaps = {'Pink-Green': 'PiYG', 'Purple-Green': 'PRGn', 'Brown-BlueGreen': 'BrBG', 'Tan-Purple': 'PuOr',
                 'Red-Black': 'RdGy', 'Red-Blue': 'RdBu', 'Red-Yellow-Blue': 'RdYlBu', 'Red-Yellow-Green': 'RdYlGn',
                 'Spectral': 'Spectral', 'Cool-Warm': 'coolwarm', 'Blue-Red': 'bwr', 'DarkBlue-DarkRed': 'seismic',
                 'Red-Purple Rainbow': 'gist_rainbow', 'Blue-Red Rainbow': 'rainbow',
                 'DarkBlue-DarkRed Rainbow': 'jet','Blue-White-Magenta':cmap1}
        self.view.color_combobox.configure(values=list(self.cmaps.keys()))
        self.view.color_combobox.configure(textvariable=self.color)
        self.view.color_combobox.set('Blue-Red')
        self.reverse = Tk.IntVar()
        self.view.reverse_button.configure(variable=self.reverse)
        self.reverse.set(0)
        self.xs = Tk.IntVar()
        self.view.xs_button.configure(variable=self.xs)
        self.xs.set(0)
        self.view.range_combobox.bind("<<ComboboxSelected>>", lambda e: self.range_selected())
        self.view.range_combobox.configure(values=['Fit to Data', 'Full Range', 'Half Range', 'Custom'])
        self.view.range_combobox.configure(textvariable=self.range)
        self.view.range_combobox.set('Fit to Data')
        self.xlim = Tk.StringVar()
        self.view.xlim_entry.configure(textvariable=self.xlim)
        self.xlim.set('50')
        self.view.cov_button.configure(command=lambda: self.map('Coverage', self))
        self.view.heat_button.configure(command=lambda: self.map('Heat', self))
        self.view.cov_save.configure(command=lambda: self.map_save('Coverage', self))
        self.view.heat_save.configure(command=lambda: self.map_save('Heat', self))
        self.view.top.protocol("WM_DELETE_WINDOW", lambda: main.on_closing(self.view.top, 'Map'))
        print('finished initiating map')

    # Updates Widgets depending on the representation selection
    def rep_selected(self):
        selection = self.representation.get()
        if selection in ['Uptake Diff','Fractional Uptake Diff']:
            self.view.state2_combobox.configure(values=self.states)
            self.view.state2_combobox['state'] = 'readonly'
            self.view.state2_label.configure(text='''State 2''')
            self.view.state2_combobox.configure(textvariable=self.state2)
            self.view.state2_combobox.set(self.states[0])
        else:
            self.view.state2_combobox.configure(values=[])
            self.view.state2_combobox['state'] = 'disabled'
            self.view.state2_label.configure(text='''''')
            self.view.state2_combobox.set(['None'])

    def range_selected(self):
        selection = self.range.get()
        if selection == 'Custom':
            self.view.min_range['state'] = 'normal'
            self.view.max_range['state'] = 'normal'
            self.view.min_range.insert(0, 0)
            self.view.max_range.insert(0, 1)
        else:
            self.view.min_range['state'] = 'disabled'
            self.view.max_range['state'] = 'disabled'

    # Classified Map object allows for creation of multiple maps
    class map():
        # Initializes plot and right click menu
        def __init__(self, map_type, parent):
            self.parent = parent
            self.map_type = map_type
            self.data_type = str(parent.representation.get())
            self.protein = str(parent.protein.get())
            self.state = str(parent.state1.get())
            self.state2 = str(parent.state2.get())
            self.exposure = str(parent.exposure.get())
            self.color = str(parent.color.get())
            self.reverse = int(parent.reverse.get())
            self.xs = int(parent.xs.get())
            self.cmaps = parent.cmaps
            self.data_range = str(parent.range.get())
            self.line_length = int(parent.xlim.get())
            self.corrected = parent.main.corrected
            self.obj = parent.main.state_obj[parent.main.csvfile]

            self.settings = {'Font':'Arial',
                             'Show Title':0,
                             'Title Size':32}
            self.peplist = {0: []}
            self.pepdict = []
            self.peptides = {}
            self.peptides[0] = []
            self.peptides[1] = []
            self.vertices = []
            self.codes = []
            self.textsize = 16
            self.xpad = 0

            # Create Figure & Axis and apply plot properties
            self.fig = plt.figure()
            self.ax = self.fig.add_axes([0.02, 0.1, 0.94, 0.9])
            self.fig.canvas.toolbar.pack_forget()
            self.fig.set_dpi(72)
            self.fig.subplots_adjust(top=1, bottom=0.0, left=0.05, right=0.95)
            self.ax.axis("off")

            self.title_val = Tk.IntVar()
            self.title_val.set(1)

            # Creating axes, ticks, numbers, and letters
            self.make_cmap()
            self.get_peptide_data()
            self.set_color_values()
            self.determine_levels()
            self.create_axes_dnx(self.ax)
            self.add_peptides_dnx(self.ax)
            self.draw(self.ax, self.fig)
            self.scrollbars()
            self.axis_colorbar()
            self.popupMenu = Tk.Menu(self.fig.canvas.get_tk_widget(), tearoff=0)
            self.popupMenu.add_command(label='Save', command=lambda: self.export())
            self.popupMenu.add_command(label='Close', command=lambda: plt.close(self.fig))
            if platform == 'darwin':
                self.fig.canvas.get_tk_widget().bind("<Button-2>", lambda event: self.popupMenu.post(event.x_root, event.y_root))
            else:
                self.fig.canvas.get_tk_widget().bind("<Button-3>", lambda event: self.popupMenu.post(event.x_root, event.y_root))
            plt.show()

        # Saves plot to file
        def export(self):
            fig = plt.figure()
            ax = fig.add_axes([0.01, 0.01, 0.98, 0.98])
            fig.canvas.toolbar.pack_forget()
            fig.set_dpi(300)
            ax.axis("off")
            self.create_axes_dnx(ax)
            self.add_peptides_dnx(ax)
            self.inline_colorbar(ax)
            self.draw(ax, fig)
            plt.draw()
            if self.data_type in ['Uptake Diff', 'Fractional Uptake Diff']:
                savefile = FileDialog.asksaveasfilename(
                    title="Choose save location",
                    initialfile=self.protein + "_" + str(self.state) + "-" + str(self.state2) + "_" + str(self.exposure)+"_Map.png",
                    initialdir=self.parent.main.dir,
                    filetypes=(("PNG", '*.png'),("JPG",'*.jpg'),("TIF", '*.tif'), ("all files", "*.*")))
            else:
                savefile = FileDialog.asksaveasfilename(
                    title="Choose save location",
                    initialfile=self.protein + "_" + str(self.state) + "_" + str(self.exposure) + "_Map.png",
                    initialdir=self.parent.main.dir,
                    filetypes=(("PNG", '*.png'),("JPG",'*.jpg'),("TIF", '*.tif'), ("all files", "*.*")))
            if split(savefile)[0] != "":
                self.parent.main.dir = split(savefile)[0]
            if str(splitext(split(savefile)[1])[1]).lower() not in ['.png','.jpg','.tif']:
                savefile = savefile + '.png'
            fig.savefig(savefile, dpi=300, transparent=True)
            plt.close()

        # Customizes colormap
        def make_cmap(self):
            num_colors = 102
            cm = plt.get_cmap(self.cmaps[self.color])
            colorlist = [cm(1. * i / num_colors) for i in range(num_colors)]
            if not self.reverse:
                colorlist = list(colorlist)
            else:
                colorlist = list(reversed(colorlist))
            self.rearranged_colors = colorlist[1:101]
            self.newmap = matplotlib.colors.LinearSegmentedColormap.from_list("newmap", self.rearranged_colors)
            self.color_values = ()
            self.cval = 0
            self.mval = 0

        # Sets the color range based on specified settings and the data range
        def set_color_values(self):
            if self.data_type == 'Uptake':
                values = [i['reluptake'] for i in self.pepdict]
                if self.data_range in ['Fit to Data', 'Full Range', 'Half Range']:
                    color_values = (min(values), max(values))
                else:
                    self.custom_min = float(self.parent.view.min_range.get())
                    self.custom_max = float(self.parent.view.max_range.get())
                    color_values = (self.custom_min, self.custom_max)
            elif self.data_type == 'SD':
                values = [i['relsd'] for i in self.pepdict]
                if self.data_range in ['Fit to Data', 'Full Range', 'Half Range']:
                    color_values = (min(values), max(values))
                else:
                    self.custom_min = float(self.parent.view.min_range.get())
                    self.custom_max = float(self.parent.view.max_range.get())
                    color_values = (self.custom_min, self.custom_max)
            elif self.data_type == 'Uptake Diff':
                values = [i['uptakediff'] for i in self.pepdict]
                if self.data_range in ['Fit to Data', 'Full Range', 'Half Range']:
                    color_values = (min(values), max(values))
                else:
                    self.custom_min = float(self.parent.view.min_range.get())
                    self.custom_max = float(self.parent.view.max_range.get())
                    color_values = (self.custom_min, self.custom_max)
            elif self.data_type == 'Fractional Uptake':
                values = [i['fracuptake'] for i in self.pepdict]
                if self.data_range in ['Fit to Data']:
                    color_values = (min(values), max(values))
                elif self.data_range in ['Full Range']:
                    color_values = (0, 1)
                elif self.data_range in ['Half Range']:
                    color_values = (0, 0.5)
                else:
                    self.custom_min = float(self.parent.view.min_range.get())
                    self.custom_max = float(self.parent.view.max_range.get())
                    color_values = (self.custom_min, self.custom_max)
            elif self.data_type == 'Fractional SD (RSD)':
                values = [i['fracsd'] for i in self.pepdict]
                if self.data_range in ['Fit to Data']:
                    color_values = (min(values), max(values))
                elif self.data_range in ['Full Range']:
                    color_values = (0, 1)
                elif self.data_range in ['Half Range']:
                    color_values = (0, 0.5)
                else:
                    self.custom_min = float(self.parent.view.min_range.get())
                    self.custom_max = float(self.parent.view.max_range.get())
                    color_values = (self.custom_min, self.custom_max)
            elif self.data_type == 'Fractional Uptake Diff':
                values = [i['fracuptakediff'] for i in self.pepdict]
                if self.data_range in ['Fit to Data']:
                    color_values = (min(values), max(values))
                elif self.data_range in ['Full Range']:
                    color_values = (-1, 1)
                elif self.data_range in ['Half Range']:
                    color_values = (-0.5, 0.5)
                else:
                    self.custom_min = float(self.parent.view.min_range.get())
                    self.custom_max = float(self.parent.view.max_range.get())
                    color_values = (self.custom_min, self.custom_max)
            else:
                print("Type Error")
            if color_values[0] == 0 and color_values[1] == 0:
                cval = 0
                mval = 0
            elif color_values[0] == 0 and color_values[1] != 0:
                cval = 0
                mval = 100 / color_values[1]
            elif color_values[1] == color_values[0]:
                print("Range Error")
                cval = 0
                mval = 0
            else:
                cval = float(100 / (-(color_values[1] / color_values[0]) + 1))
                mval = float((100 - cval) / color_values[1])
            self.color_values = color_values
            self.cval = cval
            self.mval = mval

        # Obtains peptide uptake information
        def get_peptide_data(self):
            # Doesn't work for AVE exposure
            if self.corrected == 'Yes':
                self.uptakekey = 'Uptake_corr'
                self.sdkey = 'Uptake_SD_corr'
            else:
                self.uptakekey = 'Uptake'
                self.sdkey = 'Uptake SD'
            # Get peptide list and exposure list
            if self.data_type in ['Uptake Diff', 'Fractional Uptake Diff']:
                if self.exposure == 'Ave':
                    exposures1 = [float(i) for i in self.obj.data_nest[self.protein][self.state].keys()]
                    exposures2 = [float(i) for i in self.obj.data_nest[self.protein][self.state2].keys()]
                    exposures = list(set(exposures1).intersection(exposures2))
                    exposure = min(exposures)
                else:
                    exposure = float(self.exposure)
                dict0 = sorted([i['PepID'] for i in self.obj.data_nest[self.protein][self.state][exposure] if i['Fragment'] == ''])
                dict1 = sorted([i['PepID'] for i in self.obj.data_nest[self.protein][self.state2][exposure] if i['Fragment'] == ''])
                id_list = list(set(dict0).intersection(dict1))
            else:
                if self.exposure == 'Ave':
                    exposures = [float(i) for i in self.obj.data_nest[self.protein][self.state].keys()]
                    exposure = exposures[0]
                    id_list = sorted([i['PepID'] for i in self.obj.data_nest[self.protein][self.state][exposure] if i['Fragment'] == ''])
                else:
                    exposure = float(self.exposure)
                    id_list = sorted([i['PepID'] for i in self.obj.data_nest[self.protein][self.state][exposure] if i['Fragment'] == ''])

            # Get data for all peptides
            if self.exposure == 'Ave':
                if 0 in exposures:
                    exposures.remove(0)
                for pepid in id_list:
                    elist = []
                    for e in sorted(exposures):
                        elist.extend([n for n in self.obj.data_nest[self.protein][self.state][e] if n["PepID"]==pepid])
                    reluptake = np.mean([float(i[self.uptakekey]) for i in elist])
                    relsd = np.mean([float(i[self.sdkey]) for i in elist])
                    if elist[0]['MaxUptake'] == 0:
                        fracuptake = 0
                    else:
                        fracuptake = np.mean([float(i[self.uptakekey])/float(i['MaxUptake']) for i in elist])
                    if float(elist[0][self.uptakekey]) == 0:
                        fracsd = 0
                    else:
                        fracsd = np.mean([float(i[self.sdkey])/float(i[self.uptakekey]) for i in elist])
                    if self.data_type in ['Uptake Diff','Fractional Uptake Diff']:
                        e2list = []
                        for e in sorted(exposures):
                            e2list.extend([n for n in self.obj.data_nest[self.protein][self.state2][e] if
                                 n["PepID"]==pepid])
                        uptakediff = np.mean([float(j[self.uptakekey]) - float(i[self.uptakekey]) for j in e2list for i in elist])
                        if e2list[0]['MaxUptake'] == 0:
                            fracuptakediff = 0
                        else:
                            fracuptakediff = np.mean([float(j[self.uptakekey])/float(j['MaxUptake']) -
                                                              float(i[self.uptakekey])/float(i['MaxUptake'])
                                                              for j in e2list for i in elist])
                        self.pepdict.append({'Start':elist[0]['Start'],'End':elist[0]['End'], 'Fragment':elist[0]['Fragment'],  'Modification':elist[0]['Modification'], 'MaxUptake':elist[0]['MaxUptake'],
                                            'reluptake':reluptake, 'uptakediff':uptakediff, 'fracuptake':fracuptake,
                                            'fracuptakediff':fracuptakediff, 'relsd':relsd, 'fracsd':fracsd})
                    else:
                        self.pepdict.append({'Start':elist[0]['Start'],'End':elist[0]['End'], 'Fragment':elist[0]['Fragment'], 'Modification':elist[0]['Modification'], 'MaxUptake':elist[0]['MaxUptake'],
                                            'reluptake':reluptake, 'fracuptake':fracuptake, 'relsd':relsd, 'fracsd':fracsd})
            else:
                for pepid in id_list:
                    pep = [n for n in self.obj.data_nest[self.protein][self.state][exposure] if n["PepID"]==pepid][0]
                    reluptake = float(pep[self.uptakekey])
                    relsd = float(pep[self.sdkey])
                    if pep['MaxUptake'] == 0:
                        fracuptake = 0
                    else:
                        fracuptake = float(pep[self.uptakekey]) / float(pep['MaxUptake'])
                    if float(pep[self.uptakekey]) == 0:
                        fracsd = 0
                    else:
                        fracsd = float(pep[self.sdkey]) / float(pep[self.uptakekey])
                    if self.data_type in ['Uptake Diff','Fractional Uptake Diff']:
                        pep2 = [n for n in self.obj.data_nest[self.protein][self.state2][exposure] if n["PepID"]==pepid][0]
                        uptakediff = float(pep2[self.uptakekey]) - float(pep[self.uptakekey])
                        if pep2['MaxUptake'] == 0:
                            fracuptakediff = 0
                        else:
                            fracuptakediff = (float(pep2[self.uptakekey]) - float(pep[self.uptakekey])) / float(pep['MaxUptake'])
                        self.pepdict.append({'Start':pep['Start'],'End':pep['End'], 'Fragment':pep['Fragment'], 'Modification':pep['Modification'],'MaxUptake':pep['MaxUptake'],
                                            'reluptake':reluptake, 'uptakediff':uptakediff, 'fracuptake':fracuptake,
                                            'fracuptakediff':fracuptakediff, 'relsd':relsd, 'fracsd':fracsd})
                    else:
                        self.pepdict.append({'Start':pep['Start'],'End':pep['End'], 'Fragment':pep['Fragment'], 'Modification':pep['Modification'], 'MaxUptake':pep['MaxUptake'],
                                            'reluptake':reluptake, 'fracuptake':fracuptake, 'relsd':relsd, 'fracsd':fracsd})
            if self.map_type == 'Heat':
                if self.data_type in ['Fractional Uptake Diff', 'Uptake Diff']:
                    if self.exposure == 'Ave':
                        exposures1 = [float(i) for i in
                                      self.obj.data_nest[self.protein][self.state].keys()]
                        exposures2 = [float(i) for i in self.obj.data_nest[self.protein][self.state2].keys()]
                        exposures = list(set(exposures1).intersection(exposures2))
                        exposure = min(exposures)
                        self.obj.assignData(self.protein, self.state, 'Ave', self.state2)
                    else:
                        exposure = float(self.exposure)
                        self.obj.assignData(self.protein, self.state, exposure, self.state2)
                else:

                    if self.exposure == 'Ave':
                        exposures = [float(i) for i in self.obj.data_nest[self.protein][self.state].keys()]
                        exposure = exposures[0]
                        self.obj.assignData(self.protein, self.state, 'Ave')
                    else:
                        exposure = float(self.exposure)
                        self.obj.assignData(self.protein, self.state, exposure)
                self.seq_nest2 = {}
                for item in self.obj.seq_nest:
                    if item['assigned'] == 1:
                        if (item['pepmin'], item['pepmax']) in self.seq_nest2.keys():
                            self.seq_nest2[(item['pepmin'], item['pepmax'])].append(item['pos'])
                        else:
                            self.seq_nest2[(item['pepmin'], item['pepmax'])] = [item['pos']]

            self.sequence = self.obj.sequences[self.protein]
            if not self.xs:
                self.firstpos = min([i for i,j in enumerate(self.sequence) if j != "X"])
                self.startspace = 0
            else:
                self.firstpos = min([i for i,j in enumerate(self.sequence) if j != "X"])
                self.startspace = (self.firstpos)%self.line_length
            self.sequence = self.sequence[self.firstpos:]

        # Determines the height of the map depending on the line width, sequence length, and peptide coverage
        def determine_levels(self):
            self.seqlist = list(self.sequence)
            self.seqlength = len(self.seqlist)+self.startspace
            self.levels = int(self.seqlength / self.line_length) + 1
            self.level_dict = {}
            for pep in sorted(self.pepdict, key=itemgetter('Start')):
                n = 0
                while set(self.peplist[n]).intersection(list(range(pep['Start'], pep['End'] + 1))) != set([]):
                    n = n + 1
                    if n not in self.peplist.keys():
                        self.peplist[n] = []
                self.peplist[n] += range(pep['Start'], pep['End'] + 1)
                # Added modification to peptide descriptor
                self.level_dict[(pep['Start'], pep['End'], pep['Fragment'],pep['Modification'])] = n

            if self.map_type == 'Coverage':
                self.redun_dict = {}
                self.spacing = 4
                self.ymax=2
                for level in range(0,self.levels):
                    self.redun_dict[level+1] = 0
                    range_min = (level * self.line_length) + 1
                    range_max = (level + 1) * self.line_length
                    range_list = list(range(range_min, range_max + 1))
                    for n in self.peplist.keys():
                        if set([resn-self.firstpos+self.startspace for resn in self.peplist[n] if range_min<=(resn-self.firstpos+self.startspace)<=range_max]).intersection(range_list) != set([]):
                            self.redun_dict[level+1] = n+1
                    self.ymax += (self.redun_dict[level+1] + self.spacing)
            elif self.map_type == 'Heat':
                self.redun_dict = {}
                for level in range(0,self.levels):
                    self.redun_dict[level+1] = 1
                self.spacing = 4
                self.ymax = (self.levels * (1 + self.spacing)) + 2
            else:
                print('type error')

        # Initializes axis in inverted DynamX Style
        def create_axes(self, axis):
            self.ycurlevel = self.ymax
            for n in range(0, self.levels):
                level = n + 1
                redundancy = self.redun_dict[level]
                self.ycurlevel = self.ycurlevel - (redundancy + self.spacing)
                # If first line
                if level == 1:
                    # Create a line from xpad to linelength+xpad
                    self.vertices += [(2 * (self.startspace+self.xpad), 2 * self.ycurlevel),
                                      (2 * (self.line_length + self.xpad), 2 * self.ycurlevel)]
                    self.codes += [Path.MOVETO, Path.LINETO]
                    # Create Hash marks and text
                    for m in range(self.startspace, (n + 1) * self.line_length):
                        x = m % self.line_length
                        self.vertices += [(2 * (x + self.xpad), 2 * (self.ycurlevel + 0.1)),
                                          (2 * (x + self.xpad), 2 * (self.ycurlevel - 0.1))]
                        self.codes += [Path.MOVETO, Path.LINETO]
                        if (m + 1 + self.firstpos - self.startspace) % 10 == 0:
                            # 2020-06-04 m + 1 changed to m + firstpos + 1 to start map at first coverage
                            axis.annotate(str(m - self.startspace + self.firstpos + 1),
                                          (2 * (x + self.xpad + 0.5), 2 * (self.ycurlevel + 0.5)),
                                          color='black', fontsize=self.textsize, horizontalalignment='center',
                                          verticalalignment='bottom', family='monospace')

                    for m in range(self.startspace, (n + 1) * self.line_length):
                        x = m % self.line_length
                        axis.annotate(self.seqlist[m - self.startspace], (2 * (x + self.xpad + 0.5), 2 * (self.ycurlevel - 0.5)),
                                      color='black', fontsize=self.textsize, horizontalalignment='center',
                                      verticalalignment='top',
                                      family='monospace')
                    self.vertices += [(2 * (self.line_length + self.xpad), 2 * (self.ycurlevel + 0.1)),
                                      (2 * (self.line_length + self.xpad), 2 * (self.ycurlevel - 0.1))]
                    self.codes += [Path.MOVETO, Path.LINETO]
                # If sequence is longer than the current level:
                elif (self.seqlength - self.line_length * n) > self.line_length:
                    # Create a line from xpad to linelength+xpad
                    self.vertices += [(2 * self.xpad, 2 * self.ycurlevel),
                                      (2 * (self.line_length + self.xpad), 2 * self.ycurlevel)]
                    self.codes += [Path.MOVETO, Path.LINETO]
                    # Create Hash marks and text
                    for m in range(n * self.line_length, (n + 1) * self.line_length):
                        x = m % self.line_length
                        self.vertices += [(2 * (x + self.xpad), 2 * (self.ycurlevel + 0.1)),
                                          (2 * (x + self.xpad), 2 * (self.ycurlevel - 0.1))]
                        self.codes += [Path.MOVETO, Path.LINETO]
                        if (m + 1 + self.firstpos-self.startspace) % 10 == 0:
                            # 2020-06-04 m + 1 changed to m + firstpos + 1 to start map at first coverage
                            axis.annotate(str(m + self.firstpos-self.startspace + 1), (2 * (x + self.xpad + 0.5), 2 * (self.ycurlevel + 0.5)),
                                          color='black', fontsize=self.textsize, horizontalalignment='center',
                                          verticalalignment='bottom', family='monospace')

                    for m in range(n * self.line_length, (n + 1) * self.line_length):
                        x = m % self.line_length
                        axis.annotate(self.seqlist[m-self.startspace], (2 * (x + self.xpad + 0.5), 2 * (self.ycurlevel - 0.5)),
                                      color='black', fontsize=self.textsize, horizontalalignment='center',
                                      verticalalignment='top',
                                      family='monospace')
                    self.vertices += [(2 * (self.line_length + self.xpad), 2 * (self.ycurlevel + 0.1)),
                                      (2 * (self.line_length + self.xpad), 2 * (self.ycurlevel - 0.1))]
                    self.codes += [Path.MOVETO, Path.LINETO]
                # Conditions for final line
                else:
                    self.vertices += [(2 * self.xpad, 2 * self.ycurlevel),
                                      (2 * ((self.seqlength % self.line_length) + self.xpad), 2 * self.ycurlevel)]
                    self.codes += [Path.MOVETO, Path.LINETO]
                    for m in range(n * self.line_length, self.seqlength):
                        x = m % self.line_length
                        self.vertices += [(2 * (x + self.xpad), 2 * (self.ycurlevel + 0.1)),
                                          (2 * (x + self.xpad), 2 * (self.ycurlevel - 0.1))]
                        self.codes += [Path.MOVETO, Path.LINETO]
                        if (m + 1 + self.firstpos-self.startspace) % 10 == 0:
                            # 2020-06-04 m + 1 changed to m + firstpos + 1 to start map at first coverage
                            axis.annotate(str(m + self.firstpos-self.startspace + 1), (2 * (x + self.xpad + 0.5), 2 * (self.ycurlevel + 0.5)),
                                          color='black', fontsize=self.textsize, horizontalalignment='center',
                                          verticalalignment='bottom', family='monospace')
                    for m in range(n * self.line_length, self.seqlength):
                        x = m % self.line_length
                        axis.annotate(self.seqlist[m-self.startspace], (2 * (x + self.xpad + 0.5), 2 * (self.ycurlevel - 0.5)),
                                      color='black', fontsize=self.textsize, horizontalalignment='center',
                                      verticalalignment='top',
                                      family='monospace')
                    self.vertices += [(2 * (self.seqlength % self.line_length + 1), 2 * (self.ycurlevel + 0.1)),
                                      (2 * (self.seqlength % self.line_length + 1), 2 * (self.ycurlevel - 0.1))]
                    self.codes += [Path.MOVETO, Path.LINETO]
            axis.add_patch(
                PathPatch(Path(self.vertices, self.codes), facecolor='None', edgecolor='black',
                          linewidth=1))

        # Initializes axis in DynamX Style
        def create_axes_dnx(self, axis):
            self.ycurlevel = self.ymax - 1.5
            for n in range(0, self.levels):
                level = n + 1
                redundancy = self.redun_dict[level]
                # If sequence is longer than the current level:
                if level ==1:
                    for m in range(self.startspace, (n + 1) * self.line_length):
                        x = m % self.line_length
                        axis.annotate(self.seqlist[m-self.startspace], (2 * (x + self.xpad + 0.5), 2 * (self.ycurlevel + 0.5)),
                                      color='black', fontsize=self.textsize, horizontalalignment='center',
                                      verticalalignment='center',
                                      family='monospace')
                        if (m + 1 + self.firstpos-self.startspace) % 10 == 0:
                            # 2020-06-04 m + 1 changed to m + firstpos + 1 to start map at first coverage
                            axis.annotate(str(m + self.firstpos - self.startspace + 1), (2 * (x + self.xpad + 0.5), 2 * (self.ycurlevel - 0.5)),
                                          color='black', fontsize=self.textsize*0.75, horizontalalignment='center',
                                          verticalalignment='center', family='monospace')
                elif (self.seqlength - self.line_length * n) > self.line_length:
                    for m in range(n * self.line_length, (n + 1) * self.line_length):
                        x = m % self.line_length
                        axis.annotate(self.seqlist[m-self.startspace], (2 * (x + self.xpad + 0.5), 2 * (self.ycurlevel + 0.5)),
                                      color='black', fontsize=self.textsize, horizontalalignment='center',
                                      verticalalignment='center',
                                      family='monospace')
                        if (m + 1 + self.firstpos-self.startspace) % 10 == 0:
                            # 2020-06-04 m + 1 changed to m + firstpos + 1 to start map at first coverage
                            axis.annotate(str(m + self.firstpos - self.startspace + 1), (2 * (x + self.xpad + 0.5), 2 * (self.ycurlevel - 0.5)),
                                          color='black', fontsize=self.textsize*0.75, horizontalalignment='center',
                                          verticalalignment='center', family='monospace')
                # Conditions for final line
                else:
                    for m in range(n * self.line_length, self.seqlength):
                        x = m % self.line_length
                        axis.annotate(self.seqlist[m-self.startspace], (2 * (x + self.xpad + 0.5), 2 * (self.ycurlevel + 0.5)),
                                      color='black', fontsize=self.textsize, horizontalalignment='center',
                                      verticalalignment='center',
                                      family='monospace')
                        if (m + 1 + self.firstpos) % 10 == 0:
                            # 2020-06-04 m + 1 changed to m + firstpos + 1 to start map at first coverage
                            axis.annotate(str(m + self.firstpos - self.startspace + 1), (2 * (x + self.xpad + 0.5), 2 * (self.ycurlevel - 0.5)),
                                          color='black', fontsize=self.textsize*0.75, horizontalalignment='center',
                                          verticalalignment='center', family='monospace')
                self.ycurlevel = self.ycurlevel - (redundancy + self.spacing)

        # Add Bars to the Map in inverted DynamX Style
        def add_peptides(self, axis):
            # Creating Boxes
            if self.map_type == 'Coverage':
                for pep in self.pepdict:
                    vertices = []
                    codes = []
                    n = 0
                    peplevel = self.level_dict[(pep['Start'], pep['End'], pep['Fragment'],pep['Modification'])]
                    # 2020-06-04 pep['Start'] - 1 changed to pep['Start'] - firstpos - 1
                    # 2020-06-04 pep['End'] - 1 changed to pep['End'] - firstpos - 1
                    if int((pep['Start']-self.firstpos+self.startspace-1) / self.line_length) != int((pep['End']-self.firstpos+self.startspace-1) / self.line_length):
                        level = int((pep['Start']-self.firstpos+self.startspace-1) / self.line_length) + 1
                        yheight1 = peplevel + self.ymax - (
                                sum([self.redun_dict[j] for j in range(1, level + 1)]) + level * self.spacing) + 2
                        yheight2 = peplevel + self.ymax - (
                                sum([self.redun_dict[j] for j in range(1, level + 2)]) + (level + 1) * self.spacing) + 2
                        start50 = (pep['Start']-self.firstpos+self.startspace-1) % self.line_length + self.xpad
                        stop50 = (pep['End']-self.firstpos+self.startspace-1) % self.line_length + 1
                        codes += [Path.MOVETO] + [Path.LINETO] * 3 + [Path.MOVETO] + [Path.LINETO] * 3
                        vertices += [(2 * (self.line_length + self.xpad), 2 * yheight1), (2 * start50, 2 * yheight1),
                                     (2 * start50, 2 * (yheight1 + 0.75)),
                                     (2 * (self.line_length + self.xpad), 2 * (yheight1 + 0.75)),
                                     (2 * self.xpad, 2 * yheight2), (2 * stop50, 2 * yheight2),
                                     (2 * stop50, 2 * (yheight2 + 0.75)),
                                     (2 * self.xpad, 2 * (yheight2 + 0.75))]
                    else:
                        level = int((pep['Start'] - self.firstpos+self.startspace-1) / self.line_length) + 1
                        yheight = peplevel + self.ymax - (
                                sum([self.redun_dict[j] for j in range(1, level + 1)]) + level * self.spacing) + 2
                        start50 = (pep['Start'] - self.firstpos+self.startspace-1) % self.line_length + self.xpad
                        stop50 = (pep['End'] - self.firstpos+self.startspace-1) % self.line_length + 1
                        codes += [Path.MOVETO] + [Path.LINETO] * 3 + [Path.CLOSEPOLY]
                        vertices += [(2 * start50, 2 * yheight), (2 * stop50, 2 * yheight),
                                     (2 * stop50, 2 * (yheight + 0.75)),
                                     (2 * start50, 2 * (yheight + 0.75)),
                                     (0, 0)]
                    vertices = array(vertices, float)
                    if float(pep['MaxUptake']) == 0:
                        color_index = 0
                    else:
                        if self.data_type == 'Uptake':
                            color_index = int(self.mval * float(pep['reluptake']) + self.cval)
                        elif self.data_type == 'SD':
                            color_index = int(self.mval * float(pep['relsd']) + self.cval)
                        elif self.data_type == 'Uptake Diff':
                            color_index = int(self.mval * pep['uptakediff'] + self.cval)
                        elif self.data_type == 'Fractional Uptake':
                            color_index = int(self.mval * pep['fracuptake'] + self.cval)
                        elif self.data_type == 'Fractional SD (RSD)':
                            color_index = int(self.mval * pep['fracsd'] + self.cval)
                        elif self.data_type == 'Fractional Uptake Diff':
                            color_index = int(self.mval * pep['fracuptakediff'] + self.cval)
                        else:
                            print("Type Error")
                    if color_index < 0:
                        color_index = 0
                    elif color_index >= 100:
                        color_index = 99
                    color = self.rearranged_colors[color_index]
                    axis.add_patch(
                        PathPatch(Path(vertices, codes), facecolor=color, linewidth=1, edgecolor='grey'))
            elif self.map_type == 'Heat':
                for peptide in self.seq_nest2.keys():
                    pep = [n for n in self.pepdict if (n["Start"] == peptide[0] and n['End'] == peptide[1])][0]
                    start = min(self.seq_nest2[peptide])
                    end = max(self.seq_nest2[peptide])
                    vertices = []
                    codes = []
                    self.redundancy=1
                    # 2020-06-04 pep['Start'] - 1 changed to pep['Start'] - firstpos - 1
                    # 2020-06-04 pep['End'] - 1 changed to pep['End'] - firstpos - 1
                    if int((start - self.firstpos+self.startspace-1) / self.line_length) != int((end - self.firstpos+self.startspace-1) / self.line_length):  # If peptide starts and ends on different levels
                        level = int((start - self.firstpos+self.startspace-1) / self.line_length) + 1  # Level is vertical position at various sequence fragments
                        yheight1 = self.ymax - (level * (self.redundancy + self.spacing)) + 2  # Height results from level, redundancy(highest n value) and spacing
                        yheight2 = self.ymax - ((level + 1) * (self.redundancy + self.spacing)) + 2  # Height position at second level, hence, level+1
                        start50 = (start - self.firstpos+self.startspace-1) % self.line_length + self.xpad  # Residue 50 on a line of 1-50 should start at the end of the line, not the beginning of the next
                        stop50 = (end - self.firstpos+self.startspace-1) % self.line_length + 1  # Stop includes stop residue, hence +1
                        codes += [Path.MOVETO] + [Path.LINETO] * 3 + [Path.MOVETO] + [Path.LINETO] * 3
                        vertices += [(2*(self.line_length + self.xpad), 2*yheight1), (2*start50, 2*yheight1), (2*start50, 2*(yheight1 + 2)),
                                     (2*(self.line_length + self.xpad), 2*(yheight1 + 2)),
                                     (2*self.xpad, 2*yheight2), (2*stop50, 2*yheight2), (2*stop50, 2*(yheight2 + 2)), (2*self.xpad, 2*(yheight2 + 2))]
                    else:
                        level = int((start - self.firstpos+self.startspace-1) / self.line_length) + 1
                        yheight = self.ymax - (level * (self.redundancy + self.spacing)) + 2
                        start50 = (start - self.firstpos+self.startspace-1) % self.line_length + self.xpad
                        stop50 = (end - self.firstpos+self.startspace-1) % self.line_length + 1
                        codes += [Path.MOVETO] + [Path.LINETO] * 3 + [Path.CLOSEPOLY]
                        vertices += [(2*start50, 2*yheight), (2*stop50, 2*yheight), (2*stop50, 2*(yheight + 2)),
                                     (2*start50, 2*(yheight + 2)),
                                     (0, 0)]
                    vertices = np.array(vertices, float)
                    if float(pep['MaxUptake']) == 0:
                        color_index = 0
                    else:
                        if self.data_type == 'Uptake':
                            color_index = int(self.mval * float(pep['reluptake']) + self.cval)
                        elif self.data_type == 'SD':
                            color_index = int(self.mval * float(pep['relsd']) + self.cval)
                        elif self.data_type == 'Uptake Diff':
                            color_index = int(self.mval * pep['uptakediff'] + self.cval)
                        elif self.data_type == 'Fractional Uptake':
                            color_index = int(self.mval * pep['fracuptake'] + self.cval)
                        elif self.data_type == 'Fractional SD (RSD)':
                            color_index = int(self.mval * pep['fracsd'] + self.cval)
                        elif self.data_type == 'Fractional Uptake Diff':
                            color_index = int(self.mval * pep['fracuptakediff'] + self.cval)
                        else:
                            print("Type Error")
                        if color_index < 0:
                            color_index = 0
                        elif color_index >= 100:
                            color_index = 99
                        color = self.rearranged_colors[color_index]
                    axis.add_patch(PathPatch(Path(vertices, codes), facecolor=color, edgecolor='grey'))

        # Add Bars to the Map in DynamX Style
        def add_peptides_dnx(self, axis):
            # Creating Boxes
            if self.map_type == 'Coverage':
                # 2020-06-04 pep['Start'] - 1 changed to pep['Start'] - firstpos - 1
                # 2020-06-04 pep['End'] - 1 changed to pep['End'] - firstpos - 1
                for pep in self.pepdict:
                    vertices = []
                    codes = []
                    stretch = int((pep['End'] - self.firstpos+self.startspace-1) / self.line_length) - int((pep['Start'] - self.firstpos+self.startspace-1) / self.line_length)
                    peplevel = self.level_dict[(pep['Start'], pep['End'], pep['Fragment'], pep['Modification'])]
                    if stretch > 0:
                        level = int((pep['Start'] - self.firstpos+self.startspace-1) / self.line_length) + 1
                        yheight1 = -peplevel + self.ymax - (
                                sum([self.redun_dict[j] for j in range(1, level)]) + (level) * self.spacing) + 0.5
                        start50 = (pep['Start'] - self.firstpos+self.startspace-1) % self.line_length + self.xpad
                        codes += [Path.MOVETO] + [Path.LINETO] * 3
                        vertices += [(2 * (self.line_length + self.xpad), 2 * yheight1), (2 * start50, 2 * yheight1),
                                     (2 * start50, 2 * (yheight1 + 0.75)),
                                     (2 * (self.line_length + self.xpad), 2 * (yheight1 + 0.75))]

                        for n in range(1,stretch):
                            yheight = -peplevel + self.ymax - (sum([self.redun_dict[j] for j in range(1, level + n)]) +
                                                               (level + n) * self.spacing) + 0.5
                            codes += [Path.MOVETO] + [Path.LINETO] * 3
                            vertices += [(2 * (self.line_length + self.xpad), 2 * yheight),
                                         (2 * self.xpad, 2 * yheight),
                                         (2 * self.xpad, 2 * (yheight + 0.75)),
                                         (2 * (self.line_length + self.xpad), 2 * (yheight + 0.75))]
                        yheight2 = -peplevel + self.ymax - (
                                sum([self.redun_dict[j] for j in range(1, level+stretch)]) + (level + stretch) * self.spacing) + 0.5
                        stop50 = (pep['End'] - self.firstpos+self.startspace-1) % self.line_length + 1
                        codes += [Path.MOVETO] + [Path.LINETO] * 3
                        vertices += [(2 * self.xpad, 2 * yheight2), (2 * stop50, 2 * yheight2),
                                     (2 * stop50, 2 * (yheight2 + 0.75)),
                                     (2 * self.xpad, 2 * (yheight2 + 0.75))]
                    else:
                        level = int((pep['Start'] - self.firstpos+self.startspace-1) / self.line_length) + 1
                        yheight = -peplevel + self.ymax - (
                                sum([self.redun_dict[j] for j in range(1, level)]) + level * self.spacing) + 0.5
                        start50 = (pep['Start'] - self.firstpos+self.startspace-1) % self.line_length + self.xpad
                        stop50 = (pep['End'] - self.firstpos+self.startspace-1) % self.line_length + 1
                        codes += [Path.MOVETO] + [Path.LINETO] * 3 + [Path.CLOSEPOLY]
                        vertices += [(2 * start50, 2 * yheight), (2 * stop50, 2 * yheight),
                                     (2 * stop50, 2 * (yheight + 0.75)),
                                     (2 * start50, 2 * (yheight + 0.75)),
                                     (0, 0)]
                    vertices = array(vertices, float)
                    if float(pep['MaxUptake']) == 0:
                        color_index = 0
                    else:
                        if self.data_type == 'Uptake':
                            color_index = int(self.mval * float(pep['reluptake']) + self.cval)
                        elif self.data_type == 'SD':
                            color_index = int(self.mval * float(pep['relsd']) + self.cval)
                        elif self.data_type == 'Uptake Diff':
                            color_index = int(self.mval * pep['uptakediff'] + self.cval)
                        elif self.data_type == 'Fractional Uptake':
                            color_index = int(self.mval * pep['fracuptake'] + self.cval)
                        elif self.data_type == 'Fractional SD (RSD)':
                            color_index = int(self.mval * pep['fracsd'] + self.cval)
                        elif self.data_type == 'Fractional Uptake Diff':
                            color_index = int(self.mval * pep['fracuptakediff'] + self.cval)
                        else:
                            print("Type Error")
                    if color_index < 0:
                        color_index = 0
                    elif color_index >= 100:
                        color_index = 99
                    color = self.rearranged_colors[color_index]
                    axis.add_patch(
                        PathPatch(Path(vertices, codes), facecolor=color, linewidth=1, edgecolor='grey'))
            elif self.map_type == 'Heat':
                # 2020-06-04 pep['Start'] - 1 changed to pep['Start'] - firstpos - 1
                # 2020-06-04 pep['End'] - 1 changed to pep['End'] - firstpos - 1
                for peptide in self.seq_nest2.keys():
                    pep = [n for n in self.pepdict if (n["Start"] == peptide[0] and n['End'] == peptide[1])][0]
                    start = min(self.seq_nest2[peptide])
                    end = max(self.seq_nest2[peptide])
                    vertices = []
                    codes = []
                    self.redundancy=1

                    stretch = int((end - self.firstpos+self.startspace-1) / self.line_length) - int((start - self.firstpos+self.startspace-1) / self.line_length)
                    if stretch > 0:
                        level = int((start - self.firstpos+self.startspace-1) / self.line_length) + 1
                        yheight1 = self.ymax - (level * (self.redundancy + self.spacing))+0.5
                        start50 = (start - self.firstpos+self.startspace-1) % self.line_length + self.xpad
                        codes += [Path.MOVETO] + [Path.LINETO] * 3
                        vertices += [(2*(self.line_length + self.xpad), 2*yheight1), (2*start50, 2*yheight1), (2*start50, 2*(yheight1 + 2)),
                                     (2*(self.line_length + self.xpad), 2*(yheight1 + 2))]
                        for n in range(1, stretch):
                            yheight = self.ymax - ((level + n) * (self.redundancy + self.spacing))+0.5
                            codes += [Path.MOVETO] + [Path.LINETO] * 3
                            vertices += [(2 * (self.line_length + self.xpad), 2 * yheight),
                                         (2 * self.xpad, 2 * yheight),
                                         (2 * self.xpad, 2 * (yheight + 0.75)),
                                         (2 * (self.line_length + self.xpad), 2 * (yheight + 0.75))]
                        yheight2 = self.ymax - ((level + stretch) * (self.redundancy + self.spacing)) + 0.5
                        stop50 = (end - self.firstpos+self.startspace-1) % self.line_length + 1
                        codes += [Path.MOVETO] + [Path.LINETO] * 3
                        vertices += [(2*self.xpad, 2*yheight2), (2*stop50, 2*yheight2), (2*stop50, 2*(yheight2 + 2)), (2*self.xpad, 2*(yheight2 + 2))]
                    else:
                        level = int((start - self.firstpos+self.startspace-1) / self.line_length) + 1
                        yheight = self.ymax - (level * (self.redundancy + self.spacing))+0.5
                        start50 = (start - self.firstpos+self.startspace-1) % self.line_length + self.xpad
                        stop50 = (end - self.firstpos+self.startspace-1) % self.line_length + 1
                        codes += [Path.MOVETO] + [Path.LINETO] * 3 + [Path.CLOSEPOLY]
                        vertices += [(2*start50, 2*yheight), (2*stop50, 2*yheight), (2*stop50, 2*(yheight + 2)),
                                     (2*start50, 2*(yheight + 2)),
                                     (0, 0)]
                    vertices = np.array(vertices, float)
                    if float(pep['MaxUptake']) == 0:
                        color_index = 0
                    else:
                        if self.data_type == 'Uptake':
                            color_index = int(self.mval * float(pep['reluptake']) + self.cval)
                        elif self.data_type == 'SD':
                            color_index = int(self.mval * float(pep['relsd']) + self.cval)
                        elif self.data_type == 'Uptake Diff':
                            color_index = int(self.mval * pep['uptakediff'] + self.cval)
                        elif self.data_type == 'Fractional Uptake':
                            color_index = int(self.mval * pep['fracuptake'] + self.cval)
                        elif self.data_type == 'Fractional SD (RSD)':
                            color_index = int(self.mval * pep['fracsd'] + self.cval)
                        elif self.data_type == 'Fractional Uptake Diff':
                            color_index = int(self.mval * pep['fracuptakediff'] + self.cval)
                        else:
                            print("Type Error")
                        if color_index < 0:
                            color_index = 0
                        elif color_index >= 100:
                            color_index = 99
                        color = self.rearranged_colors[color_index]
                    axis.add_patch(PathPatch(Path(vertices, codes), facecolor=color, edgecolor='grey'))

        # Add colorbar to GUI plot
        def axis_colorbar(self):
            ax_c = plt.axes([0.15, 0.065, 0.7, 0.02])
            cbar = self.fig.colorbar(matplotlib.cm.ScalarMappable(cmap=self.newmap), cax=ax_c, ticks = [0,1],
                                     cmap=self.newmap, orientation="horizontal")
            cbar.ax.set_xticklabels([round(self.color_values[0],2),round(self.color_values[1],2)])
            cbar.ax.tick_params(labelsize=20)

        # Add colorbar inline to the saved plot
        def inline_colorbar(self, axis):
            col_code = [Path.MOVETO] + [Path.LINETO] * 3 + [Path.CLOSEPOLY]
            col_vert = [(self.line_length - 25, -4), (self.line_length + 25, -4), (self.line_length + 25, 0),
                        (self.line_length - 25, 0), (self.line_length - 25, -4)]
            patch = PathPatch(Path(col_vert, col_code), facecolor='none', edgecolor='grey')
            axis.add_patch(patch)
            # Fixed 2020-06-04 for colorbar linelength issue
            im = axis.imshow([[0, 1], [0., 1]], interpolation='bilinear', aspect='auto',
                             extent=(self.line_length - 25, self.line_length + 25, -4, 0), cmap=self.newmap,
                             clip_path=patch, clip_on=True)
            axis.annotate(str(round(self.color_values[0],2)), (self.line_length - 30, -2),
                          color='black', fontsize=self.textsize, horizontalalignment='center',
                          verticalalignment='center', family='monospace')
            axis.annotate(str(round(self.color_values[1],2)), (self.line_length + 30, -2),
                          color='black', fontsize=self.textsize, horizontalalignment='center',
                          verticalalignment='center', family='monospace')
            im.set_clip_path(patch)

        # Draw map
        def draw(self, axis, figure):
            axis.set_xbound(0, self.line_length + self.xpad * 2)
            axis.set_xlim((0, self.line_length * 2))
            axis.set_ylim((-4, self.ymax * 2))
            figure.set_size_inches(12 * self.line_length / 50, 12 * self.ymax / 61)
            figure.set_dpi(72)

        # Add scrollbars to map
        def scrollbars(self):
            plt.draw()
            self.ax.set_xlim((0, 100))
            self.ax.set_ylim((self.ymax * 2 - 100, self.ymax * 2))
            self.fig.canvas.manager.window.geometry("720x720+0+0")
            ax_y = self.fig.add_axes([0.975, 0.05, 0.015, 0.9])
            scroll_y = Slider(ax_y, '', 100, self.ymax * 2, valinit=self.ymax * 2, orientation='vertical',
                              facecolor='white', edgecolor='black', linewidth=3)
            ax_x = self.fig.add_axes([0.05, 0.01, 0.9, 0.015])
            scroll_x = Slider(ax_x, '', 100, self.line_length * 2 + 2, valinit=100, orientation='horizontal',
                              facecolor='white', edgecolor='black', linewidth=3)

            def updatey(val):
                newval = scroll_y.val
                self.ax.set_ylim((newval - 100, newval))
                plt.draw()

            def updatex(val):
                newval = scroll_x.val
                self.ax.set_xlim((newval - 100, newval))
                plt.draw()

            scroll_y.on_changed(updatey)
            scroll_x.on_changed(updatex)

            self.scroll_y = scroll_y
            self.scroll_x = scroll_x
    class map_save():
        # Initializes plot and right click menu
        def __init__(self, map_type, parent):
            self.parent = parent
            self.map_type = map_type
            self.data_type = str(parent.representation.get())
            self.protein = str(parent.protein.get())
            self.state = str(parent.state1.get())
            self.state2 = str(parent.state2.get())
            self.exposure = str(parent.exposure.get())
            self.color = str(parent.color.get())
            self.reverse = int(parent.reverse.get())
            self.xs = int(parent.xs.get())
            self.cmaps = parent.cmaps
            self.data_range = str(parent.range.get())
            self.line_length = int(parent.xlim.get())
            self.corrected = parent.main.corrected
            self.obj = parent.main.state_obj[parent.main.csvfile]

            self.settings = {'Font':'Arial',
                             'Show Title':0,
                             'Title Size':32}
            self.peplist = {0: []}
            self.pepdict = []
            self.peptides = {}
            self.peptides[0] = []
            self.peptides[1] = []
            self.vertices = []
            self.codes = []
            self.textsize = 16
            self.xpad = 0

            # Create Figure & Axis and apply plot properties
            self.make_cmap()
            self.get_peptide_data()
            self.set_color_values()
            self.determine_levels()

            self.export()

        # Customizes colormap
        def make_cmap(self):
            num_colors = 102
            cm = plt.get_cmap(self.cmaps[self.color])
            colorlist = [cm(1. * i / num_colors) for i in range(num_colors)]
            if not self.reverse:
                colorlist = list(colorlist)
            else:
                colorlist = list(reversed(colorlist))
            self.rearranged_colors = colorlist[1:101]
            self.newmap = matplotlib.colors.LinearSegmentedColormap.from_list("newmap", self.rearranged_colors)
            self.color_values = ()
            self.cval = 0
            self.mval = 0

        # Sets the color range based on specified settings and the data range
        def set_color_values(self):
            if self.data_type == 'Uptake':
                values = [i['reluptake'] for i in self.pepdict]
                if self.data_range in ['Fit to Data', 'Full Range', 'Half Range']:
                    color_values = (min(values), max(values))
                else:
                    self.custom_min = float(self.parent.view.min_range.get())
                    self.custom_max = float(self.parent.view.max_range.get())
                    color_values = (self.custom_min, self.custom_max)
            elif self.data_type == 'SD':
                values = [i['relsd'] for i in self.pepdict]
                if self.data_range in ['Fit to Data', 'Full Range', 'Half Range']:
                    color_values = (min(values), max(values))
                else:
                    self.custom_min = float(self.parent.view.min_range.get())
                    self.custom_max = float(self.parent.view.max_range.get())
                    color_values = (self.custom_min, self.custom_max)
            elif self.data_type == 'Uptake Diff':
                values = [i['uptakediff'] for i in self.pepdict]
                if self.data_range in ['Fit to Data', 'Full Range', 'Half Range']:
                    color_values = (min(values), max(values))
                else:
                    self.custom_min = float(self.parent.view.min_range.get())
                    self.custom_max = float(self.parent.view.max_range.get())
                    color_values = (self.custom_min, self.custom_max)
            elif self.data_type == 'Fractional Uptake':
                values = [i['fracuptake'] for i in self.pepdict]
                if self.data_range in ['Fit to Data']:
                    color_values = (min(values), max(values))
                elif self.data_range in ['Full Range']:
                    color_values = (0, 1)
                elif self.data_range in ['Half Range']:
                    color_values = (0, 0.5)
                else:
                    self.custom_min = float(self.parent.view.min_range.get())
                    self.custom_max = float(self.parent.view.max_range.get())
                    color_values = (self.custom_min, self.custom_max)
            elif self.data_type == 'Fractional SD (RSD)':
                values = [i['fracsd'] for i in self.pepdict]
                if self.data_range in ['Fit to Data']:
                    color_values = (min(values), max(values))
                elif self.data_range in ['Full Range']:
                    color_values = (0, 1)
                elif self.data_range in ['Half Range']:
                    color_values = (0, 0.5)
                else:
                    self.custom_min = float(self.parent.view.min_range.get())
                    self.custom_max = float(self.parent.view.max_range.get())
                    color_values = (self.custom_min, self.custom_max)
            elif self.data_type == 'Fractional Uptake Diff':
                values = [i['fracuptakediff'] for i in self.pepdict]
                if self.data_range in ['Fit to Data']:
                    color_values = (min(values), max(values))
                elif self.data_range in ['Full Range']:
                    color_values = (-1, 1)
                elif self.data_range in ['Half Range']:
                    color_values = (-0.5, 0.5)
                else:
                    self.custom_min = float(self.parent.view.min_range.get())
                    self.custom_max = float(self.parent.view.max_range.get())
                    color_values = (self.custom_min, self.custom_max)
            else:
                print("Type Error")
            if color_values[0] == 0 and color_values[1] == 0:
                cval = 0
                mval = 0
            elif color_values[0] == 0 and color_values[1] != 0:
                cval = 0
                mval = 100 / color_values[1]
            elif color_values[1] == color_values[0]:
                print("Range Error")
                cval = 0
                mval = 0
            else:
                cval = float(100 / (-(color_values[1] / color_values[0]) + 1))
                mval = float((100 - cval) / color_values[1])
            self.color_values = color_values
            self.cval = cval
            self.mval = mval

        # Obtains peptide uptake information
        def get_peptide_data(self):
            # Doesn't work for AVE exposure
            if self.corrected == 'Yes':
                self.uptakekey = 'Uptake_corr'
                self.sdkey = 'Uptake_SD_corr'
            else:
                self.uptakekey = 'Uptake'
                self.sdkey = 'Uptake SD'
            # Get peptide list and exposure list
            if self.data_type in ['Uptake Diff', 'Fractional Uptake Diff']:
                if self.exposure == 'Ave':
                    exposures1 = [float(i) for i in self.obj.data_nest[self.protein][self.state].keys()]
                    exposures2 = [float(i) for i in self.obj.data_nest[self.protein][self.state2].keys()]
                    exposures = list(set(exposures1).intersection(exposures2))
                    exposure = min(exposures)
                else:
                    exposure = float(self.exposure)
                dict0 = sorted([i['PepID'] for i in self.obj.data_nest[self.protein][self.state][exposure] if
                                i['Fragment'] == ''])
                dict1 = sorted([i['PepID'] for i in self.obj.data_nest[self.protein][self.state2][exposure] if
                                i['Fragment'] == ''])
                id_list = list(set(dict0).intersection(dict1))
            else:
                if self.exposure == 'Ave':
                    exposures = [float(i) for i in self.obj.data_nest[self.protein][self.state].keys()]
                    exposure = exposures[0]
                    id_list = sorted([i['PepID'] for i in self.obj.data_nest[self.protein][self.state][exposure] if
                                      i['Fragment'] == ''])
                else:
                    exposure = float(self.exposure)
                    id_list = sorted([i['PepID'] for i in self.obj.data_nest[self.protein][self.state][exposure] if
                                      i['Fragment'] == ''])

            # Get data for all peptides
            if self.exposure == 'Ave':
                if 0 in exposures:
                    exposures.remove(0)
                for pepid in id_list:
                    elist = []
                    for e in sorted(exposures):
                        elist.extend(
                            [n for n in self.obj.data_nest[self.protein][self.state][e] if n["PepID"] == pepid])
                    reluptake = np.mean([float(i[self.uptakekey]) for i in elist])
                    relsd = np.mean([float(i[self.sdkey]) for i in elist])
                    if elist[0]['MaxUptake'] == 0:
                        fracuptake = 0
                    else:
                        fracuptake = np.mean([float(i[self.uptakekey]) / float(i['MaxUptake']) for i in elist])
                    if float(elist[0][self.uptakekey]) == 0:
                        fracsd = 0
                    else:
                        fracsd = np.mean([float(i[self.sdkey]) / float(i[self.uptakekey]) for i in elist])
                    if self.data_type in ['Uptake Diff', 'Fractional Uptake Diff']:
                        e2list = []
                        for e in sorted(exposures):
                            e2list.extend([n for n in self.obj.data_nest[self.protein][self.state2][e] if
                                           n["PepID"] == pepid])
                        uptakediff = np.mean(
                            [float(j[self.uptakekey]) - float(i[self.uptakekey]) for j in e2list for i in elist])
                        if e2list[0]['MaxUptake'] == 0:
                            fracuptakediff = 0
                        else:
                            fracuptakediff = np.mean([float(j[self.uptakekey]) / float(j['MaxUptake']) -
                                                      float(i[self.uptakekey]) / float(i['MaxUptake'])
                                                      for j in e2list for i in elist])
                        self.pepdict.append(
                            {'Start': elist[0]['Start'], 'End': elist[0]['End'], 'Fragment': elist[0]['Fragment'],
                             'Modification': elist[0]['Modification'], 'MaxUptake': elist[0]['MaxUptake'],
                             'reluptake': reluptake, 'uptakediff': uptakediff, 'fracuptake': fracuptake,
                             'fracuptakediff': fracuptakediff, 'relsd': relsd, 'fracsd': fracsd})
                    else:
                        self.pepdict.append(
                            {'Start': elist[0]['Start'], 'End': elist[0]['End'], 'Fragment': elist[0]['Fragment'],
                             'Modification': elist[0]['Modification'], 'MaxUptake': elist[0]['MaxUptake'],
                             'reluptake': reluptake, 'fracuptake': fracuptake, 'relsd': relsd, 'fracsd': fracsd})
            else:
                for pepid in id_list:
                    pep = \
                    [n for n in self.obj.data_nest[self.protein][self.state][exposure] if n["PepID"] == pepid][0]
                    reluptake = float(pep[self.uptakekey])
                    relsd = float(pep[self.sdkey])
                    if pep['MaxUptake'] == 0:
                        fracuptake = 0
                    else:
                        fracuptake = float(pep[self.uptakekey]) / float(pep['MaxUptake'])
                    if float(pep[self.uptakekey]) == 0:
                        fracsd = 0
                    else:
                        fracsd = float(pep[self.sdkey]) / float(pep[self.uptakekey])
                    if self.data_type in ['Uptake Diff', 'Fractional Uptake Diff']:
                        pep2 = \
                        [n for n in self.obj.data_nest[self.protein][self.state2][exposure] if n["PepID"] == pepid][
                            0]
                        uptakediff = float(pep2[self.uptakekey]) - float(pep[self.uptakekey])
                        if pep2['MaxUptake'] == 0:
                            fracuptakediff = 0
                        else:
                            fracuptakediff = (float(pep2[self.uptakekey]) - float(pep[self.uptakekey])) / float(
                                pep['MaxUptake'])
                        self.pepdict.append({'Start': pep['Start'], 'End': pep['End'], 'Fragment': pep['Fragment'],
                                             'Modification': pep['Modification'], 'MaxUptake': pep['MaxUptake'],
                                             'reluptake': reluptake, 'uptakediff': uptakediff,
                                             'fracuptake': fracuptake,
                                             'fracuptakediff': fracuptakediff, 'relsd': relsd, 'fracsd': fracsd})
                    else:
                        self.pepdict.append({'Start': pep['Start'], 'End': pep['End'], 'Fragment': pep['Fragment'],
                                             'Modification': pep['Modification'], 'MaxUptake': pep['MaxUptake'],
                                             'reluptake': reluptake, 'fracuptake': fracuptake, 'relsd': relsd,
                                             'fracsd': fracsd})
            if self.map_type == 'Heat':
                if self.data_type in ['Fractional Uptake Diff', 'Uptake Diff']:
                    if self.exposure == 'Ave':
                        exposures1 = [float(i) for i in
                                      self.obj.data_nest[self.protein][self.state].keys()]
                        exposures2 = [float(i) for i in self.obj.data_nest[self.protein][self.state2].keys()]
                        exposures = list(set(exposures1).intersection(exposures2))
                        exposure = min(exposures)
                        self.obj.assignData(self.protein, self.state, 'Ave', self.state2)
                    else:
                        exposure = float(self.exposure)
                        self.obj.assignData(self.protein, self.state, exposure, self.state2)
                else:

                    if self.exposure == 'Ave':
                        exposures = [float(i) for i in self.obj.data_nest[self.protein][self.state].keys()]
                        exposure = exposures[0]
                        self.obj.assignData(self.protein, self.state, 'Ave')
                    else:
                        exposure = float(self.exposure)
                        self.obj.assignData(self.protein, self.state, exposure)
                self.seq_nest2 = {}
                for item in self.obj.seq_nest:
                    if item['assigned'] == 1:
                        if (item['pepmin'], item['pepmax']) in self.seq_nest2.keys():
                            self.seq_nest2[(item['pepmin'], item['pepmax'])].append(item['pos'])
                        else:
                            self.seq_nest2[(item['pepmin'], item['pepmax'])] = [item['pos']]

            self.sequence = self.obj.sequences[self.protein]
            if not self.xs:
                self.firstpos = min([i for i, j in enumerate(self.sequence) if j != "X"])
                self.startspace = 0
            else:
                self.firstpos = min([i for i, j in enumerate(self.sequence) if j != "X"])
                self.startspace = (self.firstpos) % self.line_length
            self.sequence = self.sequence[self.firstpos:]

        # Determines the height of the map depending on the line width, sequence length, and peptide coverage
        def determine_levels(self):
            self.seqlist = list(self.sequence)
            self.seqlength = len(self.seqlist) + self.startspace
            self.levels = int(self.seqlength / self.line_length) + 1
            self.level_dict = {}
            for pep in sorted(self.pepdict, key=itemgetter('Start')):
                n = 0
                while set(self.peplist[n]).intersection(list(range(pep['Start'], pep['End'] + 1))) != set([]):
                    n = n + 1
                    if n not in self.peplist.keys():
                        self.peplist[n] = []
                self.peplist[n] += range(pep['Start'], pep['End'] + 1)
                # Added modification to peptide descriptor
                self.level_dict[(pep['Start'], pep['End'], pep['Fragment'], pep['Modification'])] = n

            if self.map_type == 'Coverage':
                self.redun_dict = {}
                self.spacing = 4
                self.ymax = 2
                for level in range(0, self.levels):
                    self.redun_dict[level + 1] = 0
                    range_min = (level * self.line_length) + 1
                    range_max = (level + 1) * self.line_length
                    range_list = list(range(range_min, range_max + 1))
                    for n in self.peplist.keys():
                        if set([resn - self.firstpos + self.startspace for resn in self.peplist[n] if
                                range_min <= (resn - self.firstpos + self.startspace) <= range_max]).intersection(
                                range_list) != set([]):
                            self.redun_dict[level + 1] = n + 1
                    self.ymax += (self.redun_dict[level + 1] + self.spacing)
            elif self.map_type == 'Heat':
                self.redun_dict = {}
                for level in range(0, self.levels):
                    self.redun_dict[level + 1] = 1
                self.spacing = 4
                self.ymax = (self.levels * (1 + self.spacing)) + 2
            else:
                print('type error')

        # Saves plot to file
        def export(self):
            fig = plt.figure()
            ax = fig.add_axes([0.01, 0.01, 0.98, 0.98])
            fig.canvas.toolbar.pack_forget()
            fig.set_dpi(300)
            ax.axis("off")
            self.create_axes_dnx(ax)
            self.add_peptides_dnx(ax)
            self.inline_colorbar(ax)
            self.draw(ax, fig)
            plt.draw()
            if self.data_type in ['Uptake Diff', 'Fractional Uptake Diff']:
                savefile = FileDialog.asksaveasfilename(
                    title="Choose save location",
                    initialfile=self.protein + "_" + str(self.state) + "-" + str(self.state2) + "_" + str(self.exposure)+"_Map.png",
                    initialdir=self.parent.main.dir,
                    filetypes=(("PNG", '*.png'), ("JPG", "*.jpg"), ("TIF", "*.tif"), ("all files", "*.*")))
            else:
                savefile = FileDialog.asksaveasfilename(
                    title="Choose save location",
                    initialfile=self.protein + "_" + str(self.state) + "_" + str(self.exposure) + "_Map.png",
                    initialdir=self.parent.main.dir,
                    filetypes=(("PNG", '*.png'), ("JPG", "*.jpg"), ("TIF", "*.tif"), ("all files", "*.*")))
            if split(savefile)[0] != "":
                self.parent.main.dir = split(savefile)[0]
            if str(splitext(split(savefile)[1])[1]).lower() not in ['.png','.jpg','.tif']:
                savefile = savefile + '.png'
            fig.savefig(savefile, dpi=300, transparent=True)

        # Initializes axis in DynamX Style
        def create_axes_dnx(self, axis):
            self.ycurlevel = self.ymax - 1.5
            for n in range(0, self.levels):
                level = n + 1
                redundancy = self.redun_dict[level]
                # If sequence is longer than the current level:
                if level == 1:
                    for m in range(self.startspace, (n + 1) * self.line_length):
                        x = m % self.line_length
                        axis.annotate(self.seqlist[m - self.startspace],
                                      (2 * (x + self.xpad + 0.5), 2 * (self.ycurlevel + 0.5)),
                                      color='black', fontsize=self.textsize, horizontalalignment='center',
                                      verticalalignment='center',
                                      family='monospace')
                        if (m + 1 + self.firstpos - self.startspace) % 10 == 0:
                            # 2020-06-04 m + 1 changed to m + firstpos + 1 to start map at first coverage
                            axis.annotate(str(m + self.firstpos - self.startspace + 1),
                                          (2 * (x + self.xpad + 0.5), 2 * (self.ycurlevel - 0.5)),
                                          color='black', fontsize=self.textsize * 0.75,
                                          horizontalalignment='center',
                                          verticalalignment='center', family='monospace')
                elif (self.seqlength - self.line_length * n) > self.line_length:
                    for m in range(n * self.line_length, (n + 1) * self.line_length):
                        x = m % self.line_length
                        axis.annotate(self.seqlist[m - self.startspace],
                                      (2 * (x + self.xpad + 0.5), 2 * (self.ycurlevel + 0.5)),
                                      color='black', fontsize=self.textsize, horizontalalignment='center',
                                      verticalalignment='center',
                                      family='monospace')
                        if (m + 1 + self.firstpos - self.startspace) % 10 == 0:
                            # 2020-06-04 m + 1 changed to m + firstpos + 1 to start map at first coverage
                            axis.annotate(str(m + self.firstpos - self.startspace + 1),
                                          (2 * (x + self.xpad + 0.5), 2 * (self.ycurlevel - 0.5)),
                                          color='black', fontsize=self.textsize * 0.75,
                                          horizontalalignment='center',
                                          verticalalignment='center', family='monospace')
                # Conditions for final line
                else:
                    for m in range(n * self.line_length, self.seqlength):
                        x = m % self.line_length
                        axis.annotate(self.seqlist[m - self.startspace],
                                      (2 * (x + self.xpad + 0.5), 2 * (self.ycurlevel + 0.5)),
                                      color='black', fontsize=self.textsize, horizontalalignment='center',
                                      verticalalignment='center',
                                      family='monospace')
                        if (m + 1 + self.firstpos) % 10 == 0:
                            # 2020-06-04 m + 1 changed to m + firstpos + 1 to start map at first coverage
                            axis.annotate(str(m + self.firstpos - self.startspace + 1),
                                          (2 * (x + self.xpad + 0.5), 2 * (self.ycurlevel - 0.5)),
                                          color='black', fontsize=self.textsize * 0.75,
                                          horizontalalignment='center',
                                          verticalalignment='center', family='monospace')
                self.ycurlevel = self.ycurlevel - (redundancy + self.spacing)

        # Add Bars to the Map in DynamX Style
        def add_peptides_dnx(self, axis):
            # Creating Boxes
            if self.map_type == 'Coverage':
                # 2020-06-04 pep['Start'] - 1 changed to pep['Start'] - firstpos - 1
                # 2020-06-04 pep['End'] - 1 changed to pep['End'] - firstpos - 1
                for pep in self.pepdict:
                    vertices = []
                    codes = []
                    stretch = int((pep['End'] - self.firstpos + self.startspace - 1) / self.line_length) - int(
                        (pep['Start'] - self.firstpos + self.startspace - 1) / self.line_length)
                    peplevel = self.level_dict[(pep['Start'], pep['End'], pep['Fragment'], pep['Modification'])]
                    if stretch > 0:
                        level = int((pep['Start'] - self.firstpos + self.startspace - 1) / self.line_length) + 1
                        yheight1 = -peplevel + self.ymax - (
                                sum([self.redun_dict[j] for j in range(1, level)]) + (level) * self.spacing) + 0.5
                        start50 = (pep[
                                       'Start'] - self.firstpos + self.startspace - 1) % self.line_length + self.xpad
                        codes += [Path.MOVETO] + [Path.LINETO] * 3
                        vertices += [(2 * (self.line_length + self.xpad), 2 * yheight1),
                                     (2 * start50, 2 * yheight1),
                                     (2 * start50, 2 * (yheight1 + 0.75)),
                                     (2 * (self.line_length + self.xpad), 2 * (yheight1 + 0.75))]

                        for n in range(1, stretch):
                            yheight = -peplevel + self.ymax - (
                                        sum([self.redun_dict[j] for j in range(1, level + n)]) +
                                        (level + n) * self.spacing) + 0.5
                            codes += [Path.MOVETO] + [Path.LINETO] * 3
                            vertices += [(2 * (self.line_length + self.xpad), 2 * yheight),
                                         (2 * self.xpad, 2 * yheight),
                                         (2 * self.xpad, 2 * (yheight + 0.75)),
                                         (2 * (self.line_length + self.xpad), 2 * (yheight + 0.75))]
                        yheight2 = -peplevel + self.ymax - (
                                sum([self.redun_dict[j] for j in range(1, level + stretch)]) + (
                                    level + stretch) * self.spacing) + 0.5
                        stop50 = (pep['End'] - self.firstpos + self.startspace - 1) % self.line_length + 1
                        codes += [Path.MOVETO] + [Path.LINETO] * 3
                        vertices += [(2 * self.xpad, 2 * yheight2), (2 * stop50, 2 * yheight2),
                                     (2 * stop50, 2 * (yheight2 + 0.75)),
                                     (2 * self.xpad, 2 * (yheight2 + 0.75))]
                    else:
                        level = int((pep['Start'] - self.firstpos + self.startspace - 1) / self.line_length) + 1
                        yheight = -peplevel + self.ymax - (
                                sum([self.redun_dict[j] for j in range(1, level)]) + level * self.spacing) + 0.5
                        start50 = (pep[
                                       'Start'] - self.firstpos + self.startspace - 1) % self.line_length + self.xpad
                        stop50 = (pep['End'] - self.firstpos + self.startspace - 1) % self.line_length + 1
                        codes += [Path.MOVETO] + [Path.LINETO] * 3 + [Path.CLOSEPOLY]
                        vertices += [(2 * start50, 2 * yheight), (2 * stop50, 2 * yheight),
                                     (2 * stop50, 2 * (yheight + 0.75)),
                                     (2 * start50, 2 * (yheight + 0.75)),
                                     (0, 0)]
                    vertices = array(vertices, float)
                    if float(pep['MaxUptake']) == 0:
                        color_index = 0
                    else:
                        if self.data_type == 'Uptake':
                            color_index = int(self.mval * float(pep['reluptake']) + self.cval)
                        elif self.data_type == 'SD':
                            color_index = int(self.mval * float(pep['relsd']) + self.cval)
                        elif self.data_type == 'Uptake Diff':
                            color_index = int(self.mval * pep['uptakediff'] + self.cval)
                        elif self.data_type == 'Fractional Uptake':
                            color_index = int(self.mval * pep['fracuptake'] + self.cval)
                        elif self.data_type == 'Fractional SD (RSD)':
                            color_index = int(self.mval * pep['fracsd'] + self.cval)
                        elif self.data_type == 'Fractional Uptake Diff':
                            color_index = int(self.mval * pep['fracuptakediff'] + self.cval)
                        else:
                            print("Type Error")
                    if color_index < 0:
                        color_index = 0
                    elif color_index >= 100:
                        color_index = 99
                    color = self.rearranged_colors[color_index]
                    axis.add_patch(
                        PathPatch(Path(vertices, codes), facecolor=color, linewidth=1, edgecolor='grey'))
            elif self.map_type == 'Heat':
                # 2020-06-04 pep['Start'] - 1 changed to pep['Start'] - firstpos - 1
                # 2020-06-04 pep['End'] - 1 changed to pep['End'] - firstpos - 1
                for peptide in self.seq_nest2.keys():
                    pep = [n for n in self.pepdict if (n["Start"] == peptide[0] and n['End'] == peptide[1])][0]
                    start = min(self.seq_nest2[peptide])
                    end = max(self.seq_nest2[peptide])
                    vertices = []
                    codes = []
                    self.redundancy = 1

                    stretch = int((end - self.firstpos + self.startspace - 1) / self.line_length) - int(
                        (start - self.firstpos + self.startspace - 1) / self.line_length)
                    if stretch > 0:
                        level = int((start - self.firstpos + self.startspace - 1) / self.line_length) + 1
                        yheight1 = self.ymax - (level * (self.redundancy + self.spacing)) + 0.5
                        start50 = (start - self.firstpos + self.startspace - 1) % self.line_length + self.xpad
                        codes += [Path.MOVETO] + [Path.LINETO] * 3
                        vertices += [(2 * (self.line_length + self.xpad), 2 * yheight1),
                                     (2 * start50, 2 * yheight1), (2 * start50, 2 * (yheight1 + 2)),
                                     (2 * (self.line_length + self.xpad), 2 * (yheight1 + 2))]
                        for n in range(1, stretch):
                            yheight = self.ymax - ((level + n) * (self.redundancy + self.spacing)) + 0.5
                            codes += [Path.MOVETO] + [Path.LINETO] * 3
                            vertices += [(2 * (self.line_length + self.xpad), 2 * yheight),
                                         (2 * self.xpad, 2 * yheight),
                                         (2 * self.xpad, 2 * (yheight + 0.75)),
                                         (2 * (self.line_length + self.xpad), 2 * (yheight + 0.75))]
                        yheight2 = self.ymax - ((level + stretch) * (self.redundancy + self.spacing)) + 0.5
                        stop50 = (end - self.firstpos + self.startspace - 1) % self.line_length + 1
                        codes += [Path.MOVETO] + [Path.LINETO] * 3
                        vertices += [(2 * self.xpad, 2 * yheight2), (2 * stop50, 2 * yheight2),
                                     (2 * stop50, 2 * (yheight2 + 2)), (2 * self.xpad, 2 * (yheight2 + 2))]
                    else:
                        level = int((start - self.firstpos + self.startspace - 1) / self.line_length) + 1
                        yheight = self.ymax - (level * (self.redundancy + self.spacing)) + 0.5
                        start50 = (start - self.firstpos + self.startspace - 1) % self.line_length + self.xpad
                        stop50 = (end - self.firstpos + self.startspace - 1) % self.line_length + 1
                        codes += [Path.MOVETO] + [Path.LINETO] * 3 + [Path.CLOSEPOLY]
                        vertices += [(2 * start50, 2 * yheight), (2 * stop50, 2 * yheight),
                                     (2 * stop50, 2 * (yheight + 2)),
                                     (2 * start50, 2 * (yheight + 2)),
                                     (0, 0)]
                    vertices = np.array(vertices, float)
                    if float(pep['MaxUptake']) == 0:
                        color_index = 0
                    else:
                        if self.data_type == 'Uptake':
                            color_index = int(self.mval * float(pep['reluptake']) + self.cval)
                        elif self.data_type == 'SD':
                            color_index = int(self.mval * float(pep['relsd']) + self.cval)
                        elif self.data_type == 'Uptake Diff':
                            color_index = int(self.mval * pep['uptakediff'] + self.cval)
                        elif self.data_type == 'Fractional Uptake':
                            color_index = int(self.mval * pep['fracuptake'] + self.cval)
                        elif self.data_type == 'Fractional SD (RSD)':
                            color_index = int(self.mval * pep['fracsd'] + self.cval)
                        elif self.data_type == 'Fractional Uptake Diff':
                            color_index = int(self.mval * pep['fracuptakediff'] + self.cval)
                        else:
                            print("Type Error")
                        if color_index < 0:
                            color_index = 0
                        elif color_index >= 100:
                            color_index = 99
                        color = self.rearranged_colors[color_index]
                    axis.add_patch(PathPatch(Path(vertices, codes), facecolor=color, edgecolor='grey'))

        # Add colorbar inline to the saved plot
        def inline_colorbar(self, axis):
            col_code = [Path.MOVETO] + [Path.LINETO] * 3 + [Path.CLOSEPOLY]
            col_vert = [(self.line_length - 25, -4), (self.line_length + 25, -4), (self.line_length + 25, 0),
                        (self.line_length - 25, 0), (self.line_length - 25, -4)]
            patch = PathPatch(Path(col_vert, col_code), facecolor='none', edgecolor='grey')
            axis.add_patch(patch)
            # Fixed 2020-06-04 for colorbar linelength issue
            im = axis.imshow([[0, 1], [0., 1]], interpolation='bilinear', aspect='auto',
                             extent=(self.line_length - 25, self.line_length + 25, -4, 0), cmap=self.newmap,
                             clip_path=patch, clip_on=True)
            axis.annotate(str(round(self.color_values[0], 2)), (self.line_length - 30, -2),
                          color='black', fontsize=self.textsize, horizontalalignment='center',
                          verticalalignment='center', family='monospace')
            axis.annotate(str(round(self.color_values[1], 2)), (self.line_length + 30, -2),
                          color='black', fontsize=self.textsize, horizontalalignment='center',
                          verticalalignment='center', family='monospace')
            im.set_clip_path(patch)

        # Draw map
        def draw(self, axis, figure):
            axis.set_xbound(0, self.line_length + self.xpad * 2)
            axis.set_xlim((0, self.line_length * 2))
            axis.set_ylim((-4, self.ymax * 2))
            figure.set_size_inches(12 * self.line_length / 50, 12 * self.ymax / 61)
            figure.set_dpi(72)


class Butterfly():
    '''
    Butterfly Plot Generator Window Controller
    Widgets
        Protein Dropdown
        State Dropdown
        Exposure Dropdown
        Wings Button - traditional butterfly plot
        Difference Button - difference between two states
    Determine Desired settings before plotting map in new window
    '''
    # Initializes Generator window and right click menu, Binds buttons, Defines Tk Variables
    def __init__(self, main):
        self.oldval = 0
        self.style = ttk.Style()
        self.main = main

        # Initialize Tk window
        self.view = View.Butterfly(main)
        main.widget_dict['Butterfly'] = self.view.top
        self.main.view.window_menu.add_command(label='Butterfly', command=lambda: self.main.focus(self.view.top))

        # Get data from root
        self.states = main.statesvar.get()
        self.proteins = main.proteinsvar.get()
        self.exposures = ['All']
        self.exposures.extend(main.exposuresvar.get())

        # Initialize Tk variables
        self.protein = Tk.StringVar()
        self.state1 = Tk.StringVar()
        self.state2 = Tk.StringVar()
        self.exposure = Tk.StringVar()

        # Load Data and Bindings into Tk
        self.view.combobox_protein.configure(values=self.proteins)
        self.view.combobox_protein.configure(textvariable=self.protein)
        self.view.combobox_protein.set(self.proteins[0])
        self.view.combobox_state1.configure(values=self.states)
        self.view.combobox_state1.configure(textvariable=self.state1)
        self.view.combobox_state1.set(self.states[0])
        self.view.combobox_state2.configure(values=self.states)
        self.view.combobox_state2.configure(textvariable=self.state2)
        self.view.combobox_exposure.configure(values=self.exposures)
        self.view.combobox_exposure.configure(textvariable=self.exposure)
        self.view.combobox_exposure.set('All')
        self.view.button_wings.configure(command=lambda: self.wings(self))
        self.view.button_diff.configure(command=lambda: self.diff(self))
        self.view.button_help.configure(command=lambda: self.main.help('6) Butterfly Plot'))
        self.view.top.protocol("WM_DELETE_WINDOW", lambda: main.on_closing(self.view.top, 'Butterfly'))

    # Traditional Butterfly mapping class
    class wings():
        # Initializes plot
        def __init__(self, parent):
            plt.ion()
            self.parent = parent
            self.settings = {'Font': 'Arial',
                             'Show Title': True,
                             'Title Size': 24,
                             'Show X Title': True,
                             'Show Y Title': True,
                             'Axis Label Size': 20}

            self.fig, self.ax = plt.subplots()

            self.getData()
            self.createFigure()
            plt.show()
            plt.get_current_fig_manager().window.geometry(str(72*10)+"x"+str(72*10)+"+0+0")
            plt.get_current_fig_manager().window.minsize(72*8, 72*8)
            self.c = self.context(self.fig.canvas.get_tk_widget(), self)

        # Obtains the data for the plot
        def getData(self):
            protein = str(self.parent.protein.get())
            state1 = str(self.parent.state1.get())
            state2 = str(self.parent.state2.get())
            obj = self.parent.main.state_obj[self.parent.main.csvfile]
            exposures = [float(i) for i in obj.data_nest[protein][state1].keys()]
            if self.parent.main.corrected == 'Yes':
                uptakekey = 'Uptake_corr'
                sdkey = 'Uptake_SD_corr'
            else:
                uptakekey = 'Uptake'
                sdkey = 'Uptake SD'

            obj.assignData(protein, state1, 0, state2)
            obj1 = obj.data_nest[protein][state1]
            obj2 = obj.data_nest[protein][state2]
            e_dict = {0:{},1:{}}
            for e in exposures:
                e_dict[0][e] = []
                e_dict[1][e] = []

            for item in obj.seq_nest:
                for e in exposures:
                    if item['assigned']:
                        pep = [i for i in obj1[e] if i['Start'] == item['pepmin'] and i['End'] == item['pepmax'] and
                               i['Fragment'] == item['fragment']][0]
                        e_dict[0][e].append((item['pos'],float(pep[uptakekey]) / float(pep['MaxUptake'])))
                        pep = [i for i in obj2[e] if i['Start'] == item['pepmin'] and i['End'] == item['pepmax'] and
                               i['Fragment'] == item['fragment']][0]
                        e_dict[1][e].append((item['pos'], float(pep[uptakekey]) / float(pep['MaxUptake'])))
            self.e_dict = e_dict
            self.protein = protein
            self.state1 = state1
            self.state2 = state2

        # Creates the plot
        def createFigure(self):
            print(self.settings)
            font = self.settings['Font']
            titlesize = self.settings['Title Size']
            if self.settings['Show Title']:
                self.ax.set_title('Butterfly Plot', fontname=font, fontsize=titlesize)
            axissize = self.settings['Axis Label Size']
            if self.settings['Show X Title']:
                self.ax.set_xlabel('Protein Sequence', fontname=font, fontsize=axissize)
            if self.settings['Show Y Title']:
                self.ax.set_ylabel('Fractional Uptake Difference (Da)', fontname=font, fontsize=axissize)
            for tick in self.ax.get_xticklabels():
                tick.set_fontname(font)
                tick.set_fontsize(16)
            for tick in self.ax.get_yticklabels():
                tick.set_fontname(font)
                tick.set_fontsize(16)
            num_colors = 102
            cm = plt.get_cmap('bwr')
            colorlist = [cm(1. * i / num_colors) for i in range(num_colors)]
            n=0
            y_val_list = list()
            exp = self.parent.exposure.get()
            exposures = list(self.e_dict[0].keys())
            if exp == 'All':
                for e in sorted(exposures):
                    c = int(1 + n * 100/(len(exposures)-1))
                    n=n+1
                    self.ax.plot(np.array([i[0] for i in self.e_dict[0][e]]), np.array([i[1] for i in self.e_dict[0][e]]), label = str(e)+' min', color=colorlist[c])
                    self.ax.plot(np.array([i[0] for i in self.e_dict[1][e]]), np.array([-i[1] for i in self.e_dict[1][e]]), color=colorlist[c])
                    y_val_list.extend([i[1] for i in self.e_dict[0][e]])
                    y_val_list.extend([i[1] for i in self.e_dict[1][e]])
                self.ax.set_ylim((-max(y_val_list),max(y_val_list)))
            else:
                e = float(exp)
                self.ax.plot(np.array([i[0] for i in self.e_dict[0][float(e)]]), np.array([i[1] for i in self.e_dict[0][float(e)]]), label = str(e)+' min', color='black')
                self.ax.plot(np.array([i[0] for i in self.e_dict[1][float(e)]]), np.array([-i[1] for i in self.e_dict[1][float(e)]]), color='black')
                y_val_list.extend([i[1] for i in self.e_dict[0][float(e)]])
                y_val_list.extend([i[1] for i in self.e_dict[1][float(e)]])
                self.ax.set_ylim((-max(y_val_list), max(y_val_list)))
            self.leg = self.ax.legend()
            self.leg.set_draggable(True)

        # Saves the plot to a file
        def export(self):
            savefile = FileDialog.asksaveasfilename(
                title="Choose save location",
                initialfile=self.protein + "_" + str(self.state1) + "-" + str(self.state2) + "_Map.png",
                initialdir=expanduser('~'),
                filetypes=(("PNG", '*.png'), ("all files", "*.*")))
            self.fig.savefig(savefile, dpi=300, transparent=True)

        # Updates the figure based on changed settings
        def update(self, key, value):
            self.settings[key] = value
            self.fig.clear()
            self.ax = self.fig.subplots()
            self.createFigure()
            plt.draw()

        # Context menu - why is this here?
        class context():
            # Initializes context menu
            def __init__(self, frame, main, *args, **kwargs):
                self.main = main
                fonts = sorted(set([f.name for f in matplotlib.font_manager.fontManager.ttflist]))

                popupMenu = Tk.Menu(frame, tearoff=0)
                fontMenu = Tk.Menu(popupMenu, tearoff=0)
                titleMenu = Tk.Menu(popupMenu, tearoff=0)
                axisMenu = Tk.Menu(popupMenu, tearoff=0)

                popupMenu.add_cascade(label="Font", menu=fontMenu)
                popupMenu.add_cascade(label="Title", menu=titleMenu)
                popupMenu.add_cascade(label="Axis", menu=axisMenu)
                for i in fonts:
                    if i == 'Arial':
                        fontMenu.add_command(label=i, command=lambda i=i: main.update('Font', i), background='green')
                    else:
                        fontMenu.add_command(label=i, command=lambda i=i: main.update('Font', i))
                self.title_val = Tk.IntVar()
                titleMenu.add_command(label="Show/Hide Title", command=lambda: main.update("Show Title", not main.settings['Show Title']))
                self.title_val.set(True)

                stMenu = Tk.Menu(titleMenu, tearoff=0)
                titleMenu.add_cascade(label="Title Size", menu=stMenu)
                for i in list(range(16, 65, 4)):
                    stMenu.add_command(label=i, command=lambda i=i: main.update('Title Size', i))

                axisMenu.add_command(label="Show/Hide X Title", command=lambda: main.update("Show X Title", not main.settings['Show X Title']))

                axisMenu.add_command(label="Show/Hide Y Title", command=lambda: main.update("Show Y Title", not main.settings['Show Y Title']))

                sMenu = Tk.Menu(axisMenu, tearoff=0)
                axisMenu.add_cascade(label="Y Axis Size", menu=sMenu)
                for i in list(range(16, 65, 4)):
                    sMenu.add_command(label=i, command=lambda i=i: main.update('Axis Label Size', i))

                popupMenu.add_command(label='Save', command=lambda: main.export())
                popupMenu.add_command(label='Close', command=lambda: plt.close(main.fig))

                if platform == 'darwin':
                    frame.bind("<Button-2>", lambda event: popupMenu.post(event.x_root, event.y_root))
                else:
                    frame.bind("<Button-3>", lambda event: popupMenu.post(event.x_root, event.y_root))
                self.main = main

    # Difference plot mapping class
    class diff():
        # Initializes plot
        def __init__(self, parent):
            plt.ion()
            self.parent = parent
            self.settings = {'Font':'Arial',
                             'Show Title':True,
                             'Title Size':24,
                             'Show X Title':True,
                             'Show Y Title':True,
                             'Axis Label Size': 24}

            self.fig, self.ax = plt.subplots()
            self.c = self.context(self.fig.canvas.get_tk_widget(), self)

            self.getData()
            self.createFigure()
            plt.show()
            plt.get_current_fig_manager().window.geometry(str(72 * 10) + "x" + str(72 * 10) + "+0+0")
            plt.get_current_fig_manager().window.minsize(72 * 8, 72 * 8)

        # Obtains the data for the plot
        def getData(self):
            protein = str(self.parent.protein.get())
            state1 = str(self.parent.state1.get())
            state2 = str(self.parent.state2.get())
            obj = self.parent.main.state_obj[self.parent.main.csvfile]
            exposures = [float(i) for i in obj.data_nest[protein][state1].keys()]
            obj.assignData(protein, state1, 0, state2)
            obj1 = obj.data_nest[protein][state1]
            obj2 = obj.data_nest[protein][state2]
            if self.parent.main.corrected == 'Yes':
                uptakekey = 'Uptake_corr'
                sdkey = 'Uptake_SD_corr'
            else:
                uptakekey = 'Uptake'
                sdkey = 'Uptake SD'
            e_dict = {0: {}, 1: {}}
            for e in exposures:
                e_dict[0][e] = []
                e_dict[1][e] = []
            for item in obj.seq_nest:
                for e in exposures:
                    if item['assigned']:
                        pep = [i for i in obj1[e] if i['Start'] == item['pepmin'] and i['End'] == item['pepmax'] and
                               i['Fragment'] == item['fragment']][0]
                        e_dict[0][e].append((item['pos'], float(pep[uptakekey]) / float(pep['MaxUptake'])))
                        pep = [i for i in obj2[e] if i['Start'] == item['pepmin'] and i['End'] == item['pepmax'] and
                               i['Fragment'] == item['fragment']][0]
                        e_dict[1][e].append((item['pos'], float(pep[uptakekey]) / float(pep['MaxUptake'])))
            self.e_dict = e_dict
            self.protein = protein
            self.state1 = state1
            self.state2 = state2

        # Creates the plot
        def createFigure(self):
            print(self.settings)
            font = self.settings['Font']
            titlesize = self.settings['Title Size']
            if self.settings['Show Title'] == 1:
                self.ax.set_title('Difference Butterfly Plot', fontname = font, fontsize = titlesize)
            axissize = self.settings['Axis Label Size']
            if self.settings['Show X Title'] == 1:
                self.ax.set_xlabel('Protein Sequence', fontname = font, fontsize = axissize)
            if self.settings['Show Y Title'] == 1:
                self.ax.set_ylabel('Fractional Uptake Difference (Da)', fontname = font, fontsize = axissize)
            for tick in self.ax.get_xticklabels():
                tick.set_fontname(font)
                tick.set_fontsize(16)
            for tick in self.ax.get_yticklabels():
                tick.set_fontname(font)
                tick.set_fontsize(16)
            num_colors = 102
            cm = plt.get_cmap('bwr')
            colorlist = [cm(1. * i / num_colors) for i in range(num_colors)]
            n = 0
            y_val_list = list()
            exp = self.parent.exposure.get()
            exposures = list(self.e_dict[0].keys())
            if exp == 'All':
                for e in sorted(exposures):
                    c = int(1 + n * 100 / (len(exposures) - 1))
                    n = n + 1
                    self.ax.plot(np.array([i[0] for i in self.e_dict[0][e]]), np.array([i[1]-j[1] for i,j in zip(self.e_dict[0][e],self.e_dict[1][e])]), label = str(e)+' min', color=colorlist[c])
                    y_val_list.extend([i[1]-j[1] for i,j in zip(self.e_dict[0][e],self.e_dict[1][e])])
                self.ax.set_ylim((min(y_val_list), max(y_val_list)))
            else:
                e = float(exp)
                self.ax.plot(np.array([i[0] for i in self.e_dict[0][e]]),
                           np.array([i[1] - j[1] for i, j in zip(self.e_dict[0][e], self.e_dict[1][e])]), label = str(e)+' min', color='black')
                y_val_list.extend([i[1] - j[1] for i, j in zip(self.e_dict[0][e], self.e_dict[1][e])])
                self.ax.set_ylim((min(y_val_list), max(y_val_list)))
            self.leg = self.ax.legend()
            self.leg.set_draggable(True)

        # Saves the plot to a file
        def export(self):
            savefile = FileDialog.asksaveasfilename(
                title="Choose save location",
                initialfile=self.protein + "_" + str(self.state1) + "-" + str(self.state2) + "_Map.png",
                initialdir=expanduser('~'),
                filetypes=(("PNG", '*.png'), ("all files", "*.*")))
            self.fig.savefig(savefile, dpi=300, transparent=True)

        # Updates the figure based on changed settings
        def update(self, key, value):
            self.settings[key] = value
            self.fig.clear()
            self.ax = self.fig.subplots()
            self.createFigure()
            plt.draw()

        # Context menu - why is this here?
        class context():
            # Initializes context menu
            def __init__(self, frame, main, *args, **kwargs):
                self.main = main
                fonts = sorted(set([f.name for f in matplotlib.font_manager.fontManager.ttflist]))

                popupMenu = Tk.Menu(frame, tearoff=0)
                fontMenu = Tk.Menu(popupMenu, tearoff=0)
                titleMenu = Tk.Menu(popupMenu, tearoff=0)
                axisMenu = Tk.Menu(popupMenu, tearoff=0)

                popupMenu.add_cascade(label="Font", menu=fontMenu)
                popupMenu.add_cascade(label="Title", menu=titleMenu)
                popupMenu.add_cascade(label="Axis", menu=axisMenu)
                for i in fonts:
                    if i == 'Arial':
                        fontMenu.add_command(label=i, command=lambda i=i: main.update('Font', i), background='green')
                    else:
                        fontMenu.add_command(label=i, command=lambda i=i: main.update('Font', i))
                self.title_val = Tk.IntVar()
                titleMenu.add_command(label="Show/Hide Title", command=lambda: main.update("Show Title", not main.settings['Show Title']))
                self.title_val.set(True)

                stMenu = Tk.Menu(titleMenu, tearoff=0)
                titleMenu.add_cascade(label="Title Size", menu=stMenu)
                for i in list(range(16, 65, 4)):
                    stMenu.add_command(label=i, command=lambda i=i: main.update('Title Size', i))

                axisMenu.add_command(label="Show/Hide X Title", command=lambda: main.update("Show X Title", not main.settings['Show X Title']))

                axisMenu.add_command(label="Show/Hide Y Title", command=lambda: main.update("Show Y Title", not main.settings['Show Y Title']))

                sMenu = Tk.Menu(axisMenu, tearoff=0)
                axisMenu.add_cascade(label="Y Axis Size", menu=sMenu)
                for i in list(range(16, 65, 4)):
                    sMenu.add_command(label=i, command=lambda i=i: main.update('Axis Label Size', i))

                popupMenu.add_command(label='Save', command=lambda: main.export())
                popupMenu.add_command(label='Close', command=lambda: plt.close(main.fig))

                if platform == 'darwin':
                    frame.bind("<Button-2>", lambda event: popupMenu.post(event.x_root, event.y_root))
                else:
                    frame.bind("<Button-3>", lambda event: popupMenu.post(event.x_root, event.y_root))
                self.main = main

class PML():
    '''
    PyMOL Script Generator Window Controller
    Widgets
        Protein, State, Exposure, Color, Range, Representation dropdown menus
    '''
    # Initializes Generator Window, Fills widgets, Defines Tk Variables, Binds Buttons
    def __init__(self, main):
        self.style = ttk.Style()
        self.main = main
        self.view = View.PML(main)
        main.widget_dict['PML'] = self.view.top
        self.main.view.window_menu.add_command(label='PML', command=lambda: self.main.focus(self.view.top))
        self.states = main.statesvar.get()
        self.proteins = main.proteinsvar.get()
        self.exposures = ['Ave']
        self.exposures.extend(main.exposuresvar.get())
        self.representation = Tk.StringVar()
        self.protein = Tk.StringVar()
        self.state1 = Tk.StringVar()
        self.state2 = Tk.StringVar()
        self.exposure = Tk.StringVar()
        self.color = Tk.StringVar()
        self.range = Tk.StringVar()

        self.view.combobox_representation.configure(values =
                                                   ['Uptake', 'SD', 'Uptake Diff',
                                                    'Fractional Uptake', 'Fractional SD (RSD)', 'Fractional Uptake Diff'])
        self.view.combobox_representation.bind("<<ComboboxSelected>>", lambda e: self.selected())
        self.view.combobox_representation.configure(textvariable = self.representation)
        self.view.combobox_protein.configure(values = self.proteins)
        self.view.combobox_protein.configure(textvariable = self.protein)
        self.view.combobox_state1.configure(values = self.states)
        self.view.combobox_state1.configure(textvariable = self.state1)
        self.view.combobox_state2.configure(values = [])
        self.view.combobox_state2.configure(textvariable = self.state2)
        self.view.combobox_state2.set(['None'])
        self.view.combobox_exposure.configure(values = self.exposures)
        self.view.combobox_exposure.configure(textvariable = self.exposure)
        self.view.button_pdb.configure(command = lambda: self.load_pdb(), state='disabled')
        self.view.combobox_color.configure(values = ['red_white_blue', 'blue_white_red', 'rainbow', 'rainbow_rev',])
        self.view.combobox_color.configure(textvariable = self.color)
        self.view.combobox_color.set('blue_white_red')
        self.view.combobox_range.configure(values = ['Fit to Data', 'Full Range', 'Half Range'])
        self.view.combobox_range.configure(textvariable = self.range)
        self.view.combobox_range.set('Half Range')
        self.view.button_save.configure(command = lambda: self.save(self.representation.get()))
        self.view.button_help.configure(command=lambda: self.main.help('7) Pymol Scripts'))
        self.view.top.protocol("WM_DELETE_WINDOW", lambda: main.on_closing(self.view.top, 'PML'))

    # Updates widgets based on representation selection
    def selected(self):
        selection = self.representation.get()
        if selection in ['Uptake Diff','Fractional Uptake Diff']:
            self.view.combobox_state2.configure(values=self.states)
            self.view.combobox_state2['state'] = 'readonly'
            self.view.label_state2.configure(text='''State 2''')
            self.view.combobox_state2.configure(textvariable=self.state2)
        else:
            self.view.combobox_state2.configure(values=[])
            self.view.combobox_state2['state'] = 'disabled'
            self.view.label_state2.configure(text='''''')
            self.view.combobox_state2.set(['None'])

    # Future function to import a pdb for the selection of chain letters, names
    def load_pdb(self):
        print('Stuff')

    # Save Pml script to file
    def save(self,type):
        if type in ['Uptake Diff','Fractional Uptake Diff']:
            name = str(type).replace(' ','')+'_'+str(self.state1.get()) + '-' + str(self.state2.get()) + '_' + str(self.protein.get()) + '_' + str(
                self.exposure.get()) + 'min'
            self.main.state_obj[self.main.csvfile].exportPymol(
                self.main.corrected, self.state1.get(), self.state2.get(), self.protein.get(), self.exposure.get(),
                self.view.entry_protein.get(), self.view.entry_chain.get(), type, self.color.get(), self.range.get(),
                FileDialog.asksaveasfilename(title="Choose save location", defaultextension=".pml", initialfile=name))

        else:
            name = str(type).replace(' ','')+'_'+str(self.state1.get())+'_'+str(self.protein.get())+'_'+str(self.exposure.get())+'min'
            self.main.state_obj[self.main.csvfile].exportPymol(
                self.main.corrected, self.state1.get(), None, self.protein.get(), self.exposure.get(),
                self.view.entry_protein.get(), self.view.entry_chain.get(), type, self.color.get(), self.range.get(),
                FileDialog.asksaveasfilename(title="Choose save location", defaultextension=".pml", initialfile=name))

class RT():
    '''
    Rentention Time Prediction Window Controller
    Widgets:
        Peptide Sequence Entry
        Calculate Button
    '''
    #Initialize window, bind button
    def __init__(self, main):
        self.main = main
        self.rtv = View.RT(main.root)
        self.rtv.button_calculate.configure(command=lambda:self.calculate())

    # Calculate new retention time
    def calculate(self):
        rt = self.main.state_obj[self.main.csvfile].predictRT(str(self.rtv.text.get()))
        self.rtv.label_RT.configure(text="RT: "+str(round(rt,3)))

class Console():
    def __init__(self):
        self.log = []
    def write(self,data):
        self.log.append(data)
    def exit(self):
        with open(environ['HOME']+'/console.txt', 'w') as console:
            for line in self.log:
                console.write(line)
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        sys.exit()



if __name__ == '__main__':
    # if not exists('/Users/Shared/.DECA/'):
    #     makedirs('/Users/Shared/.DECA/')
    # with open('/Users/Shared/.DECA/console.txt', 'w') as console:
    #     stderr = stdout = console
    que = queue.Queue()

    c = Main()

    #stderr = stdout = StdRedirector(c.console_box)
    while True:
        try:
            c.run()
            break
        except UnicodeDecodeError:
            pass

