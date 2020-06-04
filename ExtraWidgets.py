# 2020-06-04 Updated EditableTreeview to fix Back Exchange and Replicate tables
try:
    import Tkinter as tk
    import ttk
except ImportError:
    import tkinter as tk
    import tkinter.ttk as ttk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class AutoScroll(object):
    """Configure the scrollbars for a widget."""

    def __init__(self, master):
        #  Rozen. Added the try-except clauses so that this class
        #  could be used for scrolled entry widget for which vertical
        #  scrolling is not supported. 5/7/14.
        try:
            vsb = ttk.Scrollbar(master, orient='vertical', command=self.yview)
        except:
            pass
        hsb = ttk.Scrollbar(master, orient='horizontal', command=self.xview)
        # self.configure(yscrollcommand=_autoscroll(vsb),
        #    xscrollcommand=_autoscroll(hsb))
        try:
            self.configure(yscrollcommand=self._autoscroll(vsb))
        except:
            pass
        self.configure(xscrollcommand=self._autoscroll(hsb))

        self.grid(column=0, row=0, sticky='nsew')
        try:
            vsb.grid(column=1, row=0, sticky='ns')
        except:
            pass
        hsb.grid(column=0, row=1, sticky='ew')
        master.grid_columnconfigure(0, weight=1)
        master.grid_rowconfigure(0, weight=1)
        # Copy geometry methods of master  (taken from ScrolledText.py)
        methods = list(tk.Pack.__dict__.keys()) + list(tk.Grid.__dict__.keys()) + list(tk.Place.__dict__.keys())
        for meth in methods:
            if meth[0] != '_' and meth not in ('config', 'configure'):
                setattr(self, meth, getattr(master, meth))

    @staticmethod
    def _autoscroll(sbar):
        """Hide and show scrollbar as needed."""

        def wrapped(first, last):
            first, last = float(first), float(last)
            if first <= 0 and last >= 1:
                sbar.grid_remove()
            else:
                sbar.grid()
            sbar.set(first, last)

        return wrapped

    def __str__(self):
        return str(self.master)

def _create_container(func):
    """Creates a ttk Frame with a given master, and use this new frame to
    place the scrollbars and the widget."""

    def wrapped(cls, master, **kw):
        container = ttk.Frame(master)
        return func(cls, container, **kw)

    return wrapped

class ScrolledListBox(AutoScroll, tk.Listbox):
    """A standard Tkinter Text widget with scrollbars that will
    automatically show/hide as needed."""

    @_create_container
    def __init__(self, master, **kw):
        tk.Listbox.__init__(self, master, **kw)
        AutoScroll.__init__(self, master)

class ScrolledTreeView(AutoScroll, ttk.Treeview):
    """A standard ttk Treeview widget with scrollbars that will
    automatically show/hide as needed."""

    @_create_container
    def __init__(self, master, **kw):
        ttk.Treeview.__init__(self, master, **kw)
        AutoScroll.__init__(self, master)

class EditableTreeview(ttk.Treeview):  # Based on code from https://github.com/alejandroautalan/pygubu
    def __init__(self, master=None, **kw):
        ttk.Treeview.__init__(self, master, **kw)

        self._curfocus = None
        self._inplace_widgets = {}
        self._inplace_widgets_show = {}
        self._inplace_vars = {}

        self.bind('<<TreeviewSelect>>', self.__check_focus)
        self.bind('<ButtonRelease-1>', self.__check_focus)
        self.bind('<Configure>', lambda e: self.after_idle(self.__update_wnds))
        self.bind('<Double-1>', lambda e: self.inplace_entry(self.__get_display_columns()[-2], self.focus()))

    def return_event(self, col, item):
        self.__update_value(col, item)
        self.__clear_inplace_widgets()

    def __check_focus(self, event):
        """Checks if the focus has changed"""
        changed = False
        if not self._curfocus:
            changed = True
        elif self._curfocus != self.focus():
            self.__clear_inplace_widgets()
            changed = True
        newfocus = self.focus()
        if changed:
            if newfocus:
                self._curfocus = newfocus
                self.__focus(newfocus)
            self.__update_wnds()

    def __focus(self, item):
        """Called when focus item has changed"""
        cols = self.__get_display_columns()
        for col in cols:
            self.__event_info = (col, item)
            if col in self._inplace_widgets:
                w = self._inplace_widgets[col]

    def __update_wnds(self, event=None):
        if not self._curfocus:
            return
        item = self._curfocus
        cols = self.__get_display_columns()
        for col in cols:
            if col in self._inplace_widgets:
                wnd = self._inplace_widgets[col]
                bbox = ''
                if self.exists(item):
                    bbox = self.bbox(item, column=col)
                if bbox == '':
                    wnd.place_forget()
                elif col in self._inplace_widgets_show:
                    wnd.place(x=bbox[0], y=bbox[1],
                              width=bbox[2], height=bbox[3])

    def __clear_inplace_widgets(self):
        """Remove all inplace edit widgets."""
        cols = self.__get_display_columns()
        for c in cols:
            if c in self._inplace_widgets:
                widget = self._inplace_widgets[c]
                widget.place_forget()
                self._inplace_widgets_show.pop(c, None)
                # widget.destroy()
                # del self._inplace_widgets[c]

    def __get_display_columns(self):
        cols = self.cget('displaycolumns')
        show = (str(s) for s in self.cget('show'))
        if '#all' in cols:
            cols = self.cget('columns') + ('#0',)
        elif 'tree' in show:
            cols = cols + ('#0',)
        return cols

    def __get_value(self, col, item):
        if col == '#0':
            return self.item(item, 'text')
        else:
            return self.set(item, col)

    def __set_value(self, col, item, value):
        if col == '#0':
            self.item(item, text=value)
        else:
            self.set(item, col, value)
        self.__event_info = (col, item)
        self.event_generate('<<TreeviewCellEdited>>')

    def __update_value(self, col, item):
        if not self.exists(item):
            return
        value = self.__get_value(col, item)
        newvalue = self._inplace_vars[col].get()
        if value != newvalue:
            self.__set_value(col, item, newvalue)

    def list_dict(self):
        self.event_generate('<Unmap>')
        treedict = {}
        cols = self.__get_display_columns()
        for i in self.get_children():
            # 2020-06-04 Updated to fix Back Exchange and Replicate tables
            treedict[str('-'.join([str(i) for i in self.item(i)['values'][0:-1]]))] = float(self.item(i)['values'][-1])
        return treedict

    def inplace_entry(self, col, item):
        if col not in self._inplace_vars:
            self._inplace_vars[col] = tk.StringVar()
        svar = self._inplace_vars[col]
        svar.set(self.__get_value(col, item))
        if col not in self._inplace_widgets:
            self._inplace_widgets[col] = ttk.Entry(self, textvariable=svar)
        entry = self._inplace_widgets[col]
        entry.bind('<Unmap>', lambda e: self.return_event(col, item))
        entry.bind('<FocusOut>', lambda e: self.return_event(col, item))
        entry.bind("<Return>", lambda e: self.return_event(col, item))
        self._inplace_widgets_show[col] = True
        self.__update_wnds()

class ComboboxTreeview(ttk.Treeview):  # Based on code from https://github.com/alejandroautalan/pygubu
    def __init__(self, master=None, combovals=[], **kw):
        ttk.Treeview.__init__(self, master, **kw)

        self._curfocus = None
        self._inplace_widgets = {}
        self._inplace_widgets_show = {}
        self._inplace_vars = {}

        self.bind('<<TreeviewSelect>>', self.__check_focus)
        self.bind('<ButtonRelease-1>', self.__check_focus)
        self.bind('<Configure>', lambda e: self.after_idle(self.__update_wnds))
        self.bind('<Double-1>', lambda e: self.inplace_entry(self.__get_display_columns()[1], self.focus()))
        self.combovals = combovals

    def return_event(self, col, item):
        self.__update_value(col, item)
        self.__clear_inplace_widgets()

    def __check_focus(self, event):
        """Checks if the focus has changed"""
        changed = False
        if not self._curfocus:
            changed = True
        elif self._curfocus != self.focus():
            self.__clear_inplace_widgets()
            changed = True
        newfocus = self.focus()
        if changed:
            if newfocus:
                self._curfocus = newfocus
                self.__focus(newfocus)
            self.__update_wnds()

    def __focus(self, item):
        """Called when focus item has changed"""
        cols = self.__get_display_columns()
        for col in cols:
            self.__event_info = (col, item)
            if col in self._inplace_widgets:
                w = self._inplace_widgets[col]

    def __update_wnds(self, event=None):
        if not self._curfocus:
            return
        item = self._curfocus
        cols = self.__get_display_columns()
        for col in cols:
            if col in self._inplace_widgets:
                wnd = self._inplace_widgets[col]
                bbox = ''
                if self.exists(item):
                    bbox = self.bbox(item, column=col)
                if bbox == '':
                    wnd.place_forget()
                elif col in self._inplace_widgets_show:
                    wnd.place(x=bbox[0], y=bbox[1],
                              width=bbox[2], height=bbox[3])

    def __clear_inplace_widgets(self):
        """Remove all inplace edit widgets."""
        cols = self.__get_display_columns()
        for c in cols:
            if c in self._inplace_widgets:
                widget = self._inplace_widgets[c]
                widget.place_forget()
                self._inplace_widgets_show.pop(c, None)
                # widget.destroy()
                # del self._inplace_widgets[c]

    def __get_display_columns(self):
        cols = self.cget('displaycolumns')
        show = (str(s) for s in self.cget('show'))
        if '#all' in cols:
            cols = self.cget('columns') + ('#0',)
        elif 'tree' in show:
            cols = cols + ('#0',)
        return cols

    def __get_value(self, col, item):
        if col == '#0':
            return self.item(item, 'text')
        else:
            return self.set(item, col)

    def __set_value(self, col, item, value):
        if col == '#0':
            self.item(item, text=value)
        else:
            self.set(item, col, value)
        self.__event_info = (col, item)
        self.event_generate('<<TreeviewCellEdited>>')

    def __update_value(self, col, item):
        if not self.exists(item):
            return
        elif item == '':
            return
        value = self.__get_value(col, item)
        newvalue = self._inplace_vars[col].get()
        if value != newvalue:
            self.__set_value(col, item, newvalue)

    def treeview_dict(self):
        self.return_event(self.__get_display_columns()[1], self.focus())
        treedict = {}
        for i in self.get_children():
            treedict[self.item(i)['values'][0]] = self.item(i)['values'][1]
        return treedict

    def inplace_entry(self, col, item):
        if col not in self._inplace_vars:
            self._inplace_vars[col] = tk.StringVar()
        svar = self._inplace_vars[col]
        svar.set(self.__get_value(col, item))
        if col not in self._inplace_widgets:
            self._inplace_widgets[col] = ttk.Combobox(self, textvariable=svar, values=self.combovals)

        entry = self._inplace_widgets[col]
        entry.bind('<<ComboboxSelected>>', lambda e: self.return_event(col, item))
        #entry.bind('<Unmap>', lambda e: self.return_event(col, item))
        #entry.bind('<FocusOut>', lambda e: self.return_event(col, item))
        #entry.bind("<Return>", lambda e: self.return_event(col, item))
        self._inplace_widgets_show[col] = True
        self.__update_wnds()

def addscrollingfigure(figure, frame):
    def _on_mousewheel(event):
        canvas.yview_scroll(event.delta, "units")
    for child in frame.winfo_children():
        child.destroy()
    canvas = tk.Canvas(frame)
    canvas.place(relx=0, rely=0, relheight=1, relwidth=1)
    x_scrollbar = tk.Scrollbar(frame, orient=tk.HORIZONTAL)
    y_scrollbar = tk.Scrollbar(frame)
    y_scrollbar.place(relx=0.98, rely=0.1, relheight=0.8)
    x_scrollbar.place(relx=0.1, rely=0.98, relwidth=0.8)
    canvas.config(xscrollcommand=x_scrollbar.set)
    x_scrollbar.config(command=canvas.xview)
    canvas.config(yscrollcommand=y_scrollbar.set)
    canvas.bind_all("<MouseWheel>", _on_mousewheel)
    y_scrollbar.config(command=canvas.yview)
    fig_agg = FigureCanvasTkAgg(figure, canvas)
    mpl_canvas = fig_agg.get_tk_widget()
    mpl_canvas.place(relx=0.5, rely=0.5, relheight=1, relwidth=1, anchor=tk.CENTER)
    cwid = canvas.create_window(1, 1, window=mpl_canvas, anchor=tk.CENTER)
    figure.set_size_inches(figure.get_size_inches())
    wi, hi = [i * figure.dpi for i in figure.get_size_inches()]
    hi = hi + 50
    mpl_canvas.config(width=wi, height=hi, bd=0, highlightthickness=0)
    canvas.itemconfigure(cwid, width=wi, height=hi)
    canvas.config(scrollregion=canvas.bbox(tk.ALL), width=wi, height=hi, bd=0, highlightthickness=0, )
    canvas.yview_moveto(.0)
    figure.canvas.draw()
    return fig_agg

def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total:
        print()