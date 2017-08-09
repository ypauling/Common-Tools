#! /usr/local/bin/env python

##### This is the GUI portal for CRISPR target sites search #####
##### Author: Bing Yang                                     #####
##### Date: 2017-07-27                                      #####

from __future__ import print_function, division
from tkinter import *
from tkinter import messagebox
from tkinter import ttk
from tkinter.filedialog import askopenfilename
import sys
from crispr_search import *


###Function tv_column_sort(tv, col, reverse):
###    Description: A function for sorting the column in a treeview window
###    Parameters:
###      tv: the id of the treeview window
###      col: the column to be sorted
###      reverse: whether to do reverse sorting
###    Return: None


def tv_column_sort(tv: ttk.Treeview, col: str, reverse: bool) -> None:
    l = [(tv.set(k, col), k) for k in tv.get_children('')]
    try:
        l.sort(key=lambda t: int(t[0]), reverse=reverse)
    except ValueError:
        l.sort(reverse=reverse)

    for index, (val, k) in enumerate(l):
        tv.move(k, '', index)

    tv.heading(col, command=lambda: tv_column_sort(tv, col, not reverse))


### Class crisprGUI:
###    This is the GUI used to receive query string and display target sites
###    information

class crisprGUI:

    def __init__(self, master):

        self.master = master
        self.master.title('CRISPR Target Sites Search Script')
        self.text_read = ''
        self.filename_read = ''

        frame_text_input = Frame(master)
        frame_text_input.pack()

        self.label_ti = Label(frame_text_input,
        text='Enter your sequence here')
        self.label_ti.pack()
        self.text_input = Text(frame_text_input)
        self.text_input.pack()

        self.read_ti = Button(frame_text_input, text='READ',
        command=self.read_text)
        self.read_ti.pack()

        frame_file_input = Frame(master)
        frame_file_input.pack()

        self.label_fi = Label(frame_file_input,
        text='Or choose your input fasta file here')
        self.label_fi.pack(fill=X)

        self.read_fi = Button(frame_file_input, text='Choose File',
        command=self.read_file)
        self.read_fi.pack()
        self.deread_fi = Button(frame_file_input, text='Deselect File',
        command=self.deread)
        self.deread_fi.pack()

        self.label_fname = Label(frame_file_input, text='')
        self.label_fname.pack(side=BOTTOM)

        self.cal_button = Button(master, text='Search Sites',
        command=self.analysis)
        self.cal_button.pack()

        self.quit_button = Button(master,text='QUIT', command=master.quit)
        self.quit_button.pack()

    def read_text(self):
        sequence = self.text_input.get(1.0, END)
        self.text_read = sequence

    def read_file(self):
        file_name = askopenfilename(title='Choose a file', initialdir='~/')
        self.filename_read = file_name
        self.label_fname.config(text=file_name)

    def deread(self):
        self.label_fname.config(text='')

    def analysis(self):
        if self.filename_read is '' and self.text_read is '':
            messagebox.showinfo('Read Error', 'No sequence info provided!')
        else:
            if self.filename_read is not '' and self.text_read is not '':
                messagebox.showinfo('Info Conflict',
                        'Do not provide file and text at the same time!')
            else:
                if self.filename_read is not '':
                    file_reader = open(self.filename_read, 'r')
                    fasta_seq = file_reader.readline().rstrip()
                    self.display_result(fasta_seq)
                else:
                    fasta_seq = self.text_read.rstrip()
                    self.display_result(fasta_seq)

    def display_result(self, seq):

        self.pop_result = Toplevel()
        names = ['target', 'total score', 'strand', 'PAM proximal GC>50%',
                '30%<GC<80%', 'seed TTT', 'repeats', 'base 19 T',
                'base 20 C or T', 'base 20 G', 'base 16 G', 'base 16 C']
        tree = ttk.Treeview(self.pop_result, selectmode='browse')
        tree['columns'] = names[1:]

        targets_to_display = []
        targets = get_targets(seq)

        for i in targets:
            check1 = re.search('^CC.+$', i)
            check2 = re.search('^[ATCG{21}GG$]', i)

            if check1 == None and check2 != None:
                targets_to_display.append('-'+(reverse_complement(i)))
                targets_to_display.append('+'+i)
            elif check1 != None and check2 == None:
                targets_to_display.append('-'+reverse_complement(i))
            else:
                targets_to_display.append('+'+i)

        all_entries = []
        for i in targets_to_display:
            tar = i[1:24]
            total_score = PAM_proximal_GC(tar) + overall_GC(tar) + \
            seed_TTT(tar) + repeats(tar) + base_19_check(tar) + \
            undesirable_base_20(tar) + desirable_base_20(tar) + \
            undesirable_base_16(tar) + desirable_base_16(tar)
            entries = [i[1:24], str(total_score),
                    i[0], str(PAM_proximal_GC(tar)),
                    str(overall_GC(tar)), str(seed_TTT(tar)),
                    str(repeats(tar)), str(base_19_check(tar)),
                    str(undesirable_base_20(tar)), str(desirable_base_20(tar)),
                    str(undesirable_base_16(tar)), str(desirable_base_16(tar))]
            all_entries.append(entries)

        for e in all_entries:
            tree.insert('', 0, text=e[0], values=e[1:])

        tree.pack(side='left')

        for i in names[1:]:
            tree.column(i, width=75, anchor=E)
            tree.heading(i, text=i, anchor=E,
                    command=lambda i_=i: tv_column_sort(tree, i_, False))

        sbar = ttk.Scrollbar(self.pop_result, orient='vertical',
        command=tree.yview)
        sbar.pack(side='right', fill='y')

        tree.configure(yscrollcommand=sbar.set)

root = Tk()
my_gui = crisprGUI(root)
root.mainloop()
root.destroy()
