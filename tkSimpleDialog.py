#!/usr/bin/env python
# This class builds all the windows for the GUI
# Each window has its own class with different layouts


from Tkinter import *
import ttk
import os
import tkFileDialog
import tkMessageBox
import ToolTip


#colors
back_color = "#aa00ff"
emerald ='#008a00'
UI_blue='#2d89ef'
orange ='#da532c'#e8980e'#'#cd3700'
cobalt =UI_blue#'#0050ef'
steel ="#1d1d1d"#647687'
green ='#60a917'
mango = '#f09609'
azzurro = '#1ba1f2'
text_color='white'#'#ffffff'
txt_dim =15
global infolder_
global read1_
global read2_
predict_header = """
Automatic fastq analysis for sample-wise variant calling.
Each sample is analyzed and an output vcf file is produced.

Parameters are needed to locate the input folder and the output folder,
as well as the required computing resources.
"""


explain_attempt= "Namestart and Namelength define how to generate subfolders.\nFor each input-file a subfolder is produced from the sample name, startng from 'namestart' position for 'namelength' characters"
explain_attempt+="""

For example:
    Sample_name = xx_1234_xx.read1.fastq.gz
    Namestart = 4
    Namelength = 4
    Output Subolder = 1234

All samples with the same Subfolder will be grouped.
        """
prioritize_header="""
Variant annotation, filtering and prioritzation for trios
or families analyzed jointly.

Please deifne an output folder and the inheritance modes
for the filtering phase to run the tool.
"""


font_type ="Helvetica"

class Dialog(Toplevel):
    '''
    Dialog superclass on which th eannotate Predict and Prioritize windows will be modeled
    This dialog generates a new window, disables the original window buttons and move the focus
    to the newly created window.
    '''

    def __init__(self, parent, title = None):

        Toplevel.__init__(self, parent)
        self.transient(parent)

        if title:
            self.title(title)

        self.parent = parent

        self.result = None

        body = Frame(self,background = steel)
        self.initial_focus = self.body(body)
        body.pack(padx=5, pady=5)

        self.buttonbox()

        #self.grab_set()

        if not self.initial_focus:
            self.initial_focus = self

        self.protocol("WM_DELETE_WINDOW", self.cancel)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,
                                  parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        self.wait_window(self)

    def body(self, master):
        # create dialog body.  return widget that should have
        # initial focus.  this method should be overridden
        pass

    def buttonbox(self):
        # add standard button box. override if you don't want the
        # standard buttons
        box = Frame(self,background = steel)
        w = Button(box, text="Run", width=10, command=self.ok, default=ACTIVE)
        w.pack(side=LEFT, padx=5, pady=5)
        w = Button(box, text="Back", width=10, command=self.cancel)
        w.pack(side=LEFT, padx=5, pady=5)
        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)
        box.pack()

    def ok(self, event=None):

        if not self.validate():
            self.initial_focus.focus_set() # put focus back
            return'whit'
        self.withdraw()
        self.update_idletasks()
        self.apply()
        self.cancel()
        return self.result

    def cancel(self, event=None):

        # put focus back to the parent window
        self.parent.focus_set()
        self.destroy()
        return None

    def validate(self):

        return 1 # override

    def apply(self):

        pass # override


class MyDialog(Dialog):
    '''
    Just a test class : to be removed
    '''
    def file_search(self):
        self.file = tkFileDialog.askopenfile(parent=self,mode='rb',filetypes=[('Input file','*.vcf'), ('All files','*')],title='Choose a subtitle file')
        self.file=abs_path = os.path.abspath(self.file.name)
        #return abs_path

    def body(self, master):
        Label(master, text="First:",background=steel).grid(row=0,columnspan=1)
        Label(master, text="Second:",background=steel).grid(row=1,columnspan=2)
        self.e1 = Entry(master)
        self.e2 = Entry(master)
        self.e1.grid(row=0, column=1)
        self.e2.grid(row=1, column=1)
        var = IntVar()
        self.c = Checkbutton(master, text="Expand", variable=var)
        self.c.grid(row=2,column=0)
        self.listbox = Listbox(master,selectmode=MULTIPLE,height=5)
        self.listbox.insert(END, "a list entry")
        for item in ["one", "two", "three", "four"]:
            self.listbox.insert(END, item)
        self.listbox.grid(row=3,column=1)

        #Combo box button
        self.box_value = StringVar()
        self.box =ttk.Combobox(master, state='readonly',textvariable=self.box_value,background=steel)
        self.box['values'] = ('X', 'Y', 'Z')
        self.box.current(0)
        self.box.grid(row=3,column=0)

        #File selection button
        self.filename = Button(master,text='Input File',command=self.file_search,width=50)
        self.filename.grid(row=4,column=0,columnspan=2)

        return self.e1 # initial focus
    def apply(self):
        first = (self.e1.get())
        second = (self.e2.get())
        combo = self.box.get()
        fname = self.file
        self.result = first, second,combo,fname



class Annotate_window(Dialog):
    '''
    Window to gather the needed Annotate.py values : the output is a dictionary with variables to launch the command.
    '''
    def file_search(self):
        self.file = tkFileDialog.askopenfile(parent=self,mode='rb',filetypes=[('Input file','*.vcf'), ('All files','*')],title='Choose the input annotation file')
        if self.file != None:
            self.file=abs_path = os.path.abspath(self.file.name)
        #return abs_path

    def body(self, master):
        self.configure(background =steel)
        self.resizable(0,0)
        #Title
        Label(master, text="\nAnnotate Options:\n",background=orange,fg=text_color,font=(font_type, txt_dim)).grid(row=0,column=0,columnspan=3,sticky='W'+'E')


        #File selection button
        self.filename = Button(master,text='Input File',command=self.file_search,width=50,fg=text_color,font=(font_type, txt_dim)
                ,highlightbackground= steel,background=  steel,activebackground= orange,activeforeground=text_color        )
        self.filename.grid(row=1,column=0,columnspan=3)

        #Select Variant Type
        Label(master, text="Select the variant type:",background=steel,fg=text_color,font=(font_type, txt_dim)).grid(row=2,column=0,columnspan=1,sticky=W+E)
        #Combo box button
        self.var_value = StringVar()
        self.var =ttk.Combobox(master,state='readonly', textvariable=self.var_value)
        self.var.configure(background=steel)
        self.var['values'] = ('snp', 'indel', 'all')
        self.var.current(2)
        self.var.grid(row=2,column=1,columnspan=2,sticky=W+E)

        #SelectGeneMode ensGene,refGene,knownGene,all
        Label(master, text="Select the Gene Definition type:",background=steel,fg=text_color,font=(font_type, txt_dim)).grid(row=3,column=0,columnspan=1,sticky=W)
        #Combo box button
        self.gdef_value = StringVar()
        self.gdef =ttk.Combobox(master,state='readonly', textvariable=self.gdef_value,background=steel)
        self.gdef['values'] = ('ensGene','refGene','knownGene','all')
        self.gdef.current(2)
        self.gdef.grid(row=3,column=2,columnspan=2)

        #SelectGeneDefinition
        Label(master, text="Select the Genotype Mode:",background=steel,fg=text_color,font=(font_type, txt_dim)).grid(row=4,column=0,columnspan=1,sticky=W)
        #Combo box button
        self.gtype_value = StringVar()
        self.gtype =ttk.Combobox(master,state='readonly', textvariable=self.gtype_value,background=steel)
        self.gtype['values'] = ('compact','complete')
        self.gtype.current(0)
        self.gtype.grid(row=4,column=2,columnspan=2)

        #Force and Multisample tickmarks:
        self.var_f = IntVar()
        self.f = Checkbutton(master, text="Force writing?", variable=self.var_f,background=steel,fg=text_color,selectcolor=steel
                             ,highlightthickness=0,font=(font_type, txt_dim),activebackground= orange ,activeforeground=text_color    )
        self.f.grid(row=5,column=0,columnspan=1,sticky=W)

        self.var_o = IntVar()
        self.o = Checkbutton(master, text="onlyGenicAnnotation?", variable=self.var_o,background=steel,fg=text_color,selectcolor=steel,
                             highlightthickness=0,font=(font_type, txt_dim),activebackground= orange,activeforeground=text_color     )
        self.o.grid(row=5,column=2,columnspan=1,sticky=E)

        return self.filename # initial focus
    def apply(self):
        force       = (self.var_f.get())
        only_gen    = (self.var_o.get())
        vartype     = self.var.get()
        gdef        = self.gdef.get()
        gtype       = self.gtype.get()
        infile      = self.file
        self.result = {"force":force,"only_gen":only_gen,"vartype":vartype,"gdef":gdef,"infile":infile,"gtype":gtype}

    def validate(self):
        try:
            infile = self.file
            if os.path.isfile(infile):
                return 1
            else:
                raise
        except:
            tkMessageBox.showwarning("Annotate","Please Select one Input File")
        return 0

    def buttonbox(self):
        # add standard button box. override if you don't want the
        # standard buttons
        box = Frame(self,background = steel)
        w = Button(box, text="Run", width=10, command=self.ok,fg=text_color,font=(font_type, txt_dim)
                ,highlightbackground= steel,background=  steel,activebackground= orange,activeforeground=text_color         )
        w.pack(side=LEFT, padx=5, pady=5)
        w = Button(box, text="Back", width=10, command=self.cancel,fg=text_color,font=(font_type, txt_dim)
                ,highlightbackground= steel,background=  steel,activebackground= orange,activeforeground=text_color         )
        w.pack(side=LEFT, padx=5, pady=5)
        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)
        box.pack()
class Predict_window(Dialog):
    '''
    Window to gather the needed Annotate.py values : the output is a dictionary with variables to launch the command.
    '''
    def file_search(self):
        str_default='Config. File'
        self.file = tkFileDialog.askopenfile(parent=self,mode='rb',filetypes=[('All files','*')],title='Choose the  configuration file')
        if self.file != None:
            self.file =abs_path = os.path.abspath(self.file.name)
            self.filename_string.set(str_default + ' [...%s]'%self.file[-15:])
        else:
            self.file = None

        #return abs_path

    def infolder_search(self):
        str_default='Input Folder'
        try:
            self.infolder = tkFileDialog.askdirectory(parent=self)
            if self.infolder != None and os.path.isdir(self.infolder):
                self.infolder=os.path.abspath(self.infolder)
                self.infolder_string.set(str_default + ' [...%s]'%self.infolder[-15:])
                #self.infolder_string.set('In_dir : %s'%self.infolder)
        except:
            self.infolder=None


    def outfolder_search(self):
        str_default='Output Folder'
        try:
            self.outfolder = tkFileDialog.askdirectory(parent=self)
            if self.outfolder != None and os.path.isdir(self.outfolder):
                self.outfolder=os.path.abspath(self.outfolder)
                self.outfolder_string.set(self.outfolder_string.get() + ' [...%s]'%self.outfolder[-15:])
            else:
                self.outfolder=None
        except:
            pass

    def plus1(self,var,l):
        var.set(int(var.get())+1)
        l.textvariable = var
        return var.get()

    def minus1(self,var,l):
        var.set(max(1,int(var.get())-1))
        l.textvariable = var
        return var.get()



    def body(self, master):
        #Title
        cur_row= 0
        self.configure(background =steel)
        self.resizable(0,0)
        Label(master, text="eDiVa Variant Prediction:",background=emerald,fg=text_color,font=(font_type, txt_dim)).grid(row=cur_row,column=0,columnspan=3,sticky='W'+'E')
        cur_row+=1
        Label(master, text=predict_header,background=emerald,fg=text_color,font=(font_type, txt_dim-4),justify=LEFT).grid(row=cur_row,column=0,columnspan=3,sticky='W'+'E')
        cur_row+=1
        Label(master, text="\n",background=steel,fg=text_color,font=(font_type, txt_dim-6)).grid(row=cur_row,column=0,columnspan=3,sticky='W'+'E')
        cur_row+=1
        #Infolder selection button
        self.infolder_string = StringVar()
        self.infolder_string.set('Input Folder')
        self.infolder_ = Button(master,textvariable=self.infolder_string,command=self.infolder_search,width=50,background=steel,fg=text_color,
                                font=(font_type, txt_dim),activebackground= emerald,activeforeground=text_color  )
        #self.infolder_ = Button(master,text='Input directory',command=self.infolder_search,width=50,background=steel,fg=text_color,
        #                font=(font_type, txt_dim),activebackground= emerald,activeforeground=text_color  )
        self.infolder_.grid(row=cur_row,column=0,columnspan=3)
        t1 = ToolTip.ToolTip(self.infolder_, follow_mouse=1, text="Select here where the fastq files you want to analyze are located.")

        #Outfolder selection button
        cur_row+=1
        self.outfolder_string = StringVar()
        self.outfolder_string.set('Output Folder')
        self.outfolder_ = Button(master,textvariable=self.outfolder_string,command=self.outfolder_search,width=50,
                                 background=steel,fg=text_color,font=(font_type, txt_dim),activebackground= emerald,activeforeground=text_color  )
        self.outfolder_.grid(row=cur_row,column=0,columnspan=3)
        t2 = ToolTip.ToolTip(self.outfolder_, follow_mouse=1, text="Select here where you want the output to be stored.")

        #Config file selection button
        cur_row+=1
        self.filename_string = StringVar()
        self.filename_string.set('Config. File')
        self.filename = Button(master,textvariable=self.filename_string,command=self.file_search,width=50,background=steel,fg=text_color,
                               font=(font_type, txt_dim),activebackground= emerald,activeforeground=text_color  )
        self.filename.grid(row=cur_row,column=0,columnspan=3)
        t3 = ToolTip.ToolTip(self.filename, follow_mouse=1, text="Select the global configuration file is stored. It contains the paths to all tools used by eDiVa.")



        #Namestart label and buttons
        cur_row+=1
        self.var = IntVar()
        self.var.set(1)
        self.l = Label(master, textvariable = (self.var),width=20, background=steel,fg=text_color,font=(font_type, txt_dim))
        self.nl = Label(master, text = 'Namestart',width=20, background=steel,fg=text_color,font=(font_type, txt_dim-4))
        self.l.grid(row=cur_row+1,column=0,columnspan=1)#,sticky='E')
        self.nl.grid(row=cur_row,column=0,columnspan=1)#,sticky='W')
        self.plus_n = Button(master,text = '+',command=lambda : self.var.set(self.plus1(self.var,self.l) ) ,
                                 background=steel,fg=text_color,font=(font_type, txt_dim),activebackground= emerald,activeforeground=text_color  )
        self.minus_n = Button(master,text= '-',command=lambda : self.var.set(self.minus1(self.var,self.l) ),
                                 background=steel,fg=text_color,font=(font_type, txt_dim),activebackground= emerald,activeforeground=text_color  )
        self.plus_n.grid(row=cur_row+1,column=0,sticky='E')
        self.minus_n.grid(row=cur_row+1,column=0,sticky='W')
        Label(master, text="",width=10, background=steel,fg=text_color,font=(font_type, txt_dim-4)).grid(row=cur_row,column=1)
        #Namelength label and buttons
        self.var_length = IntVar()
        self.var_length.set(1)
        self.length = Label(master, textvariable = (self.var_length),width=20, background=steel,fg=text_color,font=(font_type, txt_dim))
        self.nlength = Label(master, text = 'Namelength',width=20, background=steel,fg=text_color,font=(font_type, txt_dim-4))
        self.length.grid(row=cur_row+1,column=2,columnspan=1)
        self.nlength.grid(row=cur_row,column=2,columnspan=1)
        self.plus_n_length = Button(master,text ='+',command=lambda : self.var_length.set(self.plus1(self.var_length,self.length) ),
                                 background=steel,fg=text_color,font=(font_type, txt_dim),activebackground= emerald,activeforeground=text_color  )
        self.minus_n_length = Button(master,text='-',command=lambda : self.var_length.set(self.minus1(self.var_length,self.length) ),
                                 background=steel,fg=text_color,font=(font_type, txt_dim),activebackground= emerald,activeforeground=text_color  )
        self.plus_n_length.grid(row=cur_row+1,column=2,sticky='E')
        self.minus_n_length.grid(row=cur_row+1,column=2,sticky='W')

        t2 = ToolTip.ToolTip(self.nl, follow_mouse=1, text=explain_attempt, wraplength=500)
        cur_row+=2


        ##Interactive example of name start:namelength
        #self.example = StringVar()
        #self.example_label= Label(master,textvariable= self.example,width=20, background=steel,fg=text_color,font=(font_type, txt_dim))
        #cur_row+=1

        #CPU label and buttons
        self.cpu = IntVar()
        self.cpu.set(4)
        self.lcpu = Label(master, textvariable = (self.cpu),width=20, background=steel,fg=text_color,font=(font_type, txt_dim))
        self.nlcpu = Label(master, text = 'CPU',width=20, background=steel,fg=text_color,font=(font_type, txt_dim-4))
        self.lcpu.grid(row=cur_row+1,column=0,columnspan=1)#,sticky='E')
        self.nlcpu.grid(row=cur_row,column=0,columnspan=1)#,sticky='W')
        self.plus_cpu = Button(master,text = '+',command=lambda : self.cpu.set(self.plus1(self.cpu,self.lcpu) ) ,
                                 background=steel,fg=text_color,font=(font_type, txt_dim),activebackground= emerald,activeforeground=text_color  )
        self.minus_cpu = Button(master,text= '-',command=lambda : self.cpu.set(self.minus1(self.cpu,self.lcpu) ),
                                 background=steel,fg=text_color,font=(font_type, txt_dim),activebackground= emerald,activeforeground=text_color  )
        self.plus_cpu.grid(row=cur_row+1,column=0,sticky='E')
        self.minus_cpu.grid(row=cur_row+1,column=0,sticky='W')

        #MEM GB label and buttons
        self.mem = IntVar()
        self.mem.set(12)
        self.l_mem = Label(master, textvariable = (self.mem),width=20, background=steel,fg=text_color,font=(font_type, txt_dim))
        self.n_mem = Label(master, text = 'GB Memory',width=20, background=steel,fg=text_color,font=(font_type, txt_dim-4))
        self.l_mem.grid(row=cur_row+1,column=2,columnspan=1)
        self.n_mem.grid(row=cur_row,column=2,columnspan=1)
        self.plus_mem = Button(master,text ='+',command=lambda : self.mem.set(self.plus1(self.mem,self.l_mem) ),
                                 background=steel,fg=text_color,font=(font_type, txt_dim),activebackground= emerald,activeforeground=text_color  )
        self.minus_mem = Button(master,text='-',command=lambda : self.mem.set(self.minus1(self.mem,self.l_mem) ),
                                 background=steel,fg=text_color,font=(font_type, txt_dim),activebackground= emerald,activeforeground=text_color  )
        self.plus_mem.grid(row=cur_row+1,column=2,sticky='E')
        self.minus_mem.grid(row=cur_row+1,column=2,sticky='W')

        #Select Indel Caller
        cur_row+=2
        Label(master, text="Select the Indel caller type:",background=steel,fg=text_color,font=(font_type, txt_dim)).grid(row=cur_row,column=0,columnspan=2,sticky=W)
        #Combo box button
        self.var_value = StringVar()
        self.var_indel =ttk.Combobox(master,state='readonly', textvariable=self.var_value, background=emerald,font=(font_type, txt_dim-4))
        self.var_indel['values'] = ('gatk', 'clindel', 'both')
        self.var_indel.current(0)
        self.var_indel.grid(row=cur_row,column=2,sticky=W+E)
        t2 = ToolTip.ToolTip(self.var_indel, follow_mouse=1, text="Choose the INDEL caller type: GATK or Clindel")

        #Force tick.
        cur_row+=1
        self.var_f = IntVar()
        self.var_f.set(1)
        Label(master, text="Do you want to fuse SNP and INDEL?:", background=steel,fg=text_color,font=(font_type, txt_dim)
              ).grid(row=cur_row,column=0,columnspan=2,sticky=W)
        self.f = Checkbutton(master, text="", variable=self.var_f,background=steel,fg=text_color,selectcolor=steel
                             ,highlightthickness=0,font=(font_type, txt_dim),activebackground= emerald ,activeforeground=text_color )
        self.f.grid(row=cur_row,column=2,columnspan=2)

        #Read extensions
        cur_row+=1
        self.read1 = Entry(master, background=emerald,fg=text_color,font=(font_type, txt_dim-2))
        self.read2 = Entry(master, background=emerald,fg=text_color,font=(font_type, txt_dim-2))
        self.read1.insert(0,'read1.fastq.gz')
        self.read2.insert(0,'read2.fastq.gz')
        Label(master, text="1st read extension:", background=steel,fg=text_color,font=(font_type, txt_dim)).grid(row=cur_row,column=0,columnspan=1,sticky=W)
        Label(master, text="2nd read extension:", background=steel,fg=text_color,font=(font_type, txt_dim)).grid(row=cur_row+1,column=0,columnspan=1,sticky=W)
        self.read1.grid(row=cur_row, column=1,columnspan=2,sticky=W+E)
        cur_row+=1
        self.read2.grid(row=cur_row, column=1,columnspan=2,sticky=W+E)

        #qsubname and max_coverage
        cur_row+=1
        self.qsub = Entry(master, background=emerald,fg=text_color,font=(font_type, txt_dim-2))
        self.qsub.insert(0,'predict.sh')
        self.max_cov= Entry(master, background=emerald,fg=text_color,font=(font_type, txt_dim-2))
        self.max_cov.insert(0,'400')
        Label(master, text="Qsub name:", background=steel,fg=text_color,font=(font_type, txt_dim)).grid(row=cur_row,column=0,columnspan=1,sticky=W)
        Label(master, text="maximum coverage:", background=steel,fg=text_color,font=(font_type, txt_dim)).grid(row=cur_row+1,column=0,columnspan=1,sticky=W)
        self.qsub.grid(row=cur_row, column=1,columnspan=2,sticky=W+E)
        cur_row+=1
        self.max_cov.grid(row=cur_row, column=1,columnspan=2,sticky=W+E)

        return self.infolder_ # initial focus
    def buttonbox(self):
        # add standard button box. override if you don't want the
        # standard buttons
        box = Frame(self,background = steel)
        w = Button(box, text="Run", width=10, command=self.ok,fg=text_color,font=(font_type, txt_dim)
                ,highlightbackground= steel,background=  steel,activebackground= emerald ,activeforeground=text_color        )
        w.pack(side=LEFT, padx=5, pady=5)
        w = Button(box, text="Back", width=10, command=self.cancel,fg=text_color,font=(font_type, txt_dim)
                ,highlightbackground= steel,background=  steel,activebackground= emerald,activeforeground=text_color       )
        w.pack(side=LEFT, padx=5, pady=5)
        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)
        box.pack()

    def ok(self, event=None):
        if not self.validate():
            self.initial_focus.focus_set() # put focus back
            return'whit'
        #
        self.update_idletasks()
        self.apply()
        if self.result != None:
            self.cancel()
            #self.withdraw()
            return self.result
        else:
            return None

    def apply(self):
        NewWin = qSubDetails(self.master)
        cpu_mem_params = ',-pe smp %d,-l virtual_free=%dG,-N %s ,-e %s/, -o %s/'%(self.cpu.get(),self.mem.get(),self.qsub.get(),self.outfolder,self.outfolder)
        if NewWin.result != None:
            self.result = {
                "force"     :self.var_f.get(),
                "infolder"  :self.infolder,
                "outfolder" :self.outfolder,
                "config"    :self.file,
                "namestart" :int(self.var.get()),
                "namelength":int(self.var_length.get()),
                "read1"     :self.read1.get(),
                "read2"     :self.read2.get(),
                "qsubname"  :self.qsub.get(),
                "max_cov"   :self.max_cov.get(),
                "indel"     :self.var_indel.get(),
                "mem"       :self.mem.get(),
                "cpu"       :self.cpu.get(),
                "qopts"     :NewWin.result+cpu_mem_params

                }
                        #IDEA: setup master.something with the fields I want and thats it

            infolder_ =self.infolder
            read1_    =self.read1.get()
            read2_    = self.read2.get()
            first_read  = [each for each in os.listdir(infolder_) if each.endswith(read1_)]
            second_read = [each for each in os.listdir(infolder_) if each.endswith(read2_)]
            for i in range(0,len(first_read)):
                tmp = first_read[i]
                tmp = tmp[:-len(read1_)]
                first_read[i] = tmp
            for i in range(0,len(second_read)):
                tmp = second_read[i]
                tmp = tmp[:-len(read2_)]
                second_read[i] = tmp
            self.master.sample_list = list(set(first_read) & set(second_read))

            get_samples = Sample_list_window(self.master)
            if get_samples.result != None:
                self.result['sample_list']  = get_samples.result
            else:
                self.result = None
                #print self.result
        else :
            self.result = None




    def validate(self):
        #validate infolder
        try:
            dummy = self.infolder
            if os.path.isdir(dummy):
                pass
            else:
                raise
        except:
            tkMessageBox.showwarning("Predict","Please Select one Input Folder")
            return 0
        #validate outfolder
        try:
            dummy = self.outfolder
            if os.path.isdir(dummy):
                pass
            else:
                raise
        except:
            tkMessageBox.showwarning("Predict","Please Select one Output Folder")
            return 0
        #valiate file
        try:
            infile = self.file
            if not(os.path.isfile(infile)):
                raise
        except:
            tkMessageBox.showwarning("Predict","Please Select one Input File")
            return 0
        ##check that namestart < namelength
        #if int(self.var_length.get())-int(self.var.get())>=0:
        #    pass
        #else:
        #    tkMessageBox.showwarning("Predict","Namelength must be at least than namestart")
        #    return 0
        #check that exist the first and second read ?
        if len(self.read1.get())<1 or len(self.read2.get())<1:
            tkMessageBox.showwarning("Predict","Read extensions must contain text")
            return 0
        #check qsubname
        if len(self.qsub.get())<1 :
            tkMessageBox.showwarning("Predict","Qsub must contain text")
        elif self.qsub.get()[0].isdigit():
            tkMessageBox.showwarning("Predict","Qsub must must not start with a digit")
            return 0
        #check maximum coverage exists
        if len(self.max_cov.get())<1 :
            tkMessageBox.showwarning("Predict","Maximum coverage must contain text")
            return 0
        else:
            try :
                dummy = int(self.max_cov.get())
            except:
                tkMessageBox.showwarning("Predict","Maximum coverage must be an integer")
                return 0
        return 1
class Prioritize_window(Dialog):
    '''
    Window to gather the needed Annotate.py values : the output is a dictionary with variables to launch the command.
    '''
    def plus1(self,var,l):
        var.set(int(var.get())+1)
        l.textvariable = var
        return var.get()

    def minus1(self,var,l):
        var.set(max(1,int(var.get())-1))
        l.textvariable = var
        return var.get()

    def config_file_search(self):
        str_default='Global configuration file'
        file_ = tkFileDialog.askopenfile(parent=self,mode='rb',filetypes=[('All files','*')],title='Choose the  configuration file')
        if file_ != None and file != 'None':
            file_ =abs_path = os.path.abspath(file_.name)
            self.config_file.set(file_)
            self.cfg_str.set(str_default+ '[...%s]'%self.config_file.get()[-15:])
            return None
        else:
            pass

        #return abs_path

    def family_file_search(self):
        str_default = 'Family configuration file'
        file_ = tkFileDialog.askopenfile(parent=self,mode='rb',filetypes=[('All files','*')],title='Choose the  configuration file')
        if file_ != None and file != 'None':
            file_ =abs_path = os.path.abspath(file_.name)
            self.family_file.set(file_)
            self.family_str.set(str_default + '[...%s]'%self.family_str.get()[-15:])
            return None
        else:
            pass

        #return abs_path


    def outfolder_search(self):
        str_default='Outfolder directory'
        try:
            self.outfolder = tkFileDialog.askdirectory(parent=self)

            if self.outfolder != None:
                self.outfolder=os.path.abspath(self.outfolder)
                self.outfolder_string.set(str_default + ' [...%s]'%self.outfolder[-15:])
        except:
            pass

    def gex_file_search(self):
        str_default='Gene exclusion list file'
        file_ = tkFileDialog.askopenfile(parent=self,mode='rb',filetypes=[('txt files','*.txt')],
                title='Choose the  gene exclusion list file',
                initialdir='/users/GD/tools/ediva/Resource/')
        if file_ != None and file != 'None':
            file_ =abs_path = os.path.abspath(file_.name)
            self.gex_file.set(file_)
            self.gex_str.set(str_default+ '[...%s]'%self.gex_file.get()[-15:])
            return None
        else:
            pass

    def white_list_search(self):
        str_default='White list file'
        file_ = tkFileDialog.askopenfile(parent=self,mode='rb',filetypes=[('txt files','*.txt')],
                title='Choose the  white list file')#,
                #initialdir='/users/GD/tools/ediva/Resource/')
        if file_ != None and file != 'None':
            file_ =abs_path = os.path.abspath(file_.name)
            self.white_file.set(file_)
            self.white_str.set(str_default+ '[...%s]'%self.white_file.get()[-15:])
            return None
        else:
            pass


    def ok(self, event=None):
        if not self.validate():
            self.initial_focus.focus_set() # put focus back
            return'whit'
        #self.withdraw()
        self.update_idletasks()
        self.apply()
        if self.result != None:
            self.withdraw()
            self.cancel()
            return self.result
        else:
            pass

    def body(self, master):
        self.configure(background =steel)
        self.resizable(0,0)
        cur_row =0
        #Title
        Label(master, text="eDiVa variant prioritization",background=cobalt,fg=text_color,font=(font_type, txt_dim)).grid(row=cur_row,column=0,columnspan=3,sticky='W'+'E')
        cur_row+=1
        Label(master, text=prioritize_header,background=cobalt,fg=text_color,font=(font_type, txt_dim-6)).grid(row=cur_row,column=0,columnspan=3,sticky='W'+'E')
        cur_row+=1
        Label(master, text="\n",background=steel,fg=text_color,font=(font_type, txt_dim-6)).grid(row=cur_row,column=0,columnspan=3,sticky='W'+'E')
        cur_row+=1

        #Config file selection button
        cur_row+=1
        self.config_file = StringVar()
        self.cfg_str = StringVar()
        self.cfg_str.set('Global configuration file')
        self.cfg_button = Button(master,textvariable=self.cfg_str,command=self.config_file_search,width=50,
                                 background=steel,fg=text_color,font=(font_type, txt_dim),activebackground=cobalt,activeforeground=text_color  )
        self.cfg_button.grid(row=cur_row,column=0,columnspan=3)
        ToolTip.ToolTip(self.cfg_button, follow_mouse=1, text="Select the global configuration file is stored. It contains the paths to all tools used by eDiVa.")

        #Config file selection button
        cur_row+=1
        self.family_file = StringVar()
        self.family_str = StringVar()
        self.family_str.set('Family configuration file')
        fam_button = Button(master,textvariable = self.family_str,command= self.family_file_search ,width=50,
                                 background=steel,fg=text_color,font=(font_type, txt_dim),activebackground= cobalt,activeforeground=text_color  )
        fam_button.grid(row=cur_row,column=0,columnspan=3)
        ToolTip.ToolTip(fam_button, follow_mouse=1, text="Select the family description file with the vcf location for each sample.")

        #Outfolder selection button
        cur_row+=1
        self.outfolder_string = StringVar()
        self.outfolder_string.set('Outfolder directory')
        self.outfolder_ = Button(master,textvariable=self.outfolder_string,command=self.outfolder_search,width=50,
                                 background=steel,fg=text_color,font=(font_type, txt_dim),activebackground= cobalt,activeforeground=text_color  )
        self.outfolder_.grid(row=cur_row,column=0,columnspan=3)
        ToolTip.ToolTip(self.outfolder_, follow_mouse=1, text="Define here the output folder.")


        #Gene_exclusion_list file
        cur_row+=1
        self.gex_file = StringVar()
        self.gex_str = StringVar()
        self.gex_str.set('Gene exclusion list file')
        self.gex_button = Button(master,textvariable=self.gex_str,command=self.gex_file_search,width=50,
                                 background=steel,fg=text_color,font=(font_type, txt_dim),activebackground= cobalt,activeforeground=text_color  )
        self.gex_button.grid(row=cur_row,column=0,columnspan=3)
        ToolTip.ToolTip(self.gex_button, follow_mouse=1, text="Select a .txt file with a list of genes which are surely not related with the case. a eDiVa default file is provided.")

        #White list  file
        cur_row+=1
        self.white_file = StringVar()
        self.white_str = StringVar()
        self.white_str.set('White  list file')
        self.white_button = Button(master,textvariable=self.white_str,command=self.white_list_search,width=50,
                                 background=steel,fg=text_color,font=(font_type, txt_dim),activebackground= cobalt,activeforeground=text_color  )
        self.white_button.grid(row=cur_row,column=0,columnspan=3)
        ToolTip.ToolTip(self.white_button, follow_mouse=1, text="Select a .txt file with a list of genes which are known to be related with the case. This field can be empty")



        cur_row+=1
        Label(master, text="\n",background=steel,fg=text_color,font=(font_type, txt_dim-8)).grid(row=cur_row,column=0,columnspan=3,sticky='W'+'E')
        cur_row+=1
                #Family Type [trio - family]
        cur_row+=1
        Label(master, text="Select the Family type:",background=steel,fg=text_color,font=(font_type, txt_dim)).grid(row=cur_row,column=0,columnspan=1,sticky=W)
        #Combo box button
        self.var_value = StringVar()
        fam_type =ttk.Combobox(master,state='readonly', textvariable=self.var_value,background=cobalt,font=(font_type, txt_dim-2))
        fam_type['values'] = ('family', 'trio')
        fam_type.current(0)
        fam_type.grid(row=cur_row,column=2,sticky='W')

        cur_row+=1
        Label(master, text="\n",background=steel,fg=text_color,font=(font_type, txt_dim-8)).grid(row=cur_row,column=0,columnspan=3,sticky='W'+'E')
        cur_row+=1
        #Inheritance
        cur_row+=1
        Label(master, text="Select Inheritance mode/s:",background=steel,fg=text_color,font=(font_type, txt_dim)).grid(row=cur_row,column=0,columnspan=3)
        cur_row+=1
        self.var_ih1 = IntVar()
        ih1 = Checkbutton(master, text="dominant_inherited", variable=self.var_ih1,background=steel,fg=text_color,selectcolor=steel
                             ,highlightthickness=0,font=(font_type, txt_dim-2),activebackground= cobalt ,activeforeground=text_color )
        ih1.grid(row=cur_row ,column=0)
        self.var_ih2 = IntVar()
        ih2 = Checkbutton(master, text="recessive", variable=self.var_ih2,background=steel,fg=text_color,selectcolor=steel
                             ,highlightthickness=0,font=(font_type, txt_dim-2),activebackground= cobalt ,activeforeground=text_color )
        ih2.grid(row=cur_row ,column=1)
        self.var_ih3 = IntVar()
        ih3 = Checkbutton(master, text="dominant_denovo", variable=self.var_ih3,background=steel,fg=text_color,selectcolor=steel
                             ,highlightthickness=0,font=(font_type, txt_dim-2),activebackground= cobalt ,activeforeground=text_color )
        ih3.grid(row=cur_row ,column=2)

        cur_row+=1
        self.var_ih4 = IntVar()
        ih4 = Checkbutton(master, text="Xlinked", variable=self.var_ih4,background=steel,fg=text_color,selectcolor=steel
                             ,highlightthickness=0,font=(font_type, txt_dim-2),activebackground= cobalt ,activeforeground=text_color )
        ih4.grid(row=cur_row ,column=0,columnspan=2)
        self.var_ih5 = IntVar()
        ih5 = Checkbutton(master, text="compound", variable=self.var_ih5,background=steel,fg=text_color,selectcolor=steel
                             ,highlightthickness=0,font=(font_type, txt_dim-2),activebackground= cobalt ,activeforeground=text_color )
        ih5.grid(row=cur_row ,column=1,columnspan=2)

        #Force and Multisample tickmarks:
        cur_row+=1
        Label(master, text="---------------------------------------------------",background=steel,fg=text_color).grid(row=cur_row,column=0,columnspan=3)
        self.var_f = IntVar()

        cur_row+=1
        self.f = Checkbutton(master, text="Force writing", variable=self.var_f,background=steel,fg=text_color,selectcolor=steel
                             ,highlightthickness=0,font=(font_type, txt_dim-2),activebackground= cobalt ,activeforeground=text_color )
        self.f.grid(row=cur_row,column=0,columnspan=1,sticky='W')
        ToolTip.ToolTip(self.f, follow_mouse=1, text="Select this to force overwriting in case of some output files are already present.")


        self.var_o = IntVar()
        self.o = Checkbutton(master, text="Multisample Input", variable=self.var_o,background=steel,fg=text_color,selectcolor=steel
                             ,highlightthickness=0,font=(font_type, txt_dim-2),activebackground= cobalt ,activeforeground=text_color )
        self.o.grid(row=cur_row,column=2,columnspan=1,sticky='E')
        ToolTip.ToolTip(self.o , follow_mouse=1, text="Select this if the input vcf files have already been merged in a multisample vcf file.")

        #
        ##qsubname and jobname
        cur_row+=1
        Label(master, text="---------------------------------------------------",background=steel,fg=text_color).grid(row=cur_row,column=0,columnspan=3)
        self.qsub = Entry(master, background=cobalt,fg=text_color,font=(font_type, txt_dim-2))
        self.qsub.insert(0,'prioritize.sh')
        self.jname= Entry(master, background=cobalt,fg=text_color,font=(font_type, txt_dim-2))
        self.jname.insert(0,'prioritize_job')

        cur_row+=1
        Label(master, text="Qsub name:", background=steel,fg=text_color,font=(font_type, txt_dim)).grid(row=cur_row,column=0,columnspan=1,sticky='W')
        self.qsub.grid(row=cur_row, column=2,columnspan=2,sticky='E')
        cur_row+=1
        Label(master, text="Job name:", background=steel,fg=text_color,font=(font_type, txt_dim)).grid(row=cur_row,column=0,columnspan=1,sticky='W')
        self.jname.grid(row=cur_row, column=2,columnspan=2,sticky='E')

        #CPU label and buttons
        cur_row+=1
        self.cpu = IntVar()
        self.cpu.set(3)
        self.lcpu = Label(master, textvariable = (self.cpu),width=20, background=steel,fg=text_color,font=(font_type, txt_dim))
        self.nlcpu = Label(master, text = 'CPU',width=20, background=steel,fg=text_color,font=(font_type, txt_dim-4))
        self.lcpu.grid(row=cur_row+1,column=0,columnspan=1)#,sticky='E')
        self.nlcpu.grid(row=cur_row,column=0,columnspan=1)#,sticky='W')
        self.plus_cpu = Button(master,text = '+',command=lambda : self.cpu.set(self.plus1(self.cpu,self.lcpu) ) ,
                                 background=steel,fg=text_color,font=(font_type, txt_dim),activebackground= cobalt,activeforeground=text_color  )
        self.minus_cpu = Button(master,text= '-',command=lambda : self.cpu.set(self.minus1(self.cpu,self.lcpu) ),
                                 background=steel,fg=text_color,font=(font_type, txt_dim),activebackground= cobalt,activeforeground=text_color  )
        self.plus_cpu.grid(row=cur_row+1,column=0,sticky='E')
        self.minus_cpu.grid(row=cur_row+1,column=0,sticky='W')

        #MEM GB label and buttons
        self.mem = IntVar()
        self.mem.set(20)
        self.l_mem = Label(master, textvariable = (self.mem),width=20, background=steel,fg=text_color,font=(font_type, txt_dim))
        self.n_mem = Label(master, text = 'GB Memory',width=20, background=steel,fg=text_color,font=(font_type, txt_dim-4))
        self.l_mem.grid(row=cur_row+1,column=2,columnspan=1)
        self.n_mem.grid(row=cur_row,column=2,columnspan=1)
        self.plus_mem = Button(master,text ='+',command=lambda : self.mem.set(self.plus1(self.mem,self.l_mem) ),
                                 background=steel,fg=text_color,font=(font_type, txt_dim),activebackground= cobalt,activeforeground=text_color  )
        self.minus_mem = Button(master,text='-',command=lambda : self.mem.set(self.minus1(self.mem,self.l_mem) ),
                                 background=steel,fg=text_color,font=(font_type, txt_dim),activebackground= cobalt,activeforeground=text_color  )
        self.plus_mem.grid(row=cur_row+1,column=2,sticky='E')
        self.minus_mem.grid(row=cur_row+1,column=2,sticky='W')


        return self.cfg_button
    def apply(self):
        NewWin = qSubDetails(self.master)
        cpu_mem_params = ',-pe smp %d,-l virtual_free=%dG,-N %s, -e %s/, -o %s/'%(self.cpu.get(),self.mem.get(),self.jname.get(),self.outfolder,self.outfolder)
        if NewWin.result != None:
            self.result = {
                "force"     :self.var_f.get(),
                "multisample":self.var_o.get(),
                "outfolder" :self.outfolder,
                "familytype":self.var_value.get(),
                "qsubname"  :self.qsub.get(),
                "jobname"   :self.jname.get(),
                "inheritance":self.inheritance,
                "qopts"     :NewWin.result+cpu_mem_params,
                "gex_file"  :self.gex_file.get(),
                "white_list":self.white_file.get()
                }

            if len(self.config_file.get())>0and self.config_file.get()!='None' :
                self.result["config"] = self.config_file.get()
            if len(self.family_file.get())>0 and self.family_file.get()!='None' :
                self.result["family"]=self.family_file.get()
        else:
            self.result=None


    def validate(self):
        #validate outfolder
        try:
            dummy = self.outfolder
            if os.path.isdir(dummy):
                pass
            else:
                raise
        except:
            tkMessageBox.showwarning("Prioritize","Please Select one Output Folder")
            return 0
        #valiate cfg_file
        try:
            infile = self.gex_file.get()
            if not(os.path.isfile(infile)):
                raise
        except:
            tkMessageBox.showwarning("Prioritize","Please Select one gene exclusion list file")
            return 0
        #check qsubname
        if len(self.qsub.get())<1 :
            tkMessageBox.showwarning("Prioritize","Qsub must contain text")
            return 0
        #check job name exists
        if len(self.jname.get())<1 :
            tkMessageBox.showwarning("Prioritize","Job name must contain text")
            return 0
        elif self.jname.get()[0].isdigit():
            tkMessageBox.showwarning("Prioritize","Job name must not start with a digit")
            return 0
        #Inheritance:
        or_val = self.var_ih5.get() + self.var_ih4.get() +self.var_ih3.get() +self.var_ih2.get() +self.var_ih1.get()
        if or_val == 0:
            tkMessageBox.showwarning("Prioritize","Choose at least one inheritance mode")
            return 0
        else:
            self.inheritance = list()
            if self.var_ih1.get() >0:
                self.inheritance.append("dominant_inherited")
            if self.var_ih2.get() >0:
                self.inheritance.append("recessive")
            if self.var_ih3.get() >0:
                self.inheritance.append("dominant_denovo")
            if self.var_ih4.get() >0:
                self.inheritance.append("Xlinked")
            if self.var_ih5.get() >0:
                self.inheritance.append("compound")
        return 1
    def buttonbox(self):
        # add standard button box. override if you don't want the
        # standard buttons
        box = Frame(self,background = steel)
        w = Button(box, text="Run", width=10, command=self.ok,fg=text_color,font=(font_type, txt_dim)
                ,highlightbackground= steel,background=  steel,activebackground= cobalt ,activeforeground=text_color        )
        w.pack(side=LEFT, padx=5, pady=5)
        w = Button(box, text="Back", width=10, command=self.cancel,fg=text_color,font=(font_type, txt_dim)
                ,highlightbackground= steel,background=  steel,activebackground= cobalt,activeforeground=text_color       )
        w.pack(side=LEFT, padx=5, pady=5)
        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)
        box.pack()

class qSubDetails(Dialog):
    '''
    Window to gather the needed qsub parameters: the output is a string to be added to the pipeline call
    '''
    def plus1(self,var,l):
        var.set(int(var.get())+1)
        l.textvariable = var
        return var.get()

    def minus1(self,var,l):
        var.set(max(0,int(var.get())-1))
        l.textvariable = var
        return var.get()


    def body(self, master):
        self.configure(background=steel)
        self.resizable(0,0)
        #Title
        Label(master, text="\nQsub Options:\n",background=UI_blue,fg=text_color,font=(font_type, txt_dim)).grid(row=0,column=0,columnspan=3,sticky='W'+'E')

        #Select Queue
        Label(master, text="Select the qsub queue:",background=steel,fg=text_color,font=(font_type, txt_dim)).grid(row=2,column=0,columnspan=1,sticky=W+E)
        #Combo box button
        self.var_value = StringVar()
        self.queue =ttk.Combobox(master,state='readonly', textvariable=self.var_value)
        self.queue.configure(background=steel)
        self.queue['values'] = ('short-sl65', 'long-sl65', 'mem_256','mem_512','so-el6','xe-el6','cn-el6','fk-el6','pr-el6','rg-el6','tg-el6')
        self.queue.current(1)
        self.queue.grid(row=2,column=1,columnspan=2,sticky=W+E)


        #Select the Time
        lbl_width = 20
        self.hh = IntVar()
        self.hh.set(6)
        self.l_hh = Label(master, textvariable = (self.hh),width=lbl_width, background=steel,fg=text_color,font=(font_type, txt_dim))
        n_hh = Label(master, text = 'HH',width=lbl_width, background=steel,fg=text_color,font=(font_type, txt_dim-4))
        self.l_hh.grid(row=7,column=0,columnspan=1)
        n_hh.grid(row=6,column=0,columnspan=1)
        self.plus_hh = Button(master,text ='+',command=lambda : self.hh.set(self.plus1(self.hh,self.l_hh) ),
                                 background=steel,fg=text_color,font=(font_type, txt_dim),activebackground= UI_blue,activeforeground=text_color  )
        self.minus_hh = Button(master,text='-',command=lambda : self.hh.set(self.minus1(self.hh,self.l_hh) ),
                                 background=steel,fg=text_color,font=(font_type, txt_dim),activebackground=UI_blue,activeforeground=text_color  )
        self.plus_hh.grid(row=7,column=0,sticky='E')
        self.minus_hh.grid(row=7,column=0,sticky='W')

        self.mm = IntVar()
        self.mm.set(0)
        self.l_mm = Label(master, textvariable = (self.mm),width=lbl_width, background=steel,fg=text_color,font=(font_type, txt_dim))
        n_mm = Label(master, text = 'MM',width=lbl_width, background=steel,fg=text_color,font=(font_type, txt_dim-4))
        self.l_mm.grid(row=7,column=1,columnspan=1)
        n_mm.grid(row=6,column=1,columnspan=1)
        self.plus_mm = Button(master,text ='+',command=lambda : self.mm.set(self.plus1(self.mm,self.l_mm) ),
                                 background=steel,fg=text_color,font=(font_type, txt_dim),activebackground= UI_blue,activeforeground=text_color  )
        self.minus_mm = Button(master,text='-',command=lambda : self.mm.set(self.minus1(self.mm,self.l_mm) ),
                                 background=steel,fg=text_color,font=(font_type, txt_dim),activebackground=UI_blue,activeforeground=text_color  )
        self.plus_mm.grid(row=7,column=1,sticky='E')
        self.minus_mm.grid(row=7,column=1,sticky='W')

        self.ss = IntVar()
        self.ss.set(0)
        self.l_ss = Label(master, textvariable = (self.ss),width=lbl_width, background=steel,fg=text_color,font=(font_type, txt_dim))
        n_ss = Label(master, text = 'SS',width=lbl_width, background=steel,fg=text_color,font=(font_type, txt_dim-4))
        self.l_ss.grid(row=7,column=2,columnspan=1)
        n_ss.grid(row=6,column=2,columnspan=1)
        self.plus_ss = Button(master,text ='+',command=lambda : self.ss.set(self.plus1(self.ss,self.l_ss) ),
                                 background=steel,fg=text_color,font=(font_type, txt_dim),activebackground= UI_blue,activeforeground=text_color  )
        self.minus_ss = Button(master,text='-',command=lambda : self.ss.set(self.minus1(self.ss,self.l_ss) ),
                                 background=steel,fg=text_color,font=(font_type, txt_dim),activebackground=UI_blue,activeforeground=text_color  )
        self.plus_ss.grid(row=7,column=2,sticky='E')
        self.minus_ss.grid(row=7,column=2,sticky='W')

        #Mail Options
        self.var_a= IntVar()
        self.var_a.set(1)
        self.a = Checkbutton(master, text="Abort", variable=self.var_a,background=steel,fg=text_color,selectcolor=steel
                             ,highlightthickness=0,font=(font_type, txt_dim),activebackground= UI_blue ,activeforeground=text_color    )
        self.a.grid(row=5,column=0,columnspan=1,sticky=W)

        self.var_b = IntVar()
        self.var_b.set(1)
        self.b = Checkbutton(master, text="Begin", variable=self.var_b,background=steel,fg=text_color,selectcolor=steel,
                             highlightthickness=0,font=(font_type, txt_dim),activebackground= UI_blue,activeforeground=text_color     )
        self.b.grid(row=5,column=1,columnspan=1)

        self.var_e = IntVar()
        self.var_e.set(1)
        self.e = Checkbutton(master, text="End", variable=self.var_e,background=steel,fg=text_color,selectcolor=steel,
                             highlightthickness=0,font=(font_type, txt_dim),activebackground= UI_blue,activeforeground=text_color     )
        self.e.grid(row=5,column=2,columnspan=1,sticky=E)

        self.var_custom = IntVar()
        self.var_custom.set(1)
        self.custom = Checkbutton(master, text="Custom Mail annotation", variable=self.var_custom ,background=steel,fg=text_color,selectcolor=steel,
                             highlightthickness=0,font=(font_type, txt_dim),activebackground= UI_blue,activeforeground=text_color     )
        self.custom .grid(row=4,column=1,columnspan=1,sticky=E)

        self.email = Entry(master, background=UI_blue,fg=text_color,font=(font_type, txt_dim-2))
        self.email.insert(0,'xxx@crg.es')
        Label(master, text="Mail Address:", background=steel,fg=text_color,font=(font_type, txt_dim)).grid(row=8,column=0,columnspan=1,sticky=W)
        self.email.grid(row =8, column=1)

        return self.queue # initial focus
    def apply(self):
        qsub        = self.var_value.get()
        a           = self.var_a.get()
        b           = self.var_b.get()
        e           = self.var_e.get()
        custom      = self.var_custom.get()
        hh          = self.hh.get()
        mm          = self.mm.get()
        ss          = self.ss.get()
        time        = 60*60*hh + 60*mm + ss
        email       = self.email.get()

        abe = ""
        if a==1:
            abe+='a'
        if b==1:
            abe+='b'
        if e==1:
            abe+='e'

        self.result = ("-hard,-q %s,-M %s,-l h_rt=%d"%(qsub,email,time))
        if custom ==1:
            self.result+=',-v CUSTOM_EMAIL=yes'
        if len(abe)>0:
            self.result+=',-m %s'%abe
        self.result+=",-b y, -shell y,-cwd"
        #print self.result


    def validate(self):
        return 1

    def buttonbox(self):
        # add standard button box. override if you don't want the
        # standard buttons
        box = Frame(self,background = steel)
        w = Button(box, text="Run", width=10, command=self.ok,fg=text_color,font=(font_type, txt_dim)
                ,highlightbackground= steel,background=  steel,activebackground= UI_blue,activeforeground=text_color)
        w.pack(side=LEFT, padx=5, pady=5)
        w = Button(box, text="Back", width=10, command=self.cancel,fg=text_color,font=(font_type, txt_dim)
                ,highlightbackground= steel,background=  steel,activebackground= UI_blue,activeforeground=text_color)
        w.pack(side=LEFT, padx=5, pady=5)
        self.bind("<Return>", self.ok)
        #self.bind("<Escape>", self.cancel)
        box.pack()


    def ok(self, event=None):
        if not self.validate():
            self.initial_focus.focus_set() # put focus back
            return'whit'
        self.withdraw()
        self.update_idletasks()
        self.apply()
        self.cancel()
        return self.result



class Sample_list_window(Dialog):
    '''
    Window to gather the sample_list to send to prediction

    '''


    def body(self, master):
        self.configure(background=steel)
        self.resizable(0,0)
        #Title


        samples = self.master.sample_list

        self.queue =Label(master, text="\nSamples to process:\n",background=UI_blue,fg=text_color,font=(font_type, txt_dim)).grid(row=0,column=1,columnspan=3,sticky='W'+'E')

        #From list to dictionary:
        self.sample_list = dict()
        for i in samples:
            self.sample_list[i] = 1
        self.var_value = StringVar()
        i = 1
        #Populate now with the samples checkboxes
        for sample in self.sample_list:
            i+=1
            self.sample_list[sample] = IntVar()
            self.sample_list[sample].set(1)
            l = Checkbutton(master, text=sample, variable=self.sample_list[sample] ,background=steel,fg=text_color,selectcolor=steel,
                 highlightthickness=0,font=(font_type, txt_dim),activebackground= UI_blue,activeforeground=text_color     )
            l.grid(row=i,columnspan=3,sticky='W')
        return self.queue# initial focus

    def apply(self):
        #scan the checkboxes and return the string
        return_string = list()
        keys = self.sample_list.keys()
        for k in keys:
            value = self.sample_list.get(k)
            if value.get() >0:
                return_string.append(k)
        return_string = ','.join(return_string)
        self.result = return_string
        #print self.result
        if len(return_string)==0:
            self.result=None
            return None

    def validate(self):
        return 1

    def buttonbox(self):
        # add standard button box. override if you don't want the
        # standard buttons
        box = Frame(self,background = steel)
        w = Button(box, text="Ok", width=10, command=self.ok,fg=text_color,font=(font_type, txt_dim)
                ,highlightbackground= steel,background=  steel,activebackground= UI_blue,activeforeground=text_color)
        w.pack(side=LEFT, padx=5, pady=5)
        w = Button(box, text="Back", width=10, command=self.cancel,fg=text_color,font=(font_type, txt_dim)
                ,highlightbackground= steel,background=  steel,activebackground= UI_blue,activeforeground=text_color)
        w.pack(side=LEFT, padx=5, pady=5)
        self.bind("<Return>", self.ok)
        #self.bind("<Escape>", self.cancel)
        box.pack()


    def ok(self, event=None):
        if not self.validate():
            self.initial_focus.focus_set() # put focus back
            return'whit'
        self.withdraw()
        self.update_idletasks()
        self.apply()
        self.cancel()
        return self.result






class Resume_window(Dialog):
    '''
    Window to gather the needed Annotate.py values : the output is a dictionary with variables to launch the command.
    '''
    def file_search(self):
        str_default='Config. File'
        self.file = tkFileDialog.askopenfile(parent=self,mode='rb',filetypes=[('Job_queue files','*.qlist'), ('Single pipeline files','*.pipe')],title='Choose the file to run')
        if self.file != None:
            self.file =abs_path = os.path.abspath(self.file.name)
            self.filename_string.set(str_default + ' [...%s]'%self.file[-15:])
        else:
            self.file = None

        #return abs_path

    def plus1(self,var,l):
        var.set(int(var.get())+1)
        l.textvariable = var
        return var.get()

    def minus1(self,var,l):
        var.set(max(1,int(var.get())-1))
        l.textvariable = var
        return var.get()

    def body(self, master):
        #Title
        cur_row=0
        self.configure(background =steel)
        self.resizable(0,0)
        Label(master, text="\nRelaunch job:\n",background=UI_blue,fg=text_color,font=(font_type, txt_dim)).grid(row=cur_row,column=0,columnspan=3,sticky='W'+'E')


        #Config file selection button
        cur_row+=1
        self.filename_string = StringVar()
        self.filename_string.set('Select file to run')
        self.filename = Button(master,textvariable=self.filename_string,command=self.file_search,width=50,background=steel,fg=text_color,
                               font=(font_type, txt_dim),activebackground= UI_blue,activeforeground=text_color  )
        self.filename.grid(row=cur_row,column=0,columnspan=3)


        #CPU label and buttons
        cur_row+=1
        self.cpu = IntVar()
        self.cpu.set(4)
        self.lcpu = Label(master, textvariable = (self.cpu),width=20, background=steel,fg=text_color,font=(font_type, txt_dim))
        self.nlcpu = Label(master, text = 'CPU',width=20, background=steel,fg=text_color,font=(font_type, txt_dim-4))
        self.lcpu.grid(row=cur_row+1,column=0,columnspan=1)#,sticky='E')
        self.nlcpu.grid(row=cur_row,column=0,columnspan=1)#,sticky='W')
        self.plus_cpu = Button(master,text = '+',command=lambda : self.cpu.set(self.plus1(self.cpu,self.lcpu) ) ,
                                 background=steel,fg=text_color,font=(font_type, txt_dim),activebackground= UI_blue,activeforeground=text_color  )
        self.minus_cpu = Button(master,text= '-',command=lambda : self.cpu.set(self.minus1(self.cpu,self.lcpu) ),
                                 background=steel,fg=text_color,font=(font_type, txt_dim),activebackground= UI_blue,activeforeground=text_color  )
        self.plus_cpu.grid(row=cur_row+1,column=0,sticky='E')
        self.minus_cpu.grid(row=cur_row+1,column=0,sticky='W')

        #MEM GB label and buttons
        self.mem = IntVar()
        self.mem.set(12)
        self.l_mem = Label(master, textvariable = (self.mem),width=20, background=steel,fg=text_color,font=(font_type, txt_dim))
        self.n_mem = Label(master, text = 'GB Memory',width=20, background=steel,fg=text_color,font=(font_type, txt_dim-4))
        self.l_mem.grid(row=cur_row+1,column=2,columnspan=1)
        self.n_mem.grid(row=cur_row,column=2,columnspan=1)
        self.plus_mem = Button(master,text ='+',command=lambda : self.mem.set(self.plus1(self.mem,self.l_mem) ),
                                 background=steel,fg=text_color,font=(font_type, txt_dim),activebackground= UI_blue,activeforeground=text_color  )
        self.minus_mem = Button(master,text='-',command=lambda : self.mem.set(self.minus1(self.mem,self.l_mem) ),
                                 background=steel,fg=text_color,font=(font_type, txt_dim),activebackground= UI_blue,activeforeground=text_color  )
        self.plus_mem.grid(row=cur_row+1,column=2,sticky='E')
        self.minus_mem.grid(row=cur_row+1,column=2,sticky='W')

        return self.filename # initial focus


    def buttonbox(self):
        # add standard button box. override if you don't want the
        # standard buttons
        box = Frame(self,background = steel)
        w = Button(box, text="Run", width=10, command=self.ok,fg=text_color,font=(font_type, txt_dim)
                ,highlightbackground= steel,background=  steel,activebackground= UI_blue ,activeforeground=text_color        )
        w.pack(side=LEFT, padx=5, pady=5)
        w = Button(box, text="Back", width=10, command=self.cancel,fg=text_color,font=(font_type, txt_dim)
                ,highlightbackground= steel,background=  steel,activebackground= UI_blue,activeforeground=text_color       )
        w.pack(side=LEFT, padx=5, pady=5)
        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)
        box.pack()

    def ok(self, event=None):
        if not self.validate():
            self.initial_focus.focus_set() # put focus back
            return'whit'
        #
        self.update_idletasks()
        self.apply()
        if self.result != None:
            self.cancel()
            #self.withdraw()
            return self.result
        else:
            return None

    def apply(self):
        NewWin = qSubDetails(self.master)
        cpu_mem_params = ',-pe smp %d,-l virtual_free=%dG '%(self.cpu.get(),self.mem.get())
        if NewWin.result != None:
            self.result = {
                "mem"       :self.mem.get(),
                "cpu"       :self.cpu.get(),
                "filename"  :self.file,
                "qopts"     :NewWin.result+cpu_mem_params

                }
            #print self.result
        else :
            self.result = None


    def validate(self):
        #valiate file
        try:
            infile = self.file
            if not(os.path.isfile(infile)):
                raise
        except:
            tkMessageBox.showwarning("Predict","Please Select one Input File")
            return 0
        return 1



class Kill_window(Dialog):
    '''

    '''
    def file_search(self):
        str_default='Job to be killed '
        self.file = tkFileDialog.askopenfile(parent=self,mode='rb',filetypes=[('Job_queue files','*.qlist')],title='Choose the jobs to stop')
        if self.file != None:
            self.file =abs_path = os.path.abspath(self.file.name)
            self.filename_string.set(str_default + ' [...%s]'%self.file[-15:])
        else:
            self.file = None

        #return abs_path

    def body(self, master):
        #Title
        cur_row=0
        self.configure(background =steel)
        self.resizable(0,0)
        Label(master, text="\nRelaunch job:\n",background=UI_blue,fg=text_color,font=(font_type, txt_dim)).grid(row=cur_row,column=0,columnspan=3,sticky='W'+'E')

        #Config file selection button
        cur_row+=1
        self.filename_string = StringVar()
        self.filename_string.set('Select job-queue to stop')
        self.filename = Button(master,textvariable=self.filename_string,command=self.file_search,width=50,background=steel,fg=text_color,
                               font=(font_type, txt_dim),activebackground= UI_blue,activeforeground=text_color  )
        self.filename.grid(row=cur_row,column=0,columnspan=3)
        return self.filename # initial focus


    def buttonbox(self):
        # add standard button box. override if you don't want the
        # standard buttons
        box = Frame(self,background = steel)
        w = Button(box, text="Stop", width=10, command=self.ok,fg=text_color,font=(font_type, txt_dim)
                ,highlightbackground= steel,background=  steel,activebackground= UI_blue ,activeforeground=text_color        )
        w.pack(side=LEFT, padx=5, pady=5)
        w = Button(box, text="Back", width=10, command=self.cancel,fg=text_color,font=(font_type, txt_dim)
                ,highlightbackground= steel,background=  steel,activebackground= UI_blue,activeforeground=text_color       )
        w.pack(side=LEFT, padx=5, pady=5)
        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)
        box.pack()

    def ok(self, event=None):
        if not self.validate():
            self.initial_focus.focus_set() # put focus back
            return'whit'
        #
        self.update_idletasks()
        self.apply()
        if self.result != None:
            self.cancel()
            #self.withdraw()
            return self.result
        else:
            return None

    def apply(self):
        self.result = {
            "filename" :self.file
            }
        return None

    def validate(self):
        #valiate file
        try:
            infile = self.file
            if not(os.path.isfile(infile)):
                raise
        except:
            tkMessageBox.showwarning("Predict","Please Select one Input File")
            return 0
        return 1
