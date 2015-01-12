#!/usr/bin/env python

from Tkinter import *
import ttk
import os
import tkFileDialog
import tkMessageBox


#colors
back_color = "#aa00ff"
emerald ='#008a00'
#orange ='#fa6800'
orange ='#da532c'#e8980e'#'#cd3700'
cobalt ='#0050ef'
steel ="#1d1d1d"#647687'
green ='#60a917'
mango = '#f09609'
azzurro = '#1ba1f2'
txt_dim =15

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

        self.grab_set()

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
        Label(master, text="\nAnnotate Options:\n",background=orange,fg='white',font=(font_type, txt_dim)).grid(row=0,column=0,columnspan=3,sticky='W'+'E')
        
        
        #File selection button
        self.filename = Button(master,text='Input File',command=self.file_search,width=50,fg='white',font=(font_type, txt_dim)
                ,highlightbackground= steel,background=  steel,activebackground= orange,activeforeground='white'         )        
        self.filename.grid(row=1,column=0,columnspan=3)
        
        #Select Variant Type
        Label(master, text="Select the variant type:",background=steel,fg='white',font=(font_type, txt_dim)).grid(row=2,column=0,columnspan=1,sticky=W+E)
        #Combo box button
        self.var_value = StringVar()
        self.var =ttk.Combobox(master,state='readonly', textvariable=self.var_value)
        self.var.configure(background=steel)
        self.var['values'] = ('snp', 'indel', 'all')
        self.var.current(2)
        self.var.grid(row=2,column=1,columnspan=2,sticky=W+E)
        
        #SelectGeneMode ensGene,refGene,knownGene,all
        Label(master, text="Select the Gene Definition type:",background=steel,fg='white',font=(font_type, txt_dim)).grid(row=3,column=0,columnspan=1,sticky=W)
        #Combo box button
        self.gdef_value = StringVar()
        self.gdef =ttk.Combobox(master,state='readonly', textvariable=self.gdef_value,background=steel)
        self.gdef['values'] = ('ensGene','refGene','knownGene','all')
        self.gdef.current(2)
        self.gdef.grid(row=3,column=2,columnspan=2)
        
        #SelectGeneDefinition
        Label(master, text="Select the Genotype Mode:",background=steel,fg='white',font=(font_type, txt_dim)).grid(row=4,column=0,columnspan=1,sticky=W)
        #Combo box button
        self.gtype_value = StringVar()
        self.gtype =ttk.Combobox(master,state='readonly', textvariable=self.gtype_value,background=steel)
        self.gtype['values'] = ('compact','complete')
        self.gtype.current(0)
        self.gtype.grid(row=4,column=2,columnspan=2)
        
        #Force and Multisample tickmarks:
        self.var_f = IntVar()
        self.f = Checkbutton(master, text="Force writing?", variable=self.var_f,background=steel,fg='white',selectcolor=steel
                             ,highlightthickness=0,font=(font_type, txt_dim),activebackground= orange ,activeforeground='white'    )
        self.f.grid(row=5,column=0,columnspan=1,sticky=W)
       
        self.var_o = IntVar()
        self.o = Checkbutton(master, text="onlyGenicAnnotation?", variable=self.var_o,background=steel,fg='white',selectcolor=steel,
                             highlightthickness=0,font=(font_type, txt_dim),activebackground= orange,activeforeground='white'     )
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
        w = Button(box, text="Run", width=10, command=self.ok,fg='white',font=(font_type, txt_dim)
                ,highlightbackground= steel,background=  steel,activebackground= orange,activeforeground='white'         )
        w.pack(side=LEFT, padx=5, pady=5)
        w = Button(box, text="Back", width=10, command=self.cancel,fg='white',font=(font_type, txt_dim)
                ,highlightbackground= steel,background=  steel,activebackground= orange,activeforeground='white'         )
        w.pack(side=LEFT, padx=5, pady=5)
        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)
        box.pack()    
class Predict_window(Dialog):
    '''
    Window to gather the needed Annotate.py values : the output is a dictionary with variables to launch the command.
    '''
    def file_search(self):
        self.file = tkFileDialog.askopenfile(parent=self,mode='rb',filetypes=[('All files','*')],title='Choose the  configuration file')
        if self.file != None:
            self.file =abs_path = os.path.abspath(self.file.name)
        else:
            self.file = None
        
        #return abs_path
        
    def infolder_search(self):
        try:
            self.infolder = tkFileDialog.askdirectory(parent=self)
            if self.infolder != None and os.path.isdir(self.infolder):
                self.infolder=os.path.abspath(self.infolder)
                self.infolder_string.set('In_dir : %s'%self.infolder)
        except:
            self.infolder=None
        
    
    def outfolder_search(self):
        try:
            self.outfolder = tkFileDialog.askdirectory(parent=self)
            if self.outfolder != None and os.path.isdir(self.outfolder):
                self.outfolder=os.path.abspath(self.outfolder)
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
        self.configure(background =steel)
        self.resizable(0,0)
        Label(master, text="\nPredict Options:\n",background=emerald,fg='white',font=(font_type, txt_dim)).grid(row=0,column=0,columnspan=3,sticky='W'+'E')
        
        #Infolder selection button
        self.infolder_string = StringVar()
        self.infolder_string.set('Input directory')
        self.infolder_ = Button(master,textvariable=self.infolder_string,command=self.infolder_search,width=50,background=steel,fg='white',
                                font=(font_type, txt_dim),activebackground= emerald,activeforeground='white'  )
        #self.infolder_ = Button(master,text='Input directory',command=self.infolder_search,width=50,background=steel,fg='white',
        #                font=(font_type, txt_dim),activebackground= emerald,activeforeground='white'  )  
        self.infolder_.grid(row=1,column=0,columnspan=3)
        
        #Outfolder selection button
        self.outfolder_ = Button(master,text='Outfolder directory',command=self.outfolder_search,width=50,
                                 background=steel,fg='white',font=(font_type, txt_dim),activebackground= emerald,activeforeground='white'  )        
        self.outfolder_.grid(row=2,column=0,columnspan=3)
        
        #Config file selection button
        self.filename = Button(master,text='Config. File',command=self.file_search,width=50,background=steel,fg='white',
                               font=(font_type, txt_dim),activebackground= emerald,activeforeground='white'  )        
        self.filename.grid(row=3,column=0,columnspan=3)
        
        
        #Namestart label and buttons
        self.var = IntVar()
        self.var.set(1)
        self.l = Label(master, textvariable = (self.var),width=20, background=steel,fg='white',font=(font_type, txt_dim))
        self.nl = Label(master, text = 'Namestart',width=20, background=steel,fg='white',font=(font_type, txt_dim-4))
        self.l.grid(row=5,column=0,columnspan=1)#,sticky='E')
        self.nl.grid(row=4,column=0,columnspan=1)#,sticky='W')
        self.plus_n = Button(master,text = '+',command=lambda : self.var.set(self.plus1(self.var,self.l) ) ,
                                 background=steel,fg='white',font=(font_type, txt_dim),activebackground= emerald,activeforeground='white'  )
        self.minus_n = Button(master,text= '-',command=lambda : self.var.set(self.minus1(self.var,self.l) ),
                                 background=steel,fg='white',font=(font_type, txt_dim),activebackground= emerald,activeforeground='white'  )
        self.plus_n.grid(row=5,column=0,sticky='E')
        self.minus_n.grid(row=5,column=0,sticky='W')
        Label(master, text="",width=10, background=steel,fg='white',font=(font_type, txt_dim-4)).grid(row=5,column=1)
        #Namelength label and buttons
        self.var_length = IntVar()
        self.var_length.set(1)
        self.length = Label(master, textvariable = (self.var_length),width=20, background=steel,fg='white',font=(font_type, txt_dim))
        self.nlength = Label(master, text = 'Namelength',width=20, background=steel,fg='white',font=(font_type, txt_dim-4))
        self.length.grid(row=5,column=2,columnspan=1)
        self.nlength.grid(row=4,column=2,columnspan=1)
        self.plus_n_length = Button(master,text ='+',command=lambda : self.var_length.set(self.plus1(self.var_length,self.length) ),
                                 background=steel,fg='white',font=(font_type, txt_dim),activebackground= emerald,activeforeground='white'  )
        self.minus_n_length = Button(master,text='-',command=lambda : self.var_length.set(self.minus1(self.var_length,self.length) ),
                                 background=steel,fg='white',font=(font_type, txt_dim),activebackground= emerald,activeforeground='white'  )
        self.plus_n_length.grid(row=5,column=2,sticky='E')
        self.minus_n_length.grid(row=5,column=2,sticky='W')

        #CPU label and buttons
        self.cpu = IntVar()
        self.cpu.set(4)
        self.lcpu = Label(master, textvariable = (self.cpu),width=20, background=steel,fg='white',font=(font_type, txt_dim))
        self.nlcpu = Label(master, text = 'CPU',width=20, background=steel,fg='white',font=(font_type, txt_dim-4))
        self.lcpu.grid(row=7,column=0,columnspan=1)#,sticky='E')
        self.nlcpu.grid(row=6,column=0,columnspan=1)#,sticky='W')
        self.plus_cpu = Button(master,text = '+',command=lambda : self.cpu.set(self.plus1(self.cpu,self.lcpu) ) ,
                                 background=steel,fg='white',font=(font_type, txt_dim),activebackground= emerald,activeforeground='white'  )
        self.minus_cpu = Button(master,text= '-',command=lambda : self.cpu.set(self.minus1(self.cpu,self.lcpu) ),
                                 background=steel,fg='white',font=(font_type, txt_dim),activebackground= emerald,activeforeground='white'  )
        self.plus_cpu.grid(row=7,column=0,sticky='E')
        self.minus_cpu.grid(row=7,column=0,sticky='W')

        #MEM GB label and buttons
        self.mem = IntVar()
        self.mem.set(12)
        self.l_mem = Label(master, textvariable = (self.mem),width=20, background=steel,fg='white',font=(font_type, txt_dim))
        self.n_mem = Label(master, text = 'GB Memory',width=20, background=steel,fg='white',font=(font_type, txt_dim-4))
        self.l_mem.grid(row=7,column=2,columnspan=1)
        self.n_mem.grid(row=6,column=2,columnspan=1)
        self.plus_mem = Button(master,text ='+',command=lambda : self.mem.set(self.plus1(self.mem,self.l_mem) ),
                                 background=steel,fg='white',font=(font_type, txt_dim),activebackground= emerald,activeforeground='white'  )
        self.minus_mem = Button(master,text='-',command=lambda : self.mem.set(self.minus1(self.mem,self.l_mem) ),
                                 background=steel,fg='white',font=(font_type, txt_dim),activebackground= emerald,activeforeground='white'  )
        self.plus_mem.grid(row=7,column=2,sticky='E')
        self.minus_mem.grid(row=7,column=2,sticky='W')
        
        #Select Indel Caller
        Label(master, text="Select the Indel caller type:",background=steel,fg='white',font=(font_type, txt_dim)).grid(row=8,column=0,columnspan=2,sticky=W)
        #Combo box button
        self.var_value = StringVar()
        self.var_indel =ttk.Combobox(master,state='readonly', textvariable=self.var_value, background=emerald,font=(font_type, txt_dim-4))
        self.var_indel['values'] = ('gatk', 'clindel', 'both')
        self.var_indel.current(0)
        self.var_indel.grid(row=8,column=2,sticky=W+E)

        #Force tick.
        self.var_f = IntVar()
        Label(master, text="Do you want to fuse the outputs:", background=steel,fg='white',font=(font_type, txt_dim)
              ).grid(row=9,column=0,columnspan=2,sticky=W)
        self.f = Checkbutton(master, text="Fusevariants?", variable=self.var_f,background=steel,fg='white',selectcolor=steel
                             ,highlightthickness=0,font=(font_type, txt_dim),activebackground= emerald ,activeforeground='white' )
        self.f.grid(row=9,column=2,columnspan=2)

        #Read extensions
        self.read1 = Entry(master, background=emerald,fg='white',font=(font_type, txt_dim-2))
        self.read2 = Entry(master, background=emerald,fg='white',font=(font_type, txt_dim-2))
        self.read1.insert(0,'read1.fastq.gz')
        self.read2.insert(0,'read2.fastq.gz')
        Label(master, text="1st read extension:", background=steel,fg='white',font=(font_type, txt_dim)).grid(row=10,column=0,columnspan=1,sticky=W)
        Label(master, text="2nd read extension:", background=steel,fg='white',font=(font_type, txt_dim)).grid(row=11,column=0,columnspan=1,sticky=W)
        self.read1.grid(row=10, column=1,columnspan=2,sticky=W+E)
        self.read2.grid(row=11, column=1,columnspan=2,sticky=W+E)
        
        #qsubname and max_coverage
        self.qsub = Entry(master, background=emerald,fg='white',font=(font_type, txt_dim-2))
        self.qsub.insert(0,'predict.sh')
        self.max_cov= Entry(master, background=emerald,fg='white',font=(font_type, txt_dim-2))
        self.max_cov.insert(0,'400')
        Label(master, text="Qsub name:", background=steel,fg='white',font=(font_type, txt_dim)).grid(row=12,column=0,columnspan=1,sticky=W)
        Label(master, text="maximum coverage:", background=steel,fg='white',font=(font_type, txt_dim)).grid(row=13,column=0,columnspan=1,sticky=W)
        self.qsub.grid(row=12, column=1,columnspan=2,sticky=W+E)
        self.max_cov.grid(row=13, column=1,columnspan=2,sticky=W+E)
        
        return self.infolder_ # initial focus
    def buttonbox(self):
        # add standard button box. override if you don't want the
        # standard buttons
        box = Frame(self,background = steel)
        w = Button(box, text="Run", width=10, command=self.ok,fg='white',font=(font_type, txt_dim)
                ,highlightbackground= steel,background=  steel,activebackground= emerald ,activeforeground='white'        )
        w.pack(side=LEFT, padx=5, pady=5)
        w = Button(box, text="Back", width=10, command=self.cancel,fg='white',font=(font_type, txt_dim)
                ,highlightbackground= steel,background=  steel,activebackground= emerald,activeforeground='white'       )
        w.pack(side=LEFT, padx=5, pady=5)
        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)
        box.pack()
        
    def ok(self, event=None):
        if not self.validate():
            self.initial_focus.focus_set() # put focus back
            return'whit'
        #self.withdraw()
        self.update_idletasks()
        self.apply()
        if self.result != None:
            self.cancel()
            return self.result
        else:
            pass
        
    def apply(self):
        NewWin = qSubDetails(self.master)
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
                "qopts"     :NewWin.result
                }
            print self.result
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
    def file_search(self):
        file_ = tkFileDialog.askopenfile(parent=self,mode='rb',filetypes=[('All files','*')],title='Choose the  configuration file')
        if file_ != None and file != 'None':
            file_ =abs_path = os.path.abspath(file_.name)
            return file_
        else:
            pass
        
        #return abs_path
        
        
    
    def outfolder_search(self):
        try:
            self.outfolder = tkFileDialog.askdirectory(parent=self)
                
            if self.outfolder != None:
                self.outfolder=os.path.abspath(self.outfolder)
        except:
            pass

    def ok(self, event=None):
        if not self.validate():
            self.initial_focus.focus_set() # put focus back
            return'whit'
        #self.withdraw()
        self.update_idletasks()
        self.apply()
        if self.result != None:
            self.cancel()
            return self.result
        else:
            pass
        
    def body(self, master):
        self.configure(background =steel)
        self.resizable(0,0)
        #Title
        Label(master, text="\nPrioritize Options:\n",background=cobalt,fg='white',font=(font_type, txt_dim)).grid(row=0,column=0,columnspan=3,sticky='W'+'E')
        
        #Config file selection button
        self.config_file = StringVar()
        self.cfg_button = Button(master,text='Global configuration file',command= lambda : self.config_file.set(self.file_search()) ,width=50,
                                 background=steel,fg='white',font=(font_type, txt_dim),activebackground= cobalt,activeforeground='white'  )       
        self.cfg_button.grid(row=1,column=0,columnspan=3)
        
        #Config file selection button
        self.family_file = StringVar()
        fam_button = Button(master,text='Family configuration file',command= lambda : self.family_file.set(self.file_search()) ,width=50,
                                 background=steel,fg='white',font=(font_type, txt_dim),activebackground= cobalt,activeforeground='white'  )       
        fam_button.grid(row=2,column=0,columnspan=3)
       
        #Outfolder selection button
        self.outfolder_ = Button(master,text='Outfolder directory',command=self.outfolder_search,width=50,
                                 background=steel,fg='white',font=(font_type, txt_dim),activebackground= cobalt,activeforeground='white'  )        
        self.outfolder_.grid(row=3,column=0,columnspan=3)
        
        #Family Type [trio - family]
        Label(master, text="Select the Family type:",background=steel,fg='white',font=(font_type, txt_dim)).grid(row=4,column=0,columnspan=1,sticky=W)
        #Combo box button
        self.var_value = StringVar()
        fam_type =ttk.Combobox(master,state='readonly', textvariable=self.var_value,background=cobalt,font=(font_type, txt_dim-2))
        fam_type['values'] = ('family', 'trio')
        fam_type.current(0)
        fam_type.grid(row=4,column=2,sticky='W')
        
        #Inheritance
        Label(master, text="Select Inheritance mode/s:",background=steel,fg='white',font=(font_type, txt_dim)).grid(row=5,column=0,columnspan=3)
        self.var_ih1 = IntVar()
        ih1 = Checkbutton(master, text="dominant_inherited", variable=self.var_ih1,background=steel,fg='white',selectcolor=steel
                             ,highlightthickness=0,font=(font_type, txt_dim-2),activebackground= cobalt ,activeforeground='white' )
        ih1.grid(row=6 ,column=0)
        self.var_ih2 = IntVar()
        ih2 = Checkbutton(master, text="recessive", variable=self.var_ih2,background=steel,fg='white',selectcolor=steel
                             ,highlightthickness=0,font=(font_type, txt_dim-2),activebackground= cobalt ,activeforeground='white' )
        ih2.grid(row=6 ,column=1)
        self.var_ih3 = IntVar()
        ih3 = Checkbutton(master, text="dominant_denovo", variable=self.var_ih3,background=steel,fg='white',selectcolor=steel
                             ,highlightthickness=0,font=(font_type, txt_dim-2),activebackground= cobalt ,activeforeground='white' )
        ih3.grid(row=6 ,column=2)
        self.var_ih4 = IntVar()
        ih4 = Checkbutton(master, text="Xlinked", variable=self.var_ih4,background=steel,fg='white',selectcolor=steel
                             ,highlightthickness=0,font=(font_type, txt_dim-2),activebackground= cobalt ,activeforeground='white' )
        ih4.grid(row=7 ,column=0,columnspan=2)
        self.var_ih5 = IntVar()
        ih5 = Checkbutton(master, text="compound", variable=self.var_ih5,background=steel,fg='white',selectcolor=steel
                             ,highlightthickness=0,font=(font_type, txt_dim-2),activebackground= cobalt ,activeforeground='white' )
        ih5.grid(row=7 ,column=1,columnspan=2)
        
        #Force and Multisample tickmarks:
        Label(master, text="---------------------------------------------------",background=steel,fg='white').grid(row=8,column=0,columnspan=3)
        self.var_f = IntVar()
        self.f = Checkbutton(master, text="Force writing?", variable=self.var_f,background=steel,fg='white',selectcolor=steel
                             ,highlightthickness=0,font=(font_type, txt_dim-2),activebackground= cobalt ,activeforeground='white' )
        self.f.grid(row=9,column=0,columnspan=1,sticky='W')
       
        self.var_o = IntVar()
        self.o = Checkbutton(master, text="Multisample Input?", variable=self.var_o,background=steel,fg='white',selectcolor=steel
                             ,highlightthickness=0,font=(font_type, txt_dim-2),activebackground= cobalt ,activeforeground='white' )
        self.o.grid(row=9,column=2,columnspan=1,sticky='E') 
    
        #
        ##qsubname and jobname
        Label(master, text="---------------------------------------------------",background=steel,fg='white').grid(row=10,column=0,columnspan=3)

        self.qsub = Entry(master, background=cobalt,fg='white',font=(font_type, txt_dim-2))
        self.qsub.insert(0,'prioritize.sh')
        self.jname= Entry(master, background=cobalt,fg='white',font=(font_type, txt_dim-2))
        self.jname.insert(0,'prioritize_job')
        Label(master, text="Qsub name:", background=steel,fg='white',font=(font_type, txt_dim)).grid(row=11,column=0,columnspan=1,sticky='W')
        Label(master, text="Job name:", background=steel,fg='white',font=(font_type, txt_dim)).grid(row=12,column=0,columnspan=1,sticky='W')
        self.qsub.grid(row=11, column=2,columnspan=2,sticky='E')
        self.jname.grid(row=12, column=2,columnspan=2,sticky='E')
        
        return self.cfg_button
    def apply(self):
        NewWin = qSubDetails(self.master)
        if NewWin.result != None:
            self.result = {
                "force"     :self.var_f.get(),
                "multisample":self.var_o.get(),
                "outfolder" :self.outfolder,
                "familytype":self.var_value.get(),
                "qsubname"  :self.qsub.get(),
                "jobname"   :self.jname.get(),
                "inheritance":self.inheritance,
                "qopts"     :NewWin.result
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
            tkMessageBox.showwarning("Predict","Please Select one Output Folder")
            return 0
        #check qsubname 
        if len(self.qsub.get())<1 :
            tkMessageBox.showwarning("Predict","Qsub must contain text")
            return 0
        #check job name exists
        if len(self.jname.get())<1 :
            tkMessageBox.showwarning("Predict","Job name must contain text")
            return 0
        #Inheritance:
        or_val = self.var_ih5.get() + self.var_ih4.get() +self.var_ih3.get() +self.var_ih2.get() +self.var_ih1.get()
        if or_val == 0:
            tkMessageBox.showwarning("Predict","Choose at least one inheritance mode")
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
        w = Button(box, text="Run", width=10, command=self.ok,fg='white',font=(font_type, txt_dim)
                ,highlightbackground= steel,background=  steel,activebackground= cobalt ,activeforeground='white'        )
        w.pack(side=LEFT, padx=5, pady=5)
        w = Button(box, text="Back", width=10, command=self.cancel,fg='white',font=(font_type, txt_dim)
                ,highlightbackground= steel,background=  steel,activebackground= cobalt,activeforeground='white'       )
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
        Label(master, text="\nQsub Options:\n",background=orange,fg='white',font=(font_type, txt_dim)).grid(row=0,column=0,columnspan=3,sticky='W'+'E')
        
        #Select Queue
        Label(master, text="Select the qsub queue:",background=steel,fg='white',font=(font_type, txt_dim)).grid(row=2,column=0,columnspan=1,sticky=W+E)
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
        self.l_hh = Label(master, textvariable = (self.hh),width=lbl_width, background=steel,fg='white',font=(font_type, txt_dim))
        n_hh = Label(master, text = 'HH',width=lbl_width, background=steel,fg='white',font=(font_type, txt_dim-4))
        self.l_hh.grid(row=7,column=0,columnspan=1)
        n_hh.grid(row=6,column=0,columnspan=1)
        self.plus_hh = Button(master,text ='+',command=lambda : self.hh.set(self.plus1(self.hh,self.l_hh) ),
                                 background=steel,fg='white',font=(font_type, txt_dim),activebackground= orange,activeforeground='white'  )
        self.minus_hh = Button(master,text='-',command=lambda : self.hh.set(self.minus1(self.hh,self.l_hh) ),
                                 background=steel,fg='white',font=(font_type, txt_dim),activebackground=orange,activeforeground='white'  )
        self.plus_hh.grid(row=7,column=0,sticky='E')
        self.minus_hh.grid(row=7,column=0,sticky='W')
        
        self.mm = IntVar()
        self.mm.set(0)
        self.l_mm = Label(master, textvariable = (self.mm),width=lbl_width, background=steel,fg='white',font=(font_type, txt_dim))
        n_mm = Label(master, text = 'MM',width=lbl_width, background=steel,fg='white',font=(font_type, txt_dim-4))
        self.l_mm.grid(row=7,column=1,columnspan=1)
        n_mm.grid(row=6,column=1,columnspan=1)
        self.plus_mm = Button(master,text ='+',command=lambda : self.mm.set(self.plus1(self.mm,self.l_mm) ),
                                 background=steel,fg='white',font=(font_type, txt_dim),activebackground= orange,activeforeground='white'  )
        self.minus_mm = Button(master,text='-',command=lambda : self.mm.set(self.minus1(self.mm,self.l_mm) ),
                                 background=steel,fg='white',font=(font_type, txt_dim),activebackground=orange,activeforeground='white'  )
        self.plus_mm.grid(row=7,column=1,sticky='E')
        self.minus_mm.grid(row=7,column=1,sticky='W')
        
        self.ss = IntVar()
        self.ss.set(0)
        self.l_ss = Label(master, textvariable = (self.ss),width=lbl_width, background=steel,fg='white',font=(font_type, txt_dim))
        n_ss = Label(master, text = 'SS',width=lbl_width, background=steel,fg='white',font=(font_type, txt_dim-4))
        self.l_ss.grid(row=7,column=2,columnspan=1)
        n_ss.grid(row=6,column=2,columnspan=1)
        self.plus_ss = Button(master,text ='+',command=lambda : self.ss.set(self.plus1(self.ss,self.l_ss) ),
                                 background=steel,fg='white',font=(font_type, txt_dim),activebackground= orange,activeforeground='white'  )
        self.minus_ss = Button(master,text='-',command=lambda : self.ss.set(self.minus1(self.ss,self.l_ss) ),
                                 background=steel,fg='white',font=(font_type, txt_dim),activebackground=orange,activeforeground='white'  )
        self.plus_ss.grid(row=7,column=2,sticky='E')
        self.minus_ss.grid(row=7,column=2,sticky='W')
        
        #Mail Options
        self.var_a= IntVar()
        self.var_a.set(1)
        self.a = Checkbutton(master, text="Abort", variable=self.var_a,background=steel,fg='white',selectcolor=steel
                             ,highlightthickness=0,font=(font_type, txt_dim),activebackground= orange ,activeforeground='white'    )
        self.a.grid(row=5,column=0,columnspan=1,sticky=W)
       
        self.var_b = IntVar()
        self.var_b.set(1)
        self.b = Checkbutton(master, text="Begin", variable=self.var_b,background=steel,fg='white',selectcolor=steel,
                             highlightthickness=0,font=(font_type, txt_dim),activebackground= orange,activeforeground='white'     )
        self.b.grid(row=5,column=1,columnspan=1) 
        
        self.var_e = IntVar()
        self.var_e.set(1)
        self.e = Checkbutton(master, text="End", variable=self.var_e,background=steel,fg='white',selectcolor=steel,
                             highlightthickness=0,font=(font_type, txt_dim),activebackground= orange,activeforeground='white'     )
        self.e.grid(row=5,column=2,columnspan=1,sticky=E)
        
        self.var_custom = IntVar()
        self.var_custom.set(1)
        self.custom = Checkbutton(master, text="Custom Mail annotation", variable=self.var_custom ,background=steel,fg='white',selectcolor=steel,
                             highlightthickness=0,font=(font_type, txt_dim),activebackground= orange,activeforeground='white'     )
        self.custom .grid(row=4,column=1,columnspan=1,sticky=E)
        
        self.email = Entry(master, background=orange,fg='white',font=(font_type, txt_dim-2))
        self.email.insert(0,'xxx@crg.es')
        Label(master, text="Mail Addresss:", background=steel,fg='white',font=(font_type, txt_dim)).grid(row=8,column=0,columnspan=1,sticky=W)
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
        w = Button(box, text="Run", width=10, command=self.ok,fg='white',font=(font_type, txt_dim)
                ,highlightbackground= steel,background=  steel,activebackground= orange,activeforeground='white')
        w.pack(side=LEFT, padx=5, pady=5)
        w = Button(box, text="Back", width=10, command=self.cancel,fg='white',font=(font_type, txt_dim)
                ,highlightbackground= steel,background=  steel,activebackground= orange,activeforeground='white')
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