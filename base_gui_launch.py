##!/usr/bin/env python
import Tkinter as tk
import tkSimpleDialog as tkd
import os
import subprocess

curpath= os.path.realpath(__file__)
curpath = curpath.split('/')
curpath = '/'.join(curpath[:-1])
curpath += '/'
print curpath
python_path  = '/software/so/el6.3/PythonPackages-2.7.6/bin/python'
annotate_script = os.path.abspath(curpath+'Annotate/annotate.py')
predict_script = os.path.abspath(curpath+'Predict/prediction_pipeline.py')
prioritize_script = os.path.abspath(curpath+'Prioritize/prioritization_pipeline.py')

font_type ="Helvetica"


def predict_win():
    '''
    Opens the Prediction dialog box
    '''
    #Predict.config(state='disable')
    #Annotate.config(state='disable')
    #Prioritize.config(state='disable')
    #Quit.config(state='disable')

    NewWin = tkd.Predict_window(root)
    
    if NewWin.result != None:
        commandline  = """--infolder {0[infolder]} --outfolder {0[outfolder]} --qsubname {0[qsubname]} --config  {0[config]} --max_coverage {0[max_cov]}""".format(NewWin.result)
        commandline  += """ --namestart {0[namestart]} --namelength {0[namelength]} --firstreadextension {0[read1]} --secondreadextension  {0[read2]} --indelcaller {0[indel]}""".format(NewWin.result)
        commandline += """ --cpu {0[cpu]} --mem {0[mem]} """.format(NewWin.result)
        if NewWin.result['force']==1:
            commandline += " --fusevariants "
        commandline +=""" --qoptions " {0[qopts]}" """.format(NewWin.result) #The space before { is key for the command line to work
        if NewWin.result.get('sample_list',False):
               commandline += " --sample_list "+NewWin.result.get('sample_list',False)
               
        commandline = python_path + ' ' + predict_script + ' ' + commandline + '&'
        #print commandline         
        subprocess.Popen(commandline,shell=True)
        
        
    #Predict.config(state='normal')
    #Annotate.config(state='normal')
    #Prioritize.config(state='normal')
    #Quit.config(state='normal')
    #
    #


def annotate_win():
    '''
    Opens the Annotation dialog box
    '''
    #Predict.config(state='disable')
    #Annotate.config(state='disable')
    #Prioritize.config(state='disable')
    #Quit.config(state='disable')

    NewWin= tkd.Annotate_window(root)
    if NewWin.result != None:
        commandline = """-i {0[infile]} --variantType {0[vartype]} --sampleGenotypeMode  {0[gtype]} --geneDef {0[gdef]}""".format(NewWin.result)
        
        if NewWin.result['force']==1:
            commandline += " -f "
        if NewWin.result['only_gen']==1:
            commandline += " -o"

        commandline = python_path + ' ' + annotate_script + ' ' + commandline + '&'
        #print commandline
        subprocess.Popen(commandline,shell=True)
        
    #Predict.config(state='normal')
    #Annotate.config(state='normal')
    #Prioritize.config(state='normal')
    #Quit.config(state='normal')


def prioritize_win():
    '''
    Opens the Prioritization dialog box
    '''
    #
    #Predict.config(state='disable')
    #Annotate.config(state='disable')
    #Prioritize.config(state='disable')
    #Quit.config(state='disable')


    NewWin = tkd.Prioritize_window(root)
    if NewWin.result != None:
        commandline = """--outfolder {0[outfolder]} --qsubname {0[qsubname]} --jobname  {0[jobname]} --familytype {0[familytype]} --geneexclusion {0[gex_file]}""".format(NewWin.result)
        
        if NewWin.result['force']==1:
            commandline += " --force "
        if NewWin.result['multisample']==1:
            commandline += " --multisample"  
        for el in NewWin.result['inheritance']:
            commandline += ' --inheritance %s'%el
        if NewWin.result.get('config'):
            commandline += " --config %s"%NewWin.result["config"] 
        if NewWin.result.get('family'):
            commandline += " --family %s"%NewWin.result["family"]
        
        commandline +=""" --qoptions " {0[qopts]}" """.format(NewWin.result)#The space before { is key for the command line to work
        commandline = python_path + ' ' + prioritize_script + ' ' + commandline 
        #print commandline
        subprocess.Popen(commandline,shell=True)

    #Predict.config(state='normal')
    #Annotate.config(state='normal')
    #Prioritize.config(state='normal')
    #Quit.config(state='normal')
##############################
#### MAIN
##############################
#colors
back_color = "#7e3878"#aa00ff"
emerald ='#008a00'
mango ='#da532c'#e8980e'#'#fa6800'
cobalt ='#0050ef'
steel ="#1d1d1d"#647687'
green ='#60a917'
orange = '#fa6800'#'#cd3700''#f09609'
azzurro = '#1ba1f2'

root = tk.Tk()
w = root.winfo_screenwidth()
h = root.winfo_screenheight()
wh=500 #Window height
ww=500 #Window width
bw =15 #Button width
bh=4
txt_dim =20
root.title('')
root.geometry('%dx%d+%d+%d'%(ww,wh,w/2 - ww/2,h/2 -wh/2))
root.resizable(0,0)
root.configure(background=back_color)
##Buttons
w           = tk.Label(root, text="\t\n\tWhat do you want to do?  \t\n",background=back_color,fg='white',font=(font_type, txt_dim)
                       )
Predict     = tk.Button(root,text='Predict', command=predict_win,width=bw,height=bh,fg='white',font=(font_type, txt_dim)
                ,highlightbackground= emerald,relief='flat',background= emerald,activebackground= green     )
Annotate    = tk.Button(root,text='Annotate', command=annotate_win,width=bw,height=bh,fg='white',font=(font_type, txt_dim)
                ,highlightbackground= mango,relief='flat',background=mango,activebackground=  orange      )
Prioritize  = tk.Button(root,text='Prioritize', command=prioritize_win,width=bw,height=bh,fg='white',font=(font_type, txt_dim)
                ,highlightbackground= cobalt,relief='flat',background=  cobalt,activebackground= azzurro       )
Quit        = tk.Button(root,text='Quit', command=root.destroy,width=bw,height=bh,fg='white',font=(font_type,txt_dim)
                ,highlightbackground= steel,relief='flat',background=  steel,activebackground='#76608a'       )
w.grid(row=0,rowspan=3,columnspan=4,sticky='W'+'E'+'S'+'N')
Annotate.grid(row=4,column=1,sticky='E')
Predict.grid(row=3,column=1,sticky='E')
Prioritize.grid(row=3,column=2,sticky='W')
Quit.grid(row=4,column=2,sticky='W')

root.mainloop()

# http://msdn.microsoft.com/library/windows/apps/ff402557(v=vs.105).aspx

#PURPLE
#RGB 162 0 255
#A200FF
#
#MAGENTA
#RGB 255 0 151
#FF0097
#
#TEAL
#RGB 0 171 169
#00ABA9
#
#LIME
#RGB 140 191 38
#8CBF26
#
#BROWN
#RGB 160 80 0
#A05000
#
#PINK
#RGB 230 113 184
#E671B8
#
#ORANGE
#RGB 240 150 9
#F09609
#
#BLUE
#RGB 27 161 226
#1BA1E2
#
#RED
#RGB 229 20 0
#E51400