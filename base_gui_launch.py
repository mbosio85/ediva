##!/usr/bin/env python
import Tkinter as tk
import tkSimpleDialog as tkd
import os
import subprocess
import ToolTip

curpath= os.path.realpath(__file__)
curpath = curpath.split('/')
curpath = '/'.join(curpath[:-1])
curpath += '/'
print curpath
python_path  = '/software/so/el6.3/PythonPackages-2.7.6/bin/python'
annotate_script = os.path.abspath(curpath+'Annotate/annotate.py')
predict_script = os.path.abspath(curpath+'Predict/prediction_pipeline.py')
prioritize_script = os.path.abspath(curpath+'Prioritize/prioritization_pipeline.py')
relaunch_script = os.path.abspath(curpath+'pipeline_control/qsubclass.py')
pipe_element_script = os.path.abspath(curpath+'pipeline_control/pipeline_element.py')
kill_script  = os.path.abspath(curpath+'pipeline_control/kill_jobs.py')

font_type ="Helvetica"


def predict_win():
    '''
    Opens the Prediction dialog box
    '''

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

def annotate_win():
    '''
    Opens the Annotation dialog box
    '''
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

def prioritize_win():
    '''
    Opens the Prioritization dialog box
    '''
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
        if NewWin.result.get('white_list'):
            commandline += " --white_list %s"%NewWin.result["white_list"]
        commandline +=""" --qoptions " {0[qopts]}" """.format(NewWin.result)#The space before { is key for the command line to work
        commandline = python_path + ' ' + prioritize_script + ' ' + commandline
        print commandline
        subprocess.Popen(commandline,shell=True)


def Resume_window():
    '''
    Re-launch a qsub or pipe job in case it crashed or to re-launch it
    It supports a resume option fully operational for the prioritize phase
    Less complete for the predict phase (just checks the step state without checking the needed outputs)
    '''
    NewWin = tkd.Resume_window(root)

    if NewWin.result != None:
        #print NewWin.result
        commandline =""" {0[filename]} "{0[qopts]}" """.format(NewWin.result) #The space before { is key for the command line to work
        commandline = python_path + ' ' + relaunch_script + ' ' + commandline + '&'
        print commandline

        subprocess.Popen(commandline,shell=True)


def Kill_window():
    '''
    Browse for a .qsub file which  hosts running jobs and kill them all
    '''
    NewWin = tkd.Kill_window(root)

    if NewWin.result != None:
        #print NewWin.result
        commandline =""" {0[filename]}  """.format(NewWin.result) #The space before { is key for the command line to work
        commandline = python_path + ' ' + kill_script + ' ' + commandline + '&'
        #print commandline
        subprocess.Popen(commandline,shell=True)
        #print commandline

##############################
#### MAIN
##############################
#colors
back_color = "#7e3878"#aa00ff"
emerald ='#008a00'
mango ='#da532c'#e8980e'#'#fa6800'
cobalt ='#2d89ef'#'#0050ef'
steel ="#1d1d1d"#647687'
green ='#60a917'
orange = '#fa6800'#'#cd3700''#f09609'
azzurro = '#1ba1f2'
text_color='white'#'#ffffff'

root = tk.Tk()
w = root.winfo_screenwidth()
h = root.winfo_screenheight()
wh=700 #Window height
ww=500 #Window width
bw =15 #Button width
bh=4
txt_dim =20
root.title('')
root.geometry('%dx%d+%d+%d'%(ww,wh,w/2 - ww/2,h/2 -wh/2))
root.resizable(0,0)
root.configure(background=back_color)
##Buttons
w           = tk.Label(root, text="\neDiVa\nWhat do you want to do?\n",background=back_color,fg=text_color,font=(font_type, txt_dim),
                      )
Predict     = tk.Button(root,text='Predict', command=predict_win,width=bw,height=bh,fg=text_color,font=(font_type, txt_dim)
                ,highlightbackground= emerald,relief='flat',background= emerald,activebackground= green     )
Annotate    = tk.Button(root,text='Annotate', command=annotate_win,width=bw,height=bh,fg=text_color,font=(font_type, txt_dim)
                ,highlightbackground= mango,relief='flat',background=mango,activebackground=  orange      )
Prioritize  = tk.Button(root,text='Prioritize', command=prioritize_win,width=bw,height=bh,fg=text_color,font=(font_type, txt_dim)
                ,highlightbackground=cobalt,relief='flat',background=cobalt,activebackground= azzurro       )
Quit        = tk.Button(root,text='Quit', command=root.destroy,width=bw,height=bh,fg=text_color,font=(font_type,txt_dim)
                ,highlightbackground= steel,relief='flat',background=  steel,activebackground='#76608a'       )
Rerun       = tk.Button(root,text='Rerun', command=Resume_window,width=bw,height=bh,fg=text_color,font=(font_type,txt_dim)
                ,highlightbackground= '#2C4566',relief='flat',background=  '#2C4566',activebackground='#76608a'       )
Kill        = tk.Button(root,text='Stop', command=Kill_window,width=bw,height=bh,fg=text_color,font=(font_type,txt_dim)
                ,highlightbackground= '#2c3e50',relief='flat',background=  '#2c3e50',activebackground='#76608a'       )
w.grid(row=0,rowspan=3,columnspan=3,sticky='W'+'E'+'S'+'N')
Annotate.grid(row=4,column=1,sticky='EW')
Predict.grid(row=3,column=1,sticky='EW')
Prioritize.grid(row=3,column=2,sticky='WE')
Quit.grid(row=4,column=2,sticky='WE')
Rerun.grid(row=5,column=1,sticky='EW')
Kill.grid(row=5,column=2,sticky='WE')

ToolTip.ToolTip(Predict, follow_mouse=1, text="eDiVa Varitant Prediction pipeline for sample-wise variant calling.")
ToolTip.ToolTip(Prioritize, follow_mouse=1, text="eDiVa Varitant filtering and ranking pipeline for trios-families.")
ToolTip.ToolTip(Annotate, follow_mouse=1, text="eDiVa Annotation for VCF files.")
ToolTip.ToolTip(Rerun, follow_mouse=1, text="Rerun a previously defined pipeline with new computing parameters.")
ToolTip.ToolTip(Kill, follow_mouse=1, text="Kill multiple jobs depending on a single pipeline.")



root.mainloop()
