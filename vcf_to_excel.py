import argparse
import os
import pprint

try:
    import xlsxwriter
    import xlrd
    writeXLS = True
    print 'yes XLS module imported'
except:
    writeXLS = False
    print 'No XLS module imported'
    #!/usr/bin/env python
parser = argparse.ArgumentParser(description = 'Create a file that will be run in the cluster to fully run the prioritization in eDiVa.')
parser.add_argument('-tab',  dest='tabs', required=False, action = 'append', help="""choose a tab vcf file'""")
parser.add_argument('-out',  dest='outname', required=False)

args = parser.parse_args()

try :
    outname = args.outname + '.xlsx'
    print outname
except:
    outname = 'out.xlsx'
sp = pprint.PrettyPrinter()

print args.tabs
print('Excel writing')
if writeXLS == True:
    excel_name = outname
    print outname

    xls = xlsxwriter.Workbook(excel_name)
    for tabname in args.tabs:
        dir_name = tabname.split('/')[0]
        print dir_name
        worksheet = xls.add_worksheet(dir_name)
        row = 0
        tabfile = open(os.path.abspath(tabname),'r')
        # read line by line and transform to xls
        for line in tabfile:
            data = line.split(',')
            worksheet.write_row(row, 0, data)
            #sp.pprint( data)
            row += 1
        tabfile.close()
    xls.close()
