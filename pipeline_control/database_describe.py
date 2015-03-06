#!/usr/bin/env python
#!/usr/bin/env python
import os
import subprocess
#import mysql.connector
import MySQLdb
import re
import ntpath
import struct
import hashlib
import argparse


username    = "edivacrg"
#database    = "eDiVa_public_omics"
dbhost      = "mysqlsrv-ediva.linux.crg.es"
passw       = "FD5KrT3q"
 
db = MySQLdb.connect(host=dbhost, # your host, usually localhost
user=username, # your username
passwd=passw, # your password
) # name of the data base

cur = db.cursor()

# List of databases
sql = "show databases;"
cur.execute(sql)

dbs = list()


for row in cur:
    dbs.append(row[0])

del dbs[0]
#print tables for each one
colcount=0
dbcount=0
tablecount =0
for i in dbs:
    sql = "show tables in %s;"%i
    cur.execute(sql)
    dbcount +=1
    for row in cur:
       # print "%s : %s"%(i,row[0])
        sql2 = "select * from %s.%s limit 1000"%(i,row[0])  #Select all items from a table
        cur2 = db.cursor()
        cur2.execute(sql2)
        cur3 =db.cursor()
        a =  cur2.fetchall()
        
        if len(a)>=1:
            
            print "%s : %s"%(i,row[0])
            sql3  =( "select `COLUMN_NAME` FROM `INFORMATION_SCHEMA`.`COLUMNS` "+
                "WHERE `TABLE_SCHEMA`='%s' AND `TABLE_NAME`='%s';")%(i,row[0])
            cur3.execute(sql3)
            b = cur3.fetchall()
            print("\t %d columns\n \t %s \n"%(len(b),str(b)))
            if len(b)>1:
                tablecount +=1
                colcount += len(b)
            
        
        cur3.close()  
        cur2.close()

print "%d Databases \n \t%d Tables \n \t\t%d Columns"%(dbcount,tablecount,colcount)

 #       
 #SELECT `COLUMN_NAME`
 #   -> FROM `INFORMATION_SCHEMA`.`COLUMNS`
 #   -> WHERE `TABLE_SCHEMA`='eDiVa_innoDB'
 #   -> AND `TABLE_NAME`='Table_sample';

cur.close()

db.close()
