#from tkinter import *
from Tkinter import *
import os

def show_entry_fields():
  print("First file: %s\nSecond file: %s" % (e1.get(), e2.get()))

def read_file1():
  return read_vec(e1.get())

def read_file2():
  return read_vec(e2.get())

def write_file1(filename, nrn_vec):
  write_file(filename, nrn_vec)
  
def write_file2(filename, nrn_vec):
  write_file(filename, nrn_vec)

dat_filenames=[]
def load_dat_filenames():
  print "finding dat files in ",os.getcwd()
  dir_str=os.listdir('.')
  for file in dir_str:
    if file[-4:]=='.dat':
      print file
      dat_filenames.append(file)

load_dat_filenames()

master = Tk()
Label(master, text="File 1").grid(row=0)
Label(master, text="File 2").grid(row=1)

e1 = Entry(master)
e2 = Entry(master)

e1.grid(row=0, column=1)
e2.grid(row=1, column=1)

Button(master, text='Quit', command=master.quit).grid(row=3, column=0, sticky=W, pady=4)
Button(master, text='Read File 1', command=read_file1).grid(row=4, column=1, sticky=W, pady=4)
Button(master, text='Write File 1', command=write_file1(e2.get(), nrn_vec).grid(row=5, column=1, sticky=W, pady=4)
Button(master, text='Read File 2', command=read_file2).grid(row=6, column=1, sticky=W, pady=4)
Button(master, text='Write File 2', command=write_file2).grid(row=7, column=1, sticky=W, pady=4)
Button(master, text='Show', command=show_entry_fields).grid(row=8, column=1, sticky=W, pady=4)



mainloop( )
