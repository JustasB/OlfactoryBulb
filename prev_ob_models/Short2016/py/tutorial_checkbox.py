from Tkinter import *

master = Tk()

MODES = [
        ("Monochrome", "1"),
        ("Grayscale", "L"),
        ("True color", "RGB"),
        ("Color separation", "CMYK"),
    ]

v = StringVar()
v.set("L") # initialize

for text, mode in MODES:
        b = Radiobutton(master, text=text,
                        variable=v, value=mode)
        b.pack(anchor=W)

mainloop()
