REALGLOMS_FNAME = 'realgloms.txt'
glomRealCoords = []

def loadGloms():
    # load gloms positions
    f = open(REALGLOMS_FNAME)
    line = f.readline()
    while line:
        tk = line.split()
        # every line has N glom X Y of i-glomerolous
        glomRealCoords.append([ float(tk[0]), float(tk[1]), float(tk[2]) ])
        line = f.readline()

loadGloms()



    

                    
                    
            
        
    
