CC=gcc

CFLAGS=-c -Wall -fPIC -lm -I./source_mixed

SOURCES=./source_mixed/IFnet.c
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=IFnet

all:    $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:: 
	rm -f *.o IFnet 


