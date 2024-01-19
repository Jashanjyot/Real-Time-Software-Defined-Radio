# this makefile is intended for g++ on Linux

CC = g++
CFLAGS = -c -Wall -O3 -g
# CFLAGS = -c -Wall -O3 –mcpu=cortex-a71+crypto # on Rasp pi
LDFLAGS = -pthread
INCLUDE = -I../include/
SOURCES = project.cpp stereoProcessing.cpp iofunc.cpp filter.cpp  fourier.cpp  genfunc.cpp  logfunc.cpp settings.cpp rfFrontend.cpp mono.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = project

all: $(EXECUTABLE) clean

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDE) $< -o $@

clean:
	-rm $(OBJECTS)
