OBJS = Meld2Nexus.o General.o
CC = g++
DEBUG = -g
CFLAGS = -Wall -c -m64 -O3 -funroll-loops $(DEBUG)
LFLAGS = -Wall -m64 $(DEBUG)

Meld2Nexus: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o Meld2Nexus

Meld2Nexus.o: Meld2Nexus.cpp
	$(CC) $(CFLAGS) Meld2Nexus.cpp

General.o: General.cpp General.h
	$(CC) $(CFLAGS) General.cpp

clean:
	rm -rf *.o Meld2Nexus