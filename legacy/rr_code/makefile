CC = gcc
CFLAGS = -O2 
CLIBS = -lm 

OBJS = sp_one.o white_bear.o

OBJS2 = sp_one.o white_bear_mark_II.o

OBJS3 = sp_one.o rosenfeld.o

OBJS4 = sp_one.o white_bear_q3.o

OBJS5 = sp_one.o white_bear_mark_II_q3.o

OBJS6 = sp_one.o rosenfeld_q3.o


all: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(CLIBS) -o sponewb

wb2: $(OBJS2)
	$(CC) $(CFLAGS) $(OBJS2) $(CLIBS) -o sponewb2

rf: $(OBJS3)
	$(CC) $(CFLAGS) $(OBJS3) $(CLIBS) -o sponerf

wbq3: $(OBJS4)
	$(CC) $(CFLAGS) $(OBJS4) $(CLIBS) -o sponewbq3

wb2q3: $(OBJS5)
	$(CC) $(CFLAGS) $(OBJS5) $(CLIBS) -o sponewb2q3

rfq3: $(OBJS6)
	$(CC) $(CFLAGS) $(OBJS6) $(CLIBS) -o sponerfq3

clean : 
	rm -f *.o core 
