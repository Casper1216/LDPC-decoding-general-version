CC = g++
LIBNAME = rand_gen
OBJ = random_num_gen_msg.o channel_randSeq.o
LIBDIR = ../lib/
CFLAGS = -Wall -fPIC -O2 -lm -g

all : $(OBJ)
	$(CC) -shared -Wl,-soname,lib$(LIBNAME).so -o lib$(LIBNAME).so $(OBJ)
	[ -d $(LIBDIR) ] || mkdir -p $(LIBDIR)
	mv lib$(LIBNAME).so $(LIBDIR)
random_num_gen_msg.o : random_num_gen_msg.cpp random_num_gen_msg.h random_num_gen_unit.h
	$(CC) -c $(CFLAGS) $<
channel_randSeq.o : channel_randSeq.cpp channel_randSeq.h random_num_gen_unit.h
	$(CC) -c $(CFLAGS) $<

.PHONY : clean
clean :
	$(RM) $(OBJ) $(LIBDIR)lib$(LIBNAME).so
