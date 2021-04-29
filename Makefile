CC     =gcc
CFLAGS =-O2 -march=native -Werror -Wall
LIBS   =-lpthread -lm
TGTS   =newton

.PHONY : all clean cleanall

all : $(TGTS)

% : %.c
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

clean :
	@rm -f *.o

cleanall : clean
	rm -f $(TGTS)
