#####################################################################
#                                                                   #
#               Make file for PinMoka                               #
#                                                                   #
#                         tiagobscastro@gmail.com                   #
#####################################################################

# executable name
PROG = bin/PinMoka.x
PREFIXDIR = /home/$(USER)
MAIN = PinMoka.cpp power2D.cpp cMZhao.cpp cMBhattacharya.cpp util.cpp

# .cpp lib internal in MOKA
SOURCES = ../Moka/utilities.cpp \
          ../Moka/cosmology.cpp \
          ../Moka/halo.cpp \
          ../Moka/nfwHalo.cpp \
          ../Moka/nfwLens.cpp \
          ../Moka/hernq_Halo.cpp \
          ../Moka/hernq_Lens.cpp \
          ../Moka/jaffe_Halo.cpp \
          ../Moka/jaffe_Lens.cpp \
          ../Moka/sisHalo.cpp \
          ../Moka/sisLens.cpp


# gsl, cfitsio, CCfits, fftw
LIBS = -L$(PREFIXDIR)/lib/ -lgslcblas -lgsl \
       -L$(PREFIXDIR)/lib/ -lm \
       -L$(PREFIXDIR)/lib/ -lCCfits -lcfitsio \
       -L$(PREFIXDIR)/lib/ -lfftw3_mpi -lfftw3

# gsl, cfitsio, CCfits, fftw
ALLFLAGS = -I$(PREFIXDIR)/include/ \
           -I$(PREFIXDIR)/include/gsl \
           -I$(PREFIXDIR)/include/CCfits

#
DEBUG = -g
# compiler -Wall # use this to optimize -O2
CC = gcc -std=c++11 -lstdc++ -O3
#
RM = rm -f -r
#
OBJ = $(SOURCES:.cpp=.o)
#

CFLAGS=$(ALLFLAGS) $(DEBUG)
#
CLEAR = clear

default: moka
moka:
	$(CC) -c $(SOURCES) $(CFLAGS)
	ar r libmoka.a *.o
	$(CC) $(MAIN) -L. -lmoka $(CFLAGS) -o $(PROG) $(LIBS)

main:
	$(CC) $(MAIN) -L. -lmoka $(CFLAGS) -o $(PROG) $(LIBS)

clean:
	$(RM) $(PROG) $(OBJ) *.o *a

lib:
	$(CC) -c $(SOURCES) $(CFLAGS)
	ar r libmoka.a $(OBJ)

clall:
	$(RM) $(PROG) $(OBJ) libmoka.a *~ *.o

clall0:
	$(RM) $(PROG) $(OBJ) libmoka.a *~ html *.o

new: clean default
