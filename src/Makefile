#ifndef prefix
	prefix=/opt/amj
endif

FLAGS=-Wall -g -O3 -fPIC
CPPFLAGS=$(FLAGS)

AMJIMAGING_SRC=amjImaging.C
AMJIMAGING_OBJ=$(AMJIMAGING_OBJ:.C=.o)

CPP=g++ (CPPFLAGS)

all: libamjImaging.so

libamjImaging: $(AMJIMAGING_OBJ)
	$(CPP) -shared -Wl,-soname,$@ -o $@ $^

install: libamjImaging ../include/amjImaging.H ../include/amjPImaging.H
	mkdir -p $(prefix)/lib
	mkdir -p $(prefix)/include
	cp libamjFourier.so $(prefix)/lib
	cp ../include/amjImaging.H $(prefix)/include
	cp ../include/amjPImaging.H $(prefix)/include

uninstall:
	- rm $(prefix)/lib/amjImaging.so
	- rm $(prefix)/include/amjImaging.H
	- rm $(prefix)/include/amjPImaging.H
	- rmdir $(prefix)/lib
	- rmdir $(prefix)/include

clean:
	- rm $(AMJIMAGING_OBJ) libamjImaging.so
