#ifndef prefix
	prefix=/opt/amj
#endif

FLAGS=-Wall -g -O3 -fPIC -I/opt/amj/include
CPPFLAGS=$(FLAGS)

AMJIMAGING_SRC=amjImaging.C
AMJIMAGING_OBJ=$(AMJIMAGING_SRC:.C=.o)

CPP=g++ $(CPPFLAGS)

all: libamjImaging.so

libamjImaging.so: $(AMJIMAGING_OBJ)
	$(CPP) -shared -Wl,-soname,$@ -o $@ $^

install: libamjImaging.so ../include/amjImaging.H ../include/amjPImaging.H
	mkdir -p $(prefix)/lib
	mkdir -p $(prefix)/include
	cp libamjImaging.so $(prefix)/lib
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
