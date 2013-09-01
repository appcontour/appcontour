objects:=visible.o io.o extend.o build.o
sources:=$(objects:%.o=%.c)

all: visible morse2morse

visible: $(objects)
	cc $(objects) -o visible -Wall

$(objects): %.o: %.c visible.h
	cc $*.c -c -Wall

clean:
	rm *.o visible morse2morse

morseobjects:=morse2morse.o morse2morse_io.o morse2morse_extend.o morse2morse_collect.o
morsesources:=$(morseobjects:%.o=%.c)

morse2morse: $(morseobjects)
	cc $(morseobjects) -o morse2morse -Wall

$(morseobjects): %.o: %.c morse2morse.h
	cc $*.c -c -Wall

install:
	cp visible morse2morse /usr/local/bin

uninstall:
	rm /usr/local/bin/visible
	rm /usr/local/bin/morse2morse
