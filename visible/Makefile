objects:=visible.o io.o extend.o build.o
sources:=$(objects:%.o=%.c)

visible: $(objects)
	cc $(objects) -o visible -Wall

%.o: %.c visible.h
	cc $*.c -c -Wall

clean:
	rm *.o visible
