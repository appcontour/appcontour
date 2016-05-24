name:=dowker_thistlethwaite

all: $(name)

$(name): $(name).o
	cc $(name).o -o $(name) -Wall

$(name).o: $(name).c
	cc $(name).c -c -Wall

clean:
	rm *.o $(name)

install:
	cp $(name) /usr/local/bin

uninstall:
	rm /usr/local/bin/$(name)
