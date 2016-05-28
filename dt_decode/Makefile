name:=dowker_thistlethwaite
n2:=locate_in_pak

all: $(name) $(n2)

$(name): $(name).o
	cc $(name).o -o $(name) -Wall

$(name).o: $(name).c
	cc $(name).c -c -Wall

$(n2): $(n2).o
	cc $(n2).o -o $(n2) -Wall

$(n2).o: $(n2).c
	cc $(n2).c -c -Wall

clean:
	rm *.o $(name)

install:
	cp $(name) $(n2) /usr/local/bin

uninstall:
	rm /usr/local/bin/$(name) /usr/local/bin/($n2)
