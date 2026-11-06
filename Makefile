CC = cc
CFLAGS = -std=c11 -O2 -Wall -Wextra
CFLAGS += $(shell pkg-config --cflags libcurl cairo freetype2 libtiff-4 libgeotiff proj)
LIBS = $(shell pkg-config --libs libcurl cairo freetype2 libtiff-4 libgeotiff proj) -lm

.PHONY: all clean maps

all: render_altadena

render_altadena: render_altadena.c
	$(CC) $(CFLAGS) -o $@ render_altadena.c $(LIBS)

maps: render_altadena
	./render_altadena

clean:
	rm -f render_altadena
