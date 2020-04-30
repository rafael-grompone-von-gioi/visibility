STRICT=#-fPIC -ansi -Wall -Wextra -Werror
IIO=-DNDEBUG -std=c99 -lpng -ltiff -ljpeg
OPT= -O3

visibility: visibility.c iio.c
	$(CC) $(OPT) $(STRICT) -o $@ $^ $(IIO) -lm

test: visibility
	./visibility 500 data/image*

clean:
	rm -f visibility
	rm -f 000.png 001.png 002.png 003.png 004.png
	rm -f 005.png 006.png 007.png 008.png 009.png
