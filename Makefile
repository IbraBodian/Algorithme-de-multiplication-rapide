ifeq ($(shell uname -s),Darwin)
	SHA:=shasum -a 256
else
	SHA:=sha256sum
endif

check:
	([ -f fft.gp ] && ${MAKE} check-gp) || ([ -f fft.c ] && ${MAKE} check-c)

fft.gp.out: fft.gp
	gp -fq < $< > $@

fft.c.out: fft.c
	gcc -o fft $<
	./fft > $@

check-%: .%.sum fft.%.out
	$(SHA) -c $<
