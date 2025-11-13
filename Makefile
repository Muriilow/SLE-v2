CC = gcc

CFLAGS = -O3 -march=native -mavx -fopt-info-vec -DLIKWID_PERFMON -I${LIKWID_INCLUDE}
LFLAGS = -lm -L${LIKWID_LIB} -llikwid

PROG = cgSolver
MODULES = utils	sislin utils $(PROG)
OBJS = $(addsuffix .o,$(MODULES))
SRCS = $(addsuffix .c,$(MODULES)) $(addsuffix .h,$(MODULES))

# Lista de arquivos para distribuição
DISTFILES = *.c *.h Makefile LEIAME
DISTDIR = gvso24-mpb24

.PHONY: clean purge dist all

%.o: %.c %.h utils.h
	$(CC) -c $(CFLAGS) $<

$(PROG):  $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)

clean:
	@echo "Limpando sujeira ....."
	@rm -rf core *~ *.bak

purge: clean
	@echo "Fazendo a faxina ....."
	@rm -f a.out *.o $(PROG)


dist: purge
	@echo "Gerando arquivo de distribuição ($(DISTDIR).tgz) ..."
	@ln -s . $(DISTDIR)
	@tar -chzvf $(DISTDIR).tgz $(addprefix ./$(DISTDIR)/, $(DISTFILES))
	@rm -f $(DISTDIR)
