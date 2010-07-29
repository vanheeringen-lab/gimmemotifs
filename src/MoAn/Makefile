SHELL=/bin/sh

CC=gcc
CFLAGS = -g -W  -O3
#CFLAGS = -g -W -Wall -Winline -O1

obj=obj/

ds_obj = $(obj)globals.o $(obj)ds.o $(obj)shallow.o $(obj)deep2.o \
	$(obj)helped.o $(obj)blind2.o

moan_obj=$(obj)anneal.o $(obj)sesa.o $(obj)score.o

wm_obj = $(obj)build.o $(obj)search.o $(obj)hittable.o $(obj)pssm.o \
	$(obj)timer.o $(obj)error.o $(obj)stringcon.o

bpr_obj = $(obj)kbs_Alphabet.o $(obj)kbs_Error.o $(obj)kbs_Math.o \
	$(obj)kbs_String.o $(obj)kbs_SuffixArray.o $(obj)kbs_SuffixArrayConstDStepAndPre.o


VPATH = ds_ssort:bpr

.PHONY: clean

moan: moan.c $(moan_obj) $(wm_obj) $(obj)ds_ssort.a
	$(CC) $(CFLAGS) -lm -o $@ $^


sesa: $(obj)sesa.o $(obj)hittable.o
	$(CC) $(CFLAGS) -lm -o $@ $^



# $(obj)sesa.o: sesa.c sesa.h
# 	$(CC) -c $(CFLAGS) $< -o $@

# $(obj)anneal.o: anneal.c anneal.h
# 	$(CC) -c $(CFLAGS) $< -o $@


# Archive containing the of the different algorithms
$(obj)ds_ssort.a: $(ds_obj)
	ar rc $@ $^

$(obj)bpr.a: $(bpr_obj)
	ar rc $@ $^

$(obj)wm.a: $(wm_obj) $(ds_obj) $(esa_obj)
	ar rc $@ $^


# Pattern rule for all objects files
obj/%.o: %.c *.h
	$(CC) -c $(CFLAGS) $< -o $@
 
# Clean up object files and archive files
clean: 
	rm -f *.o *.a $(obj)*.o $(obj)*.a \
	wm wm.exe \
	test1 test1.exe
	rm -rf build dist
