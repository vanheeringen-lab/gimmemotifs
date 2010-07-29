## Makefile.am
bin_PROGRAMS = gadem
gadem_SOURCES = main.c mutation.c initial_population.c selection.c \
		sort.c alloc.c read_seq.c reverse_seq.c evalue_meme.c \
		pwm_score_distr.c scan_sites.c crossover.c mask_sites.c \
		construct_pwm.c	extend_alignment.c align_sites.c \
		top_kmers.c normalization.c score_subsequence.c transform_pwm.c \
		sample_wo_replacement.c output.c copy_pwm.c check_pwm_dist.c \
		consensus.c effect_seq_len.c check_convergence.c read_pwm0.c alloc.h \
		defines.h evalue_meme.h gadem.h random.h
                   
if LARGERINT 
   gadem_SOURCES += random_64bits.c
else
   gadem_SOURCES += random_32bits.c
endif

