APPNAME = 'nora'
VERSION = '1.0.0'
CC = gcc

PREFIX=$(PWD)
COMPILED= \
bfmt72s_v3 \
nss2v_v4 \
mhseol_v3 \
gt_v3 \

SCRIPTS= \
bfmt7viewer.pl \
cfq2fa.pl \
clean_gfa.pl \
extract_gfa.pl \
fa2gfaS.pl \
fa2idfa_v2.pl \
get_blastout_in_xml.pl \
get_fasta_stats.pl \
get_longest_20x_fa.pl \
idfa2samheader.pl \
nora.pl \
ohf2fa_v2.pl \
vertical2cfq_and_blacklist.pl \
vertical2fa.pl \


all: $(COMPILED)

bfmt72s_v3: bfmt72s_v3.c
	$(CC) -Wall -O3 -g -o $@ $<

nss2v_v4: nss2v_v4.c
	$(CC) -Wall -O3 -g -o $@ $<

mhseol_v3: mhseol_v3.c b_heikki.dynamic.h mymurmurhash3.h myqsort.h
	$(CC) -Wall -O3 -g -fopenmp -o $@ $<

gt_v3: gt_v3.c
	$(CC) -Wall -O3 -g -fopenmp -o $@ $<


install: $(COMPILED) $(SCRIPTS)
	chmod 766 $^
	cp $^ $(PREFIX)/bin/

clean:
	rm $(COMPILED)

CCODES= \
bfmt72s_v3.c \
nss2v_v4.c \
mhseol_v3.c \
gt_v3.c \

HEADERS= \
b_heikki.dynamic.h \
mymurmurhash3.h \
myqsort.h \

OTHERS= \
Makefile \
vertical.vim \

tar:
	touch $(CCODES) $(HEADERS) $(SCRIPTS) $(OTHERS)
	tar -cvzf nora.$(VERSION).tar.gz $(CCODES) $(HEADERS) $(SCRIPTS) $(OTHERS)

