CC = g++
Objects = modify_genome_with_small_variants.o

#default target
all: modify_genome_with_small_variants
	@echo "compiliation done"

modify_genome_with_small_variants.o : modify_genome_with_small_variants.cpp
	$(CC) -c modify_genome_with_small_variants.cpp -o modify_genome_with_small_variants.o

bam_transform : $(Objects)
	$(CC) $(Objects) -o modify_genome_with_small_variants

#target deleting unwanted files
.PHONY: clean
clean:
	-rm -f *.o *~
