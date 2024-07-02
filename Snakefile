import os

if os.path.exists("Snakefile_versioned.sk"):
    include: "Snakefile_versioned.sk"

#This snakefile is the set of rules for convertion 
#of files and using them on motif analysing tools like
#meme-suite
#
#to start we need bigbed files and indexed genome file (fa.fai)
#for n in listofinputs

configfile: "config.yaml"
include: "config.smk"


rule all:
    input: 
        #expand("tomtom/{sample}/tomtom.{ext}" , sample=config['samples'], ext=["html", "tsv", "xml"])
        expand("tomtom/{sample}-{peak_filtering_method}/tomtom.{ext}" , sample=config['samples'], peak_filtering_method=config['peak_filtering_methods'], ext=["html", "tsv", "xml"])


#ruleorder: get_bb > get_bb_overlap

rule get_bb:
    input:
        '/path/to/file.bb'
    output:
        '{sample}-{peak_filering_method}.bb'
    shell: """
        cp {input} {output}
    """

rule get_bb_overlap:
    input:
        '/path/to/file.bb'
    output:
        '{sample}-{peak_filering_method}.bb'
    shell: """
        cp {input} {output}
    """

#convertion rules
rule bigbed_tobed:
   input:
       "{sample}.bb"
   output:
       "{sample}.bed"
   shell:
       """
       bigBedToBed {input} {output}
       """

rule bed_cutoff:
    input:
        "{sample}.bed"
    output:
        "{sample}.cut{K,\d+}.bed"
    shell: """
        bsort -k9,9gr {input} | bawk 'NR<={wildcards.K}'' > {output}
    """

rule bed_even_odd:
    input:
        "{sample}.bed"
    output:
        even="{sample}.even.bed",
        odd="{sample}.odd.bed"
    shell: """
        bawk 'NR%2==0 {{print > "{output.odd}" }} NR%2==1 {{print > "{output.even}" }}' {input}
    """


rule bed_tofasta:
   input:
       "{sample}.bed",
   output:
       "{sample}.fa"
   shell: """
      bedtools getfasta -fi {config[ref_genome]} -bed {wildcards.sample}.bed -fo {wildcards.sample}.fa
    """


#meme suite commands

rule meme_chip:
   input:
       "{sample}.fa"
   output:
       "meme_chip/{sample}/combined.meme",
       "meme_chip/{sample}/meme-chip.html",
       "meme_chip/{sample}/backround"
   shell: """
   mkdir -p meme_chip/{wildcards.sample}
   meme-chip -oc meme_chip/{wildcards.sample} -ccut 100 -dna -order 2 -meme-minw 6 -meme-maxw 15 -meme-mod zoops -meme-nmotifs {config[number_of_models]} -meme-searchsize 100000 \
   -centrimo-score 5.0 -centrimo-ethresh 10.0 {input} meme_chip/{wildcards.sample}
   """

rule name_fasta:
    input:
        "{sample}.fa"
    output:
        "{sample}.named.fa"
    shell:"""
        bawk '  BEGIN{{ count=0 }} \
                $0!~/^>/ {{print}} \
                $0~/^>/ {{\
                    count=count+1;\
                    gsub(">",">Name" count ";", $0);\
                    print;\
        }}' {input} > {output}
    """

rule shuffle_letters:
    input:
        "{sample}.fa"
    output:
        "{sample}.shuffle_letters.fa"
    shell: """
        fasta-shuffle-letters {input} -kmer 1 > {output}
    """


rule combine_pos_neg:
    input:
        pos="{sample}.fa",
        neg="{sample}.shuffle_letters.fa"
    output:
        "{sample}.pos_neg.fa"
    shell: """
    bawk '{{ gsub(">", ">positive;"); print $0 }}' {input.pos} > {output}
    bawk '{{ gsub(">", ">negative;"); print $0 }}' {input.neg} >> {output}
    """

rule fimo:
    input:
        fasta="{sample}.odd.named.pos_neg.fa",
        model="meme_chip/{sample}.even/combined.meme",
        background="meme_chip/{sample}.even/background"
    output:
        "fimo/{sample}.even_odd/fimo.tsv"
    shell: """
        mkdir -p $(dirname {output})
        fimo --parse-genomic-coord --verbosity 1 --oc $(dirname {output}) --bgfile {input.background} {input.model} {input.fasta}
    """
#fimo --parse-genomic-coord --verbosity 1 --oc meme_chip/GSE104479_PRJNA412810_rep1-idr.optimal/fimo_out_1 --bgfile meme_chip/GSE104479_PRJNA412810_rep1-idr.optimal/background --motif TGTTGCCCAGGCTGG meme_chip/GSE104479_PRJNA412810_rep1-idr.optimal/meme_out/meme.xml meme_chip/GSE104479_PRJNA412810_rep1-idr.optimal/GSE104479_PRJNA412810_rep1-idr.optimal.fa

rule fimo_named_sequence:
    input:
        "fimo/{sample}.even_odd/fimo.tsv"
    output:
        "fimo/{sample}.even_odd/fimo.named.tsv"
    shell:"""
        tr ";" "\t" < {input} > {output}
    """

rule tomtom:
   input:
       "meme_chip/{sample}/combined.meme",
   output:
       "tomtom/{sample}/tomtom.html",
       "tomtom/{sample}/tomtom.tsv",
       "tomtom/{sample}/tomtom.xml"
   shell: """
   mkdir -p tomtom/{wildcards.sample}
   tomtom -no-ssc -oc tomtom/{wildcards.sample} -min-overlap 5 -dist pearson -evalue -thresh 10.0 {config[dataset]} {input}
   """


