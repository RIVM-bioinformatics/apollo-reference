bedtools bamtobed -i data/sam/SRR9007776__vs__GCA_002759435.3__external__CM105539.1.bam  | bedtools genomecov -i stdin -g data/cladegroup-cauris/GCA_002759435.3__external__CM105539.1.fa.fal  -bg | awk '{ l=($3-$2); for (i=0;i<l;i++) { print $4 } }' | sort -g | uniq -c | awk '{ if ($1>10) { print $2"\t"$1 } }' > data/sam/SRR9007776__vs__GCA_002759435.3__external__CM105539.1.bam.genomecov.tsv

