

./convert-sam-to-softclipped-bed.py data/sam/SRR9007776__vs__GCA_002759435.3__external__CM105539.1.bam > cladeV-softclipped.bed
./convert-sam-to-softclipped-bed.py filtered-TLEN-50-1000.bam > cladeV-softclipped-TLEN-50-1000.bed
./convert-sam-to-softclipped-bed.py filtered-TLEN-ext.bam > cladeV-softclipped-TLEN-ext.bed

fal=data/cladegroup-cauris/GCA_002759435.3__external__CM105539.1.fa.fal
bedtools genomecov -i cladeV-softclipped.bed -g $fal -bg > cladeV-softclipped.bg

bedtools genomecov -i cladeV-softclipped-TLEN-50-1000.bed -g $fal -bg > cladeV-softclipped-TLEN-50-1000.bg
bedtools genomecov -i cladeV-softclipped-TLEN-ext.bed -g $fal -bg > cladeV-softclipped-TLEN-ext.bg

bg2bw cladeV-softclipped-TLEN-50-1000.bg
bg2bw cladeV-softclipped-TLEN-ext.bg

#bcftools norm -a cladeV-0-5-1.vcf -o cladeV-0-5-1-norm.vcf
bcftools norm -a -m -any cladeV-0-5-1.vcf -o cladeV-0-5-1-norm.vcf
bcftools annotate -x FORMAT/GL,FORMAT/PL cladeV-0-5-1-norm.vcf -o cladeV-0-5-1-norm-fixed.vcf
bcftools filter -i 'QUAL > 1' cladeV-0-5-1-norm-fixed.vcf -Ob -o cladeV-0-5-1-norm-fixed-filt.vcf
bgzip cladeV-0-5-1-norm-fixed-filt.vcf
bcftools index cladeV-0-5-1-norm-fixed-filt.vcf.gz

function bg2bw() {
    fal=data/cladegroup-cauris/GCA_002759435.3__external__CM105539.1.fa.fal
    bg=$1
    bw=$(echo $bg | sed 's/\.bg$/.bw/')
    ~/software/ucsc/bedGraphToBigWig $bg $fal $bw
    }


cat cladeV-softclipped.bg | awk '{ if ($4>10) { print $0 } }' > cladeV-softclipped-gte10.bg
~/software/ucsc/bedGraphToBigWig cladeV-softclipped.bg $fal cladeV-softclipped.bw
~/software/ucsc/bedGraphToBigWig cladeV-softclipped-gte10.bg $fal cladeV-softclipped-gte10.bw

rm -f cladeV-softclipped.bg
rm -f cladeV-softclipped-gte10.bg


fal=data/cladegroup-cauris/GCA_002759435.3__external__CM105539.1.fa.fal
bedtools bamtobed -i data/sam/SRR9007776__vs__GCA_002759435.3__external__CM105539.1.bam  | bedtools genomecov -i stdin -g $fal -bg > cladeV-coverage.bg
~/software/ucsc/bedGraphToBigWig cladeV-coverage.bg $fal cladeV-coverage.bw
rm -f cladeV-coverage.bg


# --binSize 1 == slooooooow
bigwigCompare -b1 cladeV-softclipped.bw -b2 cladeV-coverage.bw --operation subtract --binSize 1 --outFileFormat bedgraph -o verschil.bg


bedtools genomecov -i cladeV-softclipped.bed -g $fal -bg > cladeV-softclipped.bg
bedtools bamtobed -i data/sam/SRR9007776__vs__GCA_002759435.3__external__CM105539.1.bam  | bedtools genomecov -i stdin -g $fal -bg > cladeV-coverage.bg
bedtools unionbedg -i cladeV-coverage.bg cladeV-softclipped.bg | awk '{ if ($5>=$4) { print $0 } }' | cut -f 1-3,5 > cladeV-softclipped-filtered.bg
~/software/ucsc/bedGraphToBigWig cladeV-softclipped-filtered.bg $fal cladeV-softclipped-filtered.bw

grep "^IS" data/sam/SRR9007776__vs__GCA_002759435.3__external__CM105539.1.bam.stats  | sed 1d | sed '$d' | awk '{ for (i=0;i<$3;i++) { print $2 } }' | ~/toolkit/nsmmmmsd.py
2290208	673604278	2	7997	279.0	294.124	163.2

bam=../sam/SRR9007776__vs__GCA_002759435.3__external__CM105539.1.bam
freebayes -f GCA_002759435.3__external__CM105539.1.fa -p 2 --haplotype-length 3 --min-repeat-size 5 --min-repeat-entropy 1 $bam > cladeV-3-5-1.vcf
# see the freebayes manual
freebayes -f GCA_002759435.3__external__CM105539.1.fa -p 2 --haplotype-length -1 --min-repeat-size 5 --min-repeat-entropy 1 $bam > cladeV-M1-5-1.vcf


bam=../../filtered-TLEN-50-1000.bam
freebayes -f GCA_002759435.3__external__CM105539.1.fa -p 2 --haplotype-length 3 --min-repeat-size 5 --min-repeat-entropy 1 $bam | \
    bcftools norm -a -m -any | bcftools annotate -x FORMAT/GL,FORMAT/PL | \
    bcftools filter -i 'QUAL >= 5' -Ob -o cladeV-3-5-1-TLEN-flt.vcf &
bam=../../filtered-TLEN-ext.bam
freebayes -f GCA_002759435.3__external__CM105539.1.fa -p 2 --haplotype-length 3 --min-repeat-size 5 --min-repeat-entropy 1 $bam | \
    bcftools norm -a -m -any | bcftools annotate -x FORMAT/GL,FORMAT/PL | \
    bcftools filter -i 'QUAL >= 5' -Ob -o cladeV-3-5-1-TLEN-ext.vcf &
wait

bgzip cladeV-3-5-1-TLEN-flt.vcf
bgzip cladeV-3-5-1-TLEN-ext.vcf
tabix -p vcf cladeV-3-5-1-TLEN-flt.vcf.gz
tabix -p vcf cladeV-3-5-1-TLEN-ext.vcf.gz
bcftools isec -p output_dir -n~01 cladeV-0-5-1-norm-fixed-filt.vcf.gz cladeV-0-5-1-norm-fixed-ext.vcf.gz

samtools view -e 'tlen*tlen >= 100*100 && tlen*tlen <= 1000*1000' -b data/sam/SRR9007776__vs__GCA_002759435.3__external__CM105539.1.bam > filtered-TLEN-100-1000.bam 
samtools view -e 'tlen*tlen < 100*100 || tlen*tlen > 1000*1000' -b data/sam/SRR9007776__vs__GCA_002759435.3__external__CM105539.1.bam > filtered-TLEN-ext.bam 
samtools index filtered-TLEN-50-1000.bam
samtools index filtered-TLEN-ext.bam


./convert-sam-to-softclipped-bed.py data/sam/SRR9007776__vs__GCA_002759435.3__external__CM105539.1.bam > cladeV-softclipped.bed
./convert-sam-to-softclipped-bed.py filtered-TLEN-50-1000.bam > cladeV-softclipped-TLEN-50-1000.bed
./convert-sam-to-softclipped-bed.py filtered-TLEN-ext.bam > cladeV-softclipped-TLEN-ext.bed



fal=data/cladegroup-cauris/GCA_002759435.3__external__CM105539.1.fa.fal
bedtools bamtobed -i filtered-TLEN-50-1000.bam  | bedtools genomecov -i stdin -g $fal -bg > cladeV-coverage-TLEN-50-1000.bg
~/software/ucsc/bedGraphToBigWig cladeV-coverage-TLEN-50-1000.bg $fal cladeV-coverage-TLEN-50-1000.bw
rm -f cladeV-coverage-TLEN-50-1000.bg


bedtools genomecov -i cladeV-softclipped-TLEN-50-1000.bed -g $fal -bg > cladeV-softclipped-TLEN-50-1000.bg
bedtools bamtobed -i filtered-TLEN-50-1000.bam  | bedtools genomecov -i stdin -g $fal -bg > cladeV-coverage-TLEN-50-1000.bg

bedtools unionbedg -i cladeV-coverage-TLEN-50-1000.bg cladeV-softclipped-TLEN-50-1000.bg | awk '{ if ($5>=($4+5)) { print $0 } }' | cut -f 1-3,5 > cladeV-softclipped-filtered.bg
~/software/ucsc/bedGraphToBigWig cladeV-softclipped-filtered.bg $fal cladeV-softclipped-filtered.bw


bcftools view cladeV-0-5-1-norm-fixed-filt.vcf.gz -H | awk '{ print $1"\t"$2-1"\t"$2"\t"$6 }' | bedtools intersect -a stdin -v -b <(~/software/ucsc/bigWigToBedGraph ../../cladeV-softclipped-filtered.bw  stdout)  | awk '{ print $4 }' | ~/toolkit/nsmmmmsd.py  --header | csvlook -t -I
| n      | sum          | min     | max      | median  | mean     | sd       |
| ------ | ------------ | ------- | -------- | ------- | -------- | -------- |
| 307104 | 463167852.64 | 1.01401 | 190036.0 | 1168.14 | 1508.179 | 2207.032 |

bcftools view cladeV-0-5-1-norm-fixed-filt.vcf.gz -H | awk '{ print $1"\t"$2-1"\t"$2"\t"$6 }' | bedtools intersect -a stdin -b <(~/software/ucsc/bigWigToBedGraph ../../cladeV-softclipped-filtered.bw  stdout)  | awk '{ print $4 }' | ~/toolkit/nsmmmmsd.py  --header | csvlook -t -I
| n     | sum         | min     | max     | median  | mean     | sd        |
| ----- | ----------- | ------- | ------- | ------- | -------- | --------- |
| 13789 | 9311980.433 | 1.03059 | 41176.6 | 371.842 | 675.3195 | 1401.4081 |



# dustmasking as blacklisting doesn't work as expected ...
bcftools view cladeV-0-5-1-norm-fixed-filt.vcf.gz -H | awk '{ print $1"\t"$2-1"\t"$2"\t"$6 }' | bedtools intersect -a stdin -b <(cat GCA_002759435.3__external__CM105539.1.dusted.bed | awk '{ if ($4>=40) { print $0 } }') | \
    awk '{ print $4 }' | ~/toolkit/nsmmmmsd.py  --header | csvlook -t -I
| n    | sum          | min     | max      | median  | mean      | sd        |
| ---- | ------------ | ------- | -------- | ------- | --------- | --------- |
| 9601 | 17375196.396 | 1.24853 | 190036.0 | 714.184 | 1809.7278 | 6858.8089 |


# and now lest "larger deletetions"
bedtools subtract -a <(bedtools bamtobed -i cladeV.bam) -b <(bedtools bamtobed -splitD -i cladeV.bam | awk '{ if ($3-$2>20) { print $0 } }') | awk '{ if ($3-$2 >= 20) { print $0"\t"$3-$2 } }' > cladeV.gaps.bed

# works a bit, but not a as sharp as softclipped-filtered
bcftools view cladeV-0-5-1-norm-fixed-filt.vcf.gz -H | awk '{ print $1"\t"$2-1"\t"$2"\t"$6 }' | \
    bedtools intersect -a stdin -b <(cat cladeV.gaps.bed | awk '{ if ($4>=40) { print $0 } }') | \
    awk '{ print $4 }' | ~/toolkit/nsmmmmsd.py  --header | csvlook -t -I
| n    | sum         | min     | max     | median  | mean      | sd        |
| ---- | ----------- | ------- | ------- | ------- | --------- | --------- |
| 2732 | 2971571.545 | 1.21459 | 35095.4 | 400.274 | 1087.6909 | 3523.6658 |


cd ~/gitrepos/my-appollo-mapping/data/cladegroup-cauris
( echo "main"; bcftools view cladeV-0-5-1-norm-fixed-filt.vcf.gz -H | awk '{ print $1"\t"$2-1"\t"$2"\t"$6 }' | bedtools intersect -a stdin -v -b <(~/software/ucsc/bigWigToBedGraph ../../cladeV-softclipped-filtered.bw  stdout)  | awk '{ print $4 }'; echo "softclipped"; bcftools view cladeV-0-5-1-norm-fixed-filt.vcf.gz -H | awk '{ print $1"\t"$2-1"\t"$2"\t"$6 }' | bedtools intersect -a stdin  -b <(~/software/ucsc/bigWigToBedGraph ../../cladeV-softclipped-filtered.bw  stdout)  | awk '{ print $4 }' )  | ~/toolkit/plot-binned-histogram.py --max 10000 --multi-headered



bgzip cladeV.vcf
tabix -p vcf cladeV.vcf.gz
bcftools isec -p output_dir-cladeV-ref-11 -n~11 cladeV-0-5-1-norm-fixed-filt.vcf.gz cladeV.vcf.gz
bcftools isec -p output_dir-cladeV-ref-10 -n~10 cladeV-0-5-1-norm-fixed-filt.vcf.gz cladeV.vcf.gz


