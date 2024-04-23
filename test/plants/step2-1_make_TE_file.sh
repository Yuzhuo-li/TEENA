# Ming-an Sun, 2024/4/20

## Download TE annotation file
wget https://rapdb.dna.affrc.go.jp/download/archive/irgsp1_repeat_unit.gff.gz
gunzip irgsp1_repeat_unit.gff.gz

## processing of rice TE annotation
cat irgsp1_repeat_unit.gff | \
# get useful columns
cut -f 1,4,5,9 | \
# get desired chrom
grep "^chr" | grep -v chrC | \
# only keep TEs
grep "DNA\|LINE\|SINE\|LTR" | \
# change the attribute field
sed 's/chr0/chr/; s/ID=.*Name=//;' | \
# get TEfam and TE name; refine the TE names
perl -ne 'my ($c, $s, $e, $x) = split; my ($name, $fam) = split(/\#/, $x); $fam .= ($fam =~ /\// ? ":" : "\/"); $name =~ s/_chr.*$//; print "$c\t$s\t$e\t$fam$name\n";' > irgsp1.TE.bed
