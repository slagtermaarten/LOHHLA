#!/bin/zsh

LOHHLA_dir='/home/m.slagter/libs/LOHHLA'
example_root="$LOHHLA_dir/example-file"
old_out="$example_root/correct-example-out"
new_out="$example_root/out"

setopt RE_MATCH_PCRE

files=($old_out/**/*.mpileup(:om))
files=($old_out/**/*.bam(:om))
files=($old_out/**/*.bed(:om))
files=($old_out/**/*.sam(:om))
lt $files

for f in $files; do
  paralog="${f:s/correct-example-out/out/}"
  tput setaf 4; echo "\n${f}";
  tput setaf 1;
  cmp --silent $f $paralog || echo "files are different"
  tput setaf 2; wc -l $f $paralog | \
    perl -ne 's/$ENV{old_out}|$ENV{o_dir}//; print($_);'
done


if [ 1 = 2 ]; then
  fti=example_tumor_sorted/example_tumor_sorted.chr6region.patient.reference.hlas.sam
  fti=example_BS_GL_sorted/example_BS_GL_sorted.chr6region.patient.reference.hlas.sam
  old_file="$old_out/$fti"
  new_file="$new_out/$fti"
  lt $old_file $new_file
  vimdiff $old_file $new_file
fi
