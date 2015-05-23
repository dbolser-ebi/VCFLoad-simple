#!/bin/env perl

use strict;
use warnings;

## The list of assumptions made by this script (for speed):

## 0) No existing data in the target variation database.

## 1) One source, ID given below.

## 2) One population, ID given below.

## 3) One line per-locus, per-variation. i.e. no multiply mapped
##    variations and no sets of 'overlapping' VCF files.

## 4) Variations have IDs (we don't generate them!).

## 5) GT is the first field of FORMAT (and is unphased). This
##    positional requirement is just for speed, and could be coded to
##    match the format more cleanly.

## 6) A static hash mapping sequence names to internal IDs can be
##    provided as a simple TSV (see MAP below).



## Spew some random debugging
my $debug = 0;

## One variation can only ever be linked to one source!
my $source_id = 1;

## This /could/ be a hash linking individuals to populations...
my $population_id = 1;

## These are those currently used by import_vcf.pl with IDs from the
## schema attribs table.
my %class_attrib_id =
    qw( SNV                  2
        substitution         5
        insertion           10
        deletion            12
        sequence_alteration 18
    );



## We now use a file to map 'chromosome names' in the vcf to
## seq_region_ids in the variation database. This is because there can
## be all kinds of weird and non-one-to-one mappings between the two.

## Use something like this...
## mysql-prod-1 triticum_aestivum_core_27_80_2 -Ne '
##   SELECT name, seq_region_id FROM seq_region
##   INNER JOIN seq_region_attrib WHERE attrib_type_id = 6
## ' > seq_region_id-triticum_aestivum_core_27_80_2.mapping

my %seq_region_id;

#my $mapping_file = 'seq_region_id-hordeum_vulgare_core_25_78_2.mapping';
#my $mapping_file = 'seq_region_id-solanum_lycopersicum_core_27_80_2.mapping';
my $mapping_file = 'seq_region_id-triticum_aestivum_core_27_80_2.mapping';

open MAP, '<', $mapping_file
    or die "failed to open mapping file $mapping_file: $?\n";

while(<MAP>){
    chomp;

    my ($name, $seq_region_id) = split/\t/;

    warn "Multiple mappings for '$name'\n"
        if exists $seq_region_id{$name};

    $seq_region_id{$name} = $seq_region_id;
}

warn "got ", scalar keys %seq_region_id, " mappings\n";



## GO

die "pass a vcf at least\n"
  unless @ARGV;

my $outdir = './MySQL';



## The files we're writing
my @file =
    qw( individual
        allele_code
        genotype_code
        population_allele
        population_genotype
        individual_genotype_x
        individual_genotype_y
        variation
        variation_feature
    );

## Open all files for writing
my %file;

for my $file (@file){
    local *FILE;
    open( FILE, '>', "$outdir/$file.tsv")
      or die "failed to open file $outdir/$file.tsv : $!\n";
    $file{$file} = *FILE;
}





## OK

my $variation_id;

## If it were my schema, I'd define this per variation using a
## compound primary key. It's much more intuitive that way. But then
## again...
my $allele_id;
my $genotype_id;

## The maps we build
my %allele_code;
my %genotype_code;

while(<>){
    next if /^##/;
    chomp;

    ## Catch the individuals here, you'll need to add columns to this
    ## table as appropriate...
    if (/^#/){
        my @cols = split/\t/;
        my $i = 1;
        print { $file{individual} } $i++, "\t$_\t3\n"
            for @cols[9..$#cols];
        next;
    }

    $variation_id++;

    my ($chr, $pos, $id, $ref, $alt,
        ## None of these are (currently) used here
        undef, # qual
        undef, # filter
        undef, # info
        undef, # format
        ## But we need these!
        @sample
        ) = split/\t/;

    ## Debugging
    print
        join("\t", $chr, $pos, $id, $ref, $alt), "\n"
        if $debug > 0;

    my @alleles = split(/\,/, $alt);

    if (!exists $seq_region_id{$chr}){
        warn "cant find seq_region_id for '$chr'\n";
        next;
    }

    ## Calculate the variation 'class' (id) from the alleles
    my $class_attrib_id =
        find_class_attrib($ref, \@alleles);

    ## Print the variation
    print { $file{variation} }
    join("\t",
         $variation_id,
         $source_id,
         $id, # name
         # validation_status
         # ancestral_allele
         # flipped
         $class_attrib_id,
         # somatic
         # minor_junk
         # clinical_significance
         # evidence_attribs
        ), "\n";

    ## and the variation feature
    print { $file{variation_feature} }
    join("\t",
         $variation_id,
         $seq_region_id{$chr},
         $pos,
         $pos + length($ref) - 1,
         1, # strand
         $variation_id,
         join('/', $ref, @alleles),
         $id, # variation_name
         1, # map_weight
         # flags
         $source_id,
         # validation_status
         # consequence_types
         # variation_set_id
         $class_attrib_id,
         # somatic
         # minor_junk
         # alignment quality
         # evidence_attribs
         # clinical_significance
        ), "\n";



    ## Hash all alleles into allele codes
    make_allele_codes( [$ref, @alleles] );



    ## Start processing individual genotypes

    my $individual_id;

    ## We need these counts to calculate frequencies later
    my %population_allele_count;
    my %population_genotype_count;

    for my $sample ( @sample ){
        $individual_id++;

        next if $sample eq '.';

        ## Debugging
        print "\t", $sample, "\n"
            if $debug > 0;

        ## Hash all genotypes into genotype codes
        my @genotype =
            make_genotype_codes( [$ref, @alleles], $sample );
        die if @genotype != 2; # ug

        ## WTF, quite frankly
        my $file = $file{individual_genotype_y};
        if($class_attrib_id == 2){ ## SNV
            $file = $file{individual_genotype_x};
        }

        print $file
            join("\t",
                 $variation_id,
                 # subsnp_id
                 $genotype[0],
                 $genotype[1],
                 $individual_id
                ), "\n";

        ## Collect the population allele counts for later
        $population_allele_count{$_}++
            for @genotype;

        ## Collect the population genotype counts for later
        $population_genotype_count{ join('/', @genotype) }++;
    }



    ## Move on to the final two tables

    my $allele_count_total;
    $allele_count_total += $_
        for values %population_allele_count;

    ## Note we loop through all possible alleles, not just those we
    ## see in individuals, e.g. the ref may never be seen.
    for ($ref, @alleles){
        $allele_id++;
        print { $file{population_allele} }
        join("\t",
             $allele_id,
             $variation_id,
             # subsnp_id
             $allele_code{$_},
             $population_id,
             ($population_allele_count{$_} || 0) / $allele_count_total,
             ($population_allele_count{$_} || 0),
             # frequency_submitter_handle
            ), "\n";
    }

    my $genotype_count_total;
    $genotype_count_total += $_
        for values %population_genotype_count;

    ## Note, we only loop through those genotypes we've seen!
    for (keys %population_genotype_count){
        $genotype_id++;
        print { $file{population_genotype} }
        join("\t",
             $genotype_id,
             $variation_id,
             # subsnp_id
             $genotype_code{$_},
             $population_genotype_count{$_} / $genotype_count_total,
             $population_id,
             $population_genotype_count{$_},
            ), "\n";
    }
}

print "got ", scalar keys %allele_code,   " alleles\n";
print "got ", scalar keys %genotype_code, " genotypes\n";



## Finally, dump them here...
print { $file{allele_code} }   "$allele_code{$_}\t$_\n"
    for keys %allele_code;

for my $genotype (keys %genotype_code){
    my $haplotype_id = 0;
    for my $allele (split(/\//, $genotype)){
        $haplotype_id++;
        print { $file{genotype_code} }
        join("\t",
             $genotype_code{ $genotype },
             $allele_code{ $allele },
             $haplotype_id,
             0, # phased
            ), "\n";
    }
}

warn "DONE\n";



sub find_class_attrib {
    my $ref = shift;
    my $alt_aref = shift;

    my $ref_len = length($ref);

    my $class_attrib =
        $ref_len == 1 ? 'SNV' : 'substitution';

    my %indel;
    for my $alt (@$alt_aref){
        $indel{length($alt) <=> $ref_len}++;
    }

    if(0){}

    elsif(defined $indel{+1}){
        if(defined $indel{-1}){
            $class_attrib = 'sequence_alteration';
        }
        else{
            $class_attrib = 'insertion';
        }
    }
    elsif(defined $indel{-1}){
        $class_attrib = 'deletion';
    }

    ## Return the class attrib as an id
    return $class_attrib_id{ $class_attrib };
}

sub make_allele_codes {
    my $allele_aref = shift;

    for my $allele ( @$allele_aref ){
        next if defined( $allele_code{$allele} );
        $allele_code{ $allele } = 1 + scalar keys %allele_code;
    }
}

sub make_genotype_codes {
    my $allele_aref = shift;
    my $sample = shift;

    ## This is just the way VCF is
    my @genotype =
        @$allele_aref[ split('/', (split(':', $sample))[0]) ];
    my $genotype = join('/', @genotype);
    $genotype_code{ $genotype } = 1 + scalar keys %genotype_code
        unless defined $genotype_code{ $genotype };

    # usefull side effect:
    return @genotype;
}
