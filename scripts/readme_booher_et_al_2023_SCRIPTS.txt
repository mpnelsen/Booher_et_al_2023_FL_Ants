#1_ants_for_doug7.r

#conducts an initial query of GenBank for gene sequences of interest from valid species listed in AntWiki. Note - this searches for ALL species in list, not just those in study area.

#subset_ncbi_for_doug.r
#subsets accessions to only include those from the study area.

#2_doug_join_fl_ga_add_alts_ncbi.r
#further editing including addition of alternate accessions/taxa.

#3_multi_seq_fetcher_17nov19.py
#script to query ncbi and retrieve accessions.

#4_multi_parser_revise_17nov19.py
#parses out locus of interest from accessions previously retrieved.

#5_name_changer_17nov19.r
#manually add in COI additional sequences.

#6_blasting_ants_tax_hybrid_family2_17nov19.py
#internal blast search for quality control

#trim.hybrid.family3_17nov19.r
#make summary files/clean up after blast

#7_make_out_files_17nov19.r
#make fasta file with outgroup sequences for each locus

#8_add.outs.to.master_17nov19.r
#add outgroup accession numbers to spreadsheet

#9_concatenate_outs_and_reals_17nov19
#adds outgroup sequences to ingroup sequences for each locus

#10_initial_align_fullalign_einsi_nogappy12_defaultmacse_17nov19.r
#perform initial alignments for each clade

#11_profile_aligning_fullalign_einsi_nogappy6_short_altorg_sparseqs2_defaultmacse_correctedlength_17nov19.r
#use profile aligning to align alignments of individual clades to one another

#ultimately removed Tetramorium immigrans and Tetramorium tsushimae, reducing alignment to 195 species.  In addition, Camponotus and Crematogaster were each constrained to be monophyletic in the RAxML analysis.

#congruify195tips.r
#use congruification to timescale new phylogeny to that of Nelsen et al. (2016, PNAS)

#add_doug_strumigenys.r
#cut out Strumigenys and replace with timscaled Strumigenys phylogeny including species of interest from Booher et al. (2021 - PLoS Bio)

