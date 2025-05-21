import os
import sys
import yaml
import hail as hl
import pyspark
import pandas as pd
from gnomad.sample_qc.ancestry import assign_population_pcs, pc_project
from shutil import rmtree

####################################################
# Step0 utils #
####################################################
def clear_temp_folder(tmp_dir: str) -> None:
    if not tmp_dir.startswith("file://"):
        return
    tmp_dir = tmp_dir.replace("file://", "")
    print(f"=== Cleaning up temporary folder {tmp_dir}")
    try:
        rmtree(tmp_dir)
    except FileNotFoundError:
        print("TMP folder was cleaned up by Spark. Proceeding to the next step.")
    os.makedirs(tmp_dir, exist_ok=True)

####################################################
# Step1 load VCF #
####################################################
def load_vcfs_to_mt(indir, mt_out_file, tmp_dir):
    '''
    load VCFs and save as hail mt
    '''
    objects = hl.utils.hadoop_ls(indir)
    vcfs = [vcf["path"] for vcf in objects if (vcf["path"].startswith("file") and vcf["path"].endswith("vcf.gz"))]
    print("Loading VCFs")
    #create and save MT
    mt = hl.import_vcf(vcfs, array_elements_required=False, force_bgz=True, header_file=None)
    #separate mt
    print("Saving as hail mt")
    mt.write(mt_out_file, overwrite=True)

####################################################
# Step2 hard filters #
####################################################

def apply_hard_filters(mt_file: str, filtered_mt_file: str) -> str:
    '''
    Applies hard filters and annotates samples in the filtered set with call rate
    :param MatrixTable mt: MT containing samples to be ascertained for sex
    :param str mtdir: directory output matrix tables are written to
    :return: MatrixTable with hard filtering annotation
    :rtype: MatrixTable
    '''
    mt = hl.read_matrix_table(mt_file)
    print("Applying hard filters")
    mt = mt.filter_rows((hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1]) &
        (hl.agg.mean(mt.GT.n_alt_alleles()) / 2 > 0.001) &
        (hl.agg.fraction(hl.is_defined(mt.GT)) > 0.99))
    mt = mt.annotate_cols(callrate=hl.agg.fraction(hl.is_defined(mt.GT)))
    print("Writing to " + filtered_mt_file)
    mt.write(filtered_mt_file, overwrite=True)

####################################################
# Step3 sex imputation #
####################################################

def impute_sex(filename: str, sex_mt_file: str, outdir: str, male_threshold: float = 0.8, female_threshold: float = 0.5) -> str:
    '''
    Imputes sex, exports data, and annotates mt with this data
    :param MatrixTable mt: MT containing samples to be ascertained for sex
    :param str mtdir: directory output matrix tables are written to
    :param str outdir: directory result files are written to
    :return: MatrixTable with imputed sex annotations stashed in column annotation 'sex_check'
    :rtype: MatrixTable
    '''
    print("Imputing sex with male_threshold = " + str(male_threshold) + " and female threshold = " + str(female_threshold))
    mt = hl.read_matrix_table(filename)

    #filter to X and select unphased diploid genotypes - no need to filter to X as impute_sex will do this
    #mt1 = hl.filter_intervals(mt, [hl.parse_locus_interval('chrX')])
    mt1 = hl.split_multi_hts(mt)
    mtx_unphased = mt1.select_entries(GT=hl.unphased_diploid_gt_index_call(mt1.GT.n_alt_alleles()))
    #imput sex on the unphased diploid GTs
    sex_ht = hl.impute_sex(mtx_unphased.GT, aaf_threshold=0.05, female_threshold=female_threshold, male_threshold=male_threshold)
    #export
    outfile = "file://" + os.path.join(outdir, 'sex_annotated.sex_check.txt')
    print("Writing to " + outfile)
    sex_ht.export(outfile)
    #annotate input (all chroms) mt with imputed sex and write to file
    sex_colnames = ['f_stat', 'is_female']
    sex_ht = sex_ht.select(*sex_colnames)
    mt = mt.annotate_cols(**sex_ht[mt.col_key])
    print("Writing to " + sex_mt_file)
    mt.write(sex_mt_file, overwrite=True)

####################################################
# Step4 simple QC #
####################################################
def simple_qc(mtfile: str, bed_file: str, outdir: str):
    '''
    simple QC
    '''
    print("Running QC")
    bed_intervals = hl.import_bed(bed_file, reference_genome='GRCh38')
    mt=hl.read_matrix_table(mtfile)
    exome_mt = mt
    exome_mt = exome_mt.annotate_rows(ug_hcr=hl.is_defined(bed_intervals[exome_mt.locus]))
    # Filter for exome
    exome_mt = exome_mt.filter_rows(hl.is_defined(bed_intervals[exome_mt.locus]))
    # Filter for quality
    exome_mt = exome_mt.filter_entries((exome_mt.DP >= 20) & (exome_mt.GQ >= 20) & ((exome_mt.AD[1] / hl.sum(exome_mt.AD)) >= 0.25))
    # Run sample_qc exomes
    exome_mt = hl.sample_qc(exome_mt, name='exome_sample_qc')
    # Convert to pandas and get snps and indels
    exome_df = exome_mt.cols().to_pandas()
    # Get only some of the columns
    exome_df = exome_df[["s", "exome_sample_qc.n_snp", "exome_sample_qc.n_insertion", "exome_sample_qc.n_deletion", "exome_sample_qc.n_transversion", "exome_sample_qc.n_transition", "exome_sample_qc.r_ti_tv", "exome_sample_qc.r_het_hom_var", "exome_sample_qc.r_insertion_deletion", "exome_sample_qc.n_het"]]
    # Export
    df_file=os.path.join(outdir, 'exome_qc_filtered.csv')
    print("Writing to " + df_file)
    exome_df.to_csv(df_file, index=True)
    
    # Filter for non exomes
    non_exome_mt = mt
    non_exome_mt = non_exome_mt.annotate_rows(ug_hcr=hl.is_defined(bed_intervals[non_exome_mt.locus]))
    # Filter for rows where the locus is NOT in the exome regions
    non_exome_mt = non_exome_mt.filter_rows(~hl.is_defined(bed_intervals[non_exome_mt.locus]))
    # Filter for quality
    non_exome_mt = non_exome_mt.filter_entries((non_exome_mt.DP >= 4) & (non_exome_mt.GQ >= 15))
    # Filter for vaf
    non_exome_mt = non_exome_mt.filter_entries((non_exome_mt.AD[1] / hl.sum(non_exome_mt.AD)) >= 0.25)
    # Run sample_qc non_exomes
    non_exome_mt = hl.sample_qc(non_exome_mt, name='non_exome_sample_qc')
    # Convert to pandas
    non_exome_df = non_exome_mt.cols().to_pandas()
    # Subset for the columns you want
    non_exome_df = non_exome_df[["s", "non_exome_sample_qc.n_snp", "non_exome_sample_qc.n_insertion", "non_exome_sample_qc.n_deletion", "non_exome_sample_qc.n_transversion", "non_exome_sample_qc.n_transition", "non_exome_sample_qc.r_ti_tv", "non_exome_sample_qc.r_het_hom_var", "non_exome_sample_qc.r_insertion_deletion", "non_exome_sample_qc.n_het"]]
    # Export
    df_file=os.path.join(outdir, 'non_exome_qc_filtered.csv')
    print("Writing to " + df_file)
    non_exome_df.to_csv(df_file, index=True)

####################################################
# Step5 1kg processing #
####################################################

def filter_mt_autosome_biallelic_snvs(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Keeps only bi-allelic autosome SNVs
    """
    mt = mt.filter_rows(mt.locus.in_autosome())
    # split multiallelic variants and remove them
    mt = hl.split_multi_hts(mt)  # this shouldn't do anything as only biallelic sites are used
    mt = mt.filter_rows(mt.was_split is True, keep=False)
    # keep only SNVs
    mt = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1]))
    return mt

def filter_vars_for_quality(
    mt: hl.MatrixTable, af_threshold: float, call_rate_threshold: float, hwe_threshold: float
) -> hl.MatrixTable:
    """
    Keeps "good" variants passing the thresholds for QC metrics
    """
    mt_vqc = hl.variant_qc(mt, name="variant_QC_Hail")
    mt_vqc_filtered = mt_vqc.filter_rows(
        (mt_vqc.variant_QC_Hail.call_rate >= call_rate_threshold)
        & (mt_vqc.variant_QC_Hail.AF[1] >= af_threshold)
        & (mt_vqc.variant_QC_Hail.p_value_hwe >= hwe_threshold)
    )
    return mt_vqc_filtered

def remove_ld_regions(mt: hl.MatrixTable, long_range_ld_file: str) -> hl.MatrixTable:
    # remove long ld regions
    long_range_ld_to_exclude = hl.import_bed(long_range_ld_file, reference_genome="GRCh38")
    mt = mt.filter_rows(hl.is_defined(long_range_ld_to_exclude[mt.locus]), keep=False)
    return mt

def remove_palindromes(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Removes palindromic SNVs from the matrix table
    A palindromic SNP (also known as an "ambiguous" SNP) is a SNP
    in which the possible alleles for that SNP are the same alleles
    that would pair with each other in the double helix structure.
    e.g., C/G on forward is G/C on the reverse
    https://mr-dictionary.mrcieu.ac.uk/term/palindrome/

    To get rid of these SNVs we remove all paired combinations: A<->T, G<->C
    """
    mt_non_pal = mt.filter_rows((mt.alleles[0] == "G") & (mt.alleles[1] == "C"), keep=False)
    mt_non_pal = mt_non_pal.filter_rows((mt_non_pal.alleles[0] == "C") & (mt_non_pal.alleles[1] == "G"), keep=False)
    mt_non_pal = mt_non_pal.filter_rows((mt_non_pal.alleles[0] == "A") & (mt_non_pal.alleles[1] == "T"), keep=False)
    mt_non_pal = mt_non_pal.filter_rows((mt_non_pal.alleles[0] == "T") & (mt_non_pal.alleles[1] == "A"), keep=False)
    return mt_non_pal

def filter_matrix_for_ldprune(mt: hl.MatrixTable, long_range_ld_file: str, call_rate_threshold: float = 0.99, af_threshold: float = 0.05, hwe_threshold: float = 1e-5,) -> hl.MatrixTable:
    """
    This function filters samples and variants to perform LD pruning,
    pruning of related samples

    The function is taken from the step 2/3-annotate_and_filter
    In future step 2/3 should be refactored to use this function
    """
    mt = filter_mt_autosome_biallelic_snvs(mt)
    mt_vqc_filtered = filter_vars_for_quality(mt, af_threshold, call_rate_threshold, hwe_threshold)
    mt_vqc_filtered = remove_ld_regions(mt_vqc_filtered, long_range_ld_file)
    mt_non_pal = remove_palindromes(mt_vqc_filtered)
    return mt_non_pal

def create_1kg_mt(vcf_indir: str, kg_pop_file: str, kg_unprocessed_mt_file: str):
    """
    Create matrixtable of 1kg data
    :param str vcf_indir: the directory with 1KG VCF files
    :param str kg_pop_file: Assignes superpopulations
    """
    print(f"Loading VCFs from {vcf_indir}")
    objects = hl.utils.hadoop_ls(vcf_indir)
    vcfs = [vcf["path"] for vcf in objects if (vcf["path"].startswith("file") and vcf["path"].endswith("vcf.gz"))]
    # create MT
    kg_unprocessed_mt = hl.import_vcf(vcfs, array_elements_required=False, force_bgz=True)
    # Annotating known populations
    #kg_pop_file = path_spark(kg_pop_file)
    cohorts_pop = hl.import_table(kg_pop_file, delimiter="\t").key_by("Sample name")
    kg_unprocessed_mt = kg_unprocessed_mt.annotate_cols(
        known_pop=cohorts_pop[kg_unprocessed_mt.s]["Superpopulation code"]
    )
    print("Writing to " + kg_unprocessed_mt_file)
    kg_unprocessed_mt.write(kg_unprocessed_mt_file, overwrite=True)

def kg_filter_and_ldprune(kg_unprocessed_mt_file: hl.MatrixTable, long_range_ld_file: str, pruned_kg_file: str):
    """
    Filter and prune the 1kg data
    :param kg_unprocessed_mt: The KG MT to filter and prune
    :param long_range_ld_file: The long range LD file
    :param call_rate_threshold: The call rate threshold
    :param af_threshold: The allele frequency threshold
    :param hwe_threshold: The HWE threshold
    :param r2_threshold: Squared correlation threshold: https://hail.is/docs/0.2/methods/genetics.html#hail.methods.ld_prune
    :return: The filtered and pruned 1kg MT
    """
    #long_range_ld_file = path_spark(long_range_ld_file)
    kg_unprocessed_mt=hl.read_matrix_table(kg_unprocessed_mt_file)
    call_rate_threshold=0.99
    af_threshold=0.05
    hwe_threshold=0.00001#1e-5
    r2_threshold=0.2
    print("Filtering and pruning of 1kg matrix")
    # Filtering for good variations to make LD prune
    kg_mt_filtered = filter_matrix_for_ldprune(
        kg_unprocessed_mt, long_range_ld_file, call_rate_threshold, af_threshold, hwe_threshold
    )
    # LD pruning - removing variation regions that are related to each other
    pruned_kg_ht = hl.ld_prune(kg_mt_filtered.GT, r2=r2_threshold)
    pruned_kg_mt = kg_mt_filtered.filter_rows(hl.is_defined(pruned_kg_ht[kg_mt_filtered.row_key]))
    pruned_kg_mt = pruned_kg_mt.select_entries(GT=hl.unphased_diploid_gt_index_call(pruned_kg_mt.GT.n_alt_alleles()))
    print("Writing to " + pruned_kg_file)
    pruned_kg_mt.write(pruned_kg_file, overwrite=True)

def run_pc_relate( pruned_kg_file: str, relatedness_ht_file: str, scores_file: str, samples_to_remove_file : str):
    """
    Runs PC relate on pruned MT
    :param str pruned_mt: matrixtable to prune
    :param str relatedness_ht_file: relatedness ht file
    :param str scores_file: file to wtire scores ht
    :param int pca_components:  the number of principal components
    :param dict hl_pc_related_kwargs: kwargs to pass to HL PC relate
    """
    #relatedness_ht_file = path_spark(relatedness_ht_file)
    #scores_file = path_spark(scores_file)
    pruned_mt=hl.read_matrix_table(pruned_kg_file)
    pca_components=10
    kin_threshold=0.125

    print("Running PC relate")
    eig, scores, _ = hl.hwe_normalized_pca(pruned_mt.GT, k=pca_components, compute_loadings=False)
    scores.write(scores_file, overwrite=True)

    print("Calculating relatedness (this step usually takes a while)")
    relatedness_ht = hl.pc_relate(pruned_mt.GT, scores_expr=scores[pruned_mt.col_key].scores, min_individual_maf=0.05, block_size=4096, min_kinship=0.05, statistics="kin2")
    relatedness_ht.write(relatedness_ht_file, overwrite=True)
    # prune individuals to be left with unrelated - creates a table containing one column - samples to remove
    pairs = relatedness_ht.filter(relatedness_ht["kin"] > kin_threshold)
    related_samples_to_remove = hl.maximal_independent_set(pairs.i, pairs.j, keep=False)
    print("Writing to " + samples_to_remove_file)
    related_samples_to_remove.write(samples_to_remove_file, overwrite=True)

def kg_remove_related_samples(kg_unprocessed_mt_file: str, samples_to_remove_file: str, kg_mt_file:str):
    kg_mt=hl.read_matrix_table(kg_unprocessed_mt_file)
    related_samples_to_remove=hl.read_table(samples_to_remove_file)
    variants, samples = kg_mt.count()
    print(f"Loaded form initial table: {samples} samples, {variants} variants.")
    print("Removing related samples")
    kg_mt_remove_related = kg_mt.filter_cols(hl.is_defined(related_samples_to_remove[kg_mt.col_key]), keep=False)
    variants, samples = kg_mt_remove_related.count()
    print(f"Remains after removing related samples: {samples} samples, {variants} variants.")
    print("Writing to " + kg_mt_file)
    kg_mt_remove_related.write(kg_mt_file, overwrite=True)

def get_kg_mt(kg_vcf_dir, long_range_ld_file, pops_file, mtdir, kg_mt_file):
    print("Creating 1kg matrix")
    kg_unprocessed_mt_file=os.path.join(mtdir, "kg_unprocessed.mt")
    create_1kg_mt(kg_vcf_dir, pops_file, kg_unprocessed_mt_file)
    
    pruned_kg_file=os.path.join(mtdir, "kg_pruned.mt")
    kg_filter_and_ldprune(kg_unprocessed_mt_file, long_range_ld_file, pruned_kg_file)
    
    relatedness_ht_file=os.path.join(mtdir, "kg_relatedness.ht")
    scores_file=os.path.join(mtdir, "kg_pruned.pca_scores")
    samples_to_remove_file=os.path.join(mtdir, "sample_to_remove.ht")
    run_pc_relate(pruned_kg_file, relatedness_ht_file, scores_file, samples_to_remove_file)

    kg_remove_related_samples(kg_unprocessed_mt_file, samples_to_remove_file, kg_mt_file)

####################################################
# Step6 PCA #
####################################################

def read_and_filter (mt_file):
    mt = hl.read_matrix_table(mt_file)
    mt = mt.filter_rows(mt.locus.in_autosome())
    mt = hl.split_multi_hts(mt)#this shouldn't do anything as only biallelic sites are used
    mt = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1]))
    mt = mt.filter_rows(mt.was_split==True, keep=False)
    return mt

def filter_matrix (mt: hl.MatrixTable, long_range_ld_file: str):
    #use only autosomes
    mt=mt.filter_rows(mt.locus.in_autosome())
    #split multiallelic variants and remove them
    mt = hl.split_multi_hts(mt)#this shouldn't do anything as only biallelic sites are used
    mt = mt.filter_rows(mt.was_split==True, keep=False)
    #keep only SNPs
    mt = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1]))
    #keep good variants using hail variant_qc and thre filters
    mt_vqc = hl.variant_qc(mt, name='variant_QC_Hail')
    mt_vqc_filtered = mt_vqc.filter_rows(
        (mt_vqc.variant_QC_Hail.call_rate >= 0.99) &
        (mt_vqc.variant_QC_Hail.AF[1] >= 0.05) &
        (mt_vqc.variant_QC_Hail.p_value_hwe >= 1e-5)
    )
    #remove long ld regions
    long_range_ld_to_exclude = hl.import_bed(long_range_ld_file, reference_genome='GRCh38')
    mt_vqc_filtered = mt_vqc_filtered.filter_rows(hl.is_defined(long_range_ld_to_exclude[mt_vqc_filtered.locus]), keep=False)
    #remove palindromes
    mt_non_pal = mt_vqc_filtered.filter_rows((mt_vqc_filtered.alleles[0] == "G") & (mt_vqc_filtered.alleles[1] == "C"), keep=False)
    mt_non_pal = mt_non_pal.filter_rows((mt_non_pal.alleles[0] == "C") & (mt_non_pal.alleles[1] == "G"), keep=False)
    mt_non_pal = mt_non_pal.filter_rows((mt_non_pal.alleles[0] == "A") & (mt_non_pal.alleles[1] == "T"), keep=False)
    mt_non_pal = mt_non_pal.filter_rows((mt_non_pal.alleles[0] == "T") & (mt_non_pal.alleles[1] == "A"), keep=False)
    return mt_non_pal

def merge_and_prun (mt: hl.MatrixTable, kg_mt: hl.MatrixTable, long_range_ld_file: str, pops_file: str, mtdir: str, pruned_mt_file: str):
    print("Merging and pruning matrices before PCA")
    #filter both matrices
    mt_filtered=filter_matrix(mt, long_range_ld_file)
    kg_filtered=filter_matrix(kg_mt, long_range_ld_file)
    #removing and adding needed entries to replicate filtered_mt_file structure 
    mt_filtered=mt_filtered.drop('AD', 'DP', 'GQ', 'MIN_DP', 'PGT', 'PID', 'PL', 'PS', 'SB', 'callrate', 'f_stat', 'is_female', 'RGQ')
    kg_filtered = kg_filtered.select_entries(kg_filtered.GT)
    #merging matrices
    mt_filtered = mt_filtered.annotate_cols(known_pop = hl.null(hl.tstr))
    mt_merged = mt_filtered.union_cols(kg_filtered)
    #saving the merged matrix
    merged_mt_file=os.path.join(mtdir, "merged_with_1kg_filtered.mt")
    print("Writing to " + merged_mt_file)
    mt_merged=mt_merged.checkpoint(merged_mt_file, overwrite=True)
    #prunning
    pruned_ht = hl.ld_prune(mt_merged.GT, r2=0.2)
    pruned_mt = mt_merged.filter_rows(hl.is_defined(pruned_ht[mt_merged.row_key]))
    pruned_mt = pruned_mt.select_entries(GT=hl.unphased_diploid_gt_index_call(pruned_mt.GT.n_alt_alleles()))
    #annotating known poppulations
    cohorts_pop = hl.import_table(pops_file, delimiter="\t").key_by('Sample name')
    pruned_mt = pruned_mt.annotate_cols(known_pop=cohorts_pop[pruned_mt.s]['Superpopulation code'])
    #saving matrix
    print("Writing to " + pruned_mt_file)
    pruned_mt.write(pruned_mt_file, overwrite=True)

def run_pca(filtered_mt_file: str, pca_scores_file: str, mtdir: str):
    print("Running PCA")
    #read filtered matrix
    mt = hl.read_matrix_table(filtered_mt_file)

    #divide matrix to make a projection
    mt_kg = mt.filter_cols(hl.is_defined(mt.known_pop))
    mt_study = mt.filter_cols(hl.is_missing(mt.known_pop))

    #PCA for 1000 Genomes
    pca_evals, pca_scores, pca_loadings = hl.hwe_normalized_pca(mt_kg.GT, k=10, compute_loadings=True)
    pca_scores = pca_scores.annotate(known_pop=mt_kg.cols()[pca_scores.s].known_pop)
    pca_af_ht = mt_kg.annotate_rows(pca_af=hl.agg.mean(mt_kg.GT.n_alt_alleles()) / 2).rows()
    pca_loadings = pca_loadings.annotate(pca_af=pca_af_ht[pca_loadings.key].pca_af)
    #saving files
    pca_1kg_scores_file = os.path.join(mtdir, "pca_scores_1Kg.ht")
    pca_1kg_loadings_file = os.path.join(mtdir, "pca_loadings_1Kg.ht")
    print("Writing to " + pca_1kg_scores_file)
    pca_scores.write(pca_1kg_scores_file, overwrite=True)
    print("Writing to " + pca_1kg_loadings_file)
    pca_loadings.write(pca_1kg_loadings_file, overwrite=True)
    pca_scores = pca_scores.drop(pca_scores.known_pop)
    #projection of samples on precomputed PCs and combining of two PCA_scores tables
    projection_PCA_scores=pc_project(mt_study, pca_loadings, loading_location='loadings', af_location='pca_af')
    union_PCA_scores=pca_scores.union(projection_PCA_scores)
    union_PCA_scores = union_PCA_scores.annotate(known_pop=mt.cols()[union_PCA_scores.s].known_pop)
    #saving union_PCA_scores
    print("Writing to " + pca_scores_file)
    union_PCA_scores.write(pca_scores_file, overwrite=True)

def run_prun_and_pca(filtered_mt_file: str, kg_mt_file: str, long_range_ld_file: str, pops_file: str, bed_file: str, mtdir: str, pca_scores_file: str):
    '''
    Run PCA before population prediction
    :param str filtered_mt_file: merged birth cohort wes and 1kg MT file annotated with pops and filtered
    :param str pca_sores_file: PCA scores HT file
    :param str pca_loadings_file: PCA scores HT file
    :param str pca_evals_file: PCA scores HT file
    '''
    mt = read_and_filter(filtered_mt_file)
    kg_mt = read_and_filter(kg_mt_file)
    kg_mt = kg_mt.select_entries(kg_mt.GT)
    bed_intervals = hl.import_bed(bed_file, reference_genome='GRCh38')
    
    #kg_mt = kg_mt.annotate_rows(ug_hcr=hl.is_defined(bed_intervals[kg_mt.locus]))
    # Filter for exome
    kg_mt = kg_mt.filter_rows(hl.is_defined(bed_intervals[kg_mt.locus]))

    pruned_mt_file=os.path.join(mtdir, "pruned_with_1kg_filtered.mt")
    merge_and_prun(mt, kg_mt, long_range_ld_file, pops_file, mtdir, pruned_mt_file)

    run_pca(pruned_mt_file, pca_scores_file, mtdir)

####################################################
# Step7 Population assignment #
####################################################

def predict_pops(pca_scores_file: str, mtdir: str, pop_ht_tsv: str):
    '''
    Predict populations from PCA scores
    :param str pca_sores_file: PCA scores HT file
    :param str pop_ht_file: predicted population HT file
    :param str pop_ht_tsv: population tsv file
    '''
    print("Population assignment")
    pca_scores = hl.read_table(pca_scores_file)
    known_col = "known_pop"
    pop_ht, pop_clf = assign_population_pcs(pca_scores, pca_scores.scores, known_col=known_col, n_estimators=100, prop_train=0.8, min_prob=0.5)
    pop_ht_file=os.path.join(mtdir, "pop_assignments.ht")
    print("Writing to " + pop_ht_file)
    pop_ht.write(pop_ht_file, overwrite=True)
    #convert to pandas and put in only pops files, add excluded sample back
    pop_ht_df = pop_ht.to_pandas()
    print("Writing to " + pop_ht_tsv)
    pop_ht_df.to_csv(pop_ht_tsv, sep = '\t')

def main():
    #set up input variables
    if len(sys.argv) < 2:
        print("Usage: PYSPARK_DRIVER_PYTHON=/home/ubuntu/venv/bin/python spark-submit run_wes_qc.py <config_file>")
        sys.exit(1)

    input_yaml = sys.argv[1]

    with open(input_yaml, 'r') as y:
        inputs = yaml.load(y, Loader=yaml.FullLoader)

    mtdir = "file://" + inputs['matrixtables_lustre_dir']
    outdir = inputs['results_dir']
    #control_list = inputs['control_samples']
    #initialise hail
    tmp_dir = "file://" + inputs["tmp_dir"]
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    #load VCFs
    import_vcf_dir = "file://" + inputs['gatk_import_lustre_dir']
    mt_unprocessed_file =os.path.join(mtdir, "gatk_unprocessed.mt")
    if not inputs['skip_vcf_load']:
        load_vcfs_to_mt(import_vcf_dir, mt_unprocessed_file, tmp_dir)
    
    #apply hard fitlers
    filtered_mt_file =os.path.join(mtdir, "mt_hard_filters_annotated.mt")
    if not inputs['skip_hard_filters']:
        apply_hard_filters(mt_unprocessed_file, filtered_mt_file)

    #impute sex
    sex_mt_file =os.path.join(mtdir, "mt_sex_annotated.mt")
    if not inputs['skip_impute_sex']:
        impute_sex(filtered_mt_file, sex_mt_file, outdir, male_threshold=0.79, female_threshold=0.55)
    
    #simple qc
    bed_file = "file://" + inputs['hcr_bed']
    if not inputs['skip_simple_qc']:
        simple_qc(mt_unprocessed_file, bed_file, outdir)

    #make 1000 genome matrix
    long_range_ld_file="file://" + inputs['long_range_ld']    
    pops_file="file://" + inputs['kg_pop']
    kg_mt_file=os.path.join(mtdir, "kg_filtered.mt")
    kg_vcf_dir="file://" + inputs['kg_vcf']
    if not inputs['skip_1kg_processing']:
        get_kg_mt(kg_vcf_dir, long_range_ld_file, pops_file, mtdir, kg_mt_file)

    #PCA
    pca_scores_file=os.path.join(mtdir, "pca_scores_after_pruning.ht")
    if not inputs['skip_pca']:
        run_prun_and_pca(sex_mt_file, kg_mt_file, long_range_ld_file, pops_file, bed_file, mtdir, pca_scores_file)

    #predict populations
    pop_ht_tsv = os.path.join(outdir, "pop_assignemtnts.tsv")
    if not inputs['skip_assignments']:
        predict_pops(pca_scores_file, mtdir, pop_ht_tsv)

    clear_temp_folder(tmp_dir)

if __name__ == '__main__':
    main() 
