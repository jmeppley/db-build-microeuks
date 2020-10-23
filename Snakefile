# Snakemake file to import fasta and tax files to generate and classify database for use in QIIME2
configfile: "config.yaml"

import re

#----SET VARIABLES----#

SCRATCHDIR = config["scratch"]
NAME = config["db_name"]
REGION = config["db_region"]
VERSION = config["db_version"]
PRIMER_F = config["primerF"]
PRIMER_R = config["primerR"]

#----INPUT FILES----#

# path or url for reference fasta files
INPUTDB = config["ref_db"]
# path or url to tax file, or extract from fasta header if not supplied
INPUTTAX = config.get("ref_tax", \
                      f"{SCRATCHDIR}/{NAME}/{NAME}_{VERSION}.tax")


#----DEFAULT CONFIG VALUES----#

config.setdefault('p_trunc_len', 0)
config.setdefault('max_threads', 20)

#----DEFINE RULES----#

rule all:
  """ generates the classify.qza file and runs a test using the input sequences """
  input:
    tax	= f"{SCRATCHDIR}/{NAME}/{REGION}-{NAME}_{VERSION}-outputtax-TEST.qza",
    taxviz = f"{SCRATCHDIR}/{NAME}/{REGION}-{NAME}_{VERSION}-outputtax-TEST.qzv",

rule no_test:
  """ Just generates the classif.qza database """
  input:
    classifier = f"{SCRATCHDIR}/{NAME}/{REGION}-{NAME}_{VERSION}-classifier.qza",

rule tax_from_header:
  """ extract tax file (an ID -> lineage map) from fasta header (Sivla style) """
  input:
    INPUTDB
  output:
    INPUTTAX
  shell: "grep "^>" {input} \
          | perl -pe 's/^>(\\S+)\\s+/\\1\\t/' \
          > {output}"

rule import_db:
  input:
    inputdb = INPUTDB
  output:
    dbartifact = f"{SCRATCHDIR}/{NAME}/{NAME}_{VERSION}.qza"
  log:
    f"{SCRATCHDIR}/{NAME}/logs/db-import.log"
  shell:
   """
   qiime tools import \
      --type 'FeatureData[Sequence]' \
      --input-path {input.inputdb} \
      --output-path {output.dbartifact}
   """

rule import_tax:
  input:
    tax = INPUTTAX
  output:
    dbtaxartifact = f"{SCRATCHDIR}/{NAME}/{NAME}_{VERSION}_tax.qza"
  log:
    f"{SCRATCHDIR}/{NAME}/logs/db-import-tax.log"
  shell:
    """
    qiime tools import \
      --type 'FeatureData[Taxonomy]' \
      --input-format HeaderlessTSVTaxonomyFormat \
      --input-path {input.tax} \
      --output-path {output.dbtaxartifact}
    """

rule subset_region:
  input:
    dbartifact = f"{SCRATCHDIR}/{NAME}/{NAME}_{VERSION}.qza"
  output:
    select = f"{SCRATCHDIR}/{NAME}/{REGION}-{NAME}_{VERSION}.qza"
  log:
    f"{SCRATCHDIR}/{NAME}/logs/db-subset-region.log"
  params:
    ptrunc = "--p-trunc-len {}".format(config['p_trunc_len']) \
             if int(config['p_trunc_len']) > 0 \
             else ""
  shell:
    """
    qiime feature-classifier extract-reads \
      --i-sequences {input.dbartifact} \
      --p-f-primer {config[primerF]} \
      --p-r-primer {config[primerR]} \
      {params.ptrunc} \
      --o-reads {output.select}
    """

rule classify:
  input:
    dbtaxartifact = f"{SCRATCHDIR}/{NAME}/{NAME}_{VERSION}_tax.qza",
    select = f"{SCRATCHDIR}/{NAME}/{REGION}-{NAME}_{VERSION}.qza"
  output:
    classifier = f"{SCRATCHDIR}/{NAME}/{REGION}-{NAME}_{VERSION}-classifier.qza"
  log:
    f"{SCRATCHDIR}/{NAME}/logs/db-subset-region-classifier.log"
  shell:
    """
    qiime feature-classifier fit-classifier-naive-bayes \
      --i-reference-reads {input.select} \
      --i-reference-taxonomy {input.dbtaxartifact} \
      --o-classifier {output.classifier}
    """

rule qc_check:
  input:
    classifier = f"{SCRATCHDIR}/{NAME}/{REGION}-{NAME}_{VERSION}-classifier.qza",
    dbartifact = f"{SCRATCHDIR}/{NAME}/{NAME}_{VERSION}.qza"
  output:
    tax = f"{SCRATCHDIR}/{NAME}/{REGION}-{NAME}_{VERSION}-outputtax-TEST.qza"
  log:
    f"{SCRATCHDIR}/{NAME}/logs/db-subset-region-testclassifier.log"
  threads: config['max_threads']
  shell:
    """
    qiime feature-classifier classify-sklearn \
      --i-classifier {input.classifier} \
      --i-reads {input.dbartifact} \
      --p-n-jobs {threads} \
      --o-classification {output.tax}
    """

rule gen_taxviz:
  input:
    tax = f"{SCRATCHDIR}/{NAME}/{REGION}-{NAME}_{VERSION}-outputtax-TEST.qza"
  output:
    taxviz = f"{SCRATCHDIR}/{NAME}/{REGION}-{NAME}_{VERSION}-outputtax-TEST.qzv"
  log:
    f"{SCRATCHDIR}/{NAME}/logs/db-subset-region-testclassifierViz.log"
  shell:
    """
    qiime metadata tabulate \
      --m-input-file {input.tax} \
      --o-visualization {output.taxviz}
    """
