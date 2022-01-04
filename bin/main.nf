#!/usr/bin/env nextflow

def helpMessage() {
  log.info"""
  -----------------------------------------------------------------------
                          BRAF_KRAS TEST DATA PIPELINE
  -----------------------------------------------------------------------
  Usage:

  nextflow run main.nf -profile singularity \
                       --refDir <path/to/refs> \
                       --bamCsv <path/to/bams> \
                       --bedInp <path/to/bed> \
                       --runID text_tag \

    -profile        [str]       Configuration profile
                                (required: singularity)

    --runID         [str]       Name for run, used to tag outputs

    --bedInp        [file]      Path to BED file of genes to get result VCF from

    --refDir        [file]      Path of dir in which reference data are held;
                                this should be created by download-references.nf
                                and contain dir "GRCh38"

    --bamCsv        [file]      CSV format, headers as sampleCsv except read1,
                                read2 are swapped for bam which is sent to
                                duplicate marking

    --multiqcConfig [str]       Config file for multiqc
                                (default: bin/somatic_n-of-1.multiQC_config.yaml)

    """.stripIndent()
}

if (params.help) exit 0, helpMessage()

//Test Mandatory Arguments

if(params.bamCsv == null){
  exit 1, "Please include --bamCsv, see --help for format"
}

if(!Channel.from(params.runID, checkIfExists: true)){
    exit 1, "Please include --runID <your_runID>"
}

if(!Channel.from(params.refDir, checkIfExists: true)){
  exit 1, "Please include --refDir <path> see github.com/brucemoran/somatic_n-of-1/ for how to run download-references.nf"
}

//Global Variables based on input
params.outDir = "output"
params.assembly = "GRCh38"
params.seqLevel = "WGS"
params.seqlevel = "wgs"
params.germline = false

//Java task memory allocation via task.memory
javaTaskmem = { it.replace(" GB", "g") }

//Reference data as value channels and reusable therefore
reference = [
    grchvers: false,
    fa: false,
    fai: false,
    dict: false,
    bwa: false,
    hc_dbs: false,
    dbsnp: false,
    gridss: false,
    pcgrbase: false,
    intlist: false,
    seqlevel: false,
    bbres: false
]

reference.grchvers  = Channel.fromPath("${params.refDir}/${params.assembly}/pcgr/data/*", type: 'dir').getVal()
reference.fa = Channel.value(file(params.genomes[params.assembly].fa))
reference.fai = Channel.value(file(params.genomes[params.assembly].fai))
reference.dict = Channel.value(file(params.genomes[params.assembly].dict))
reference.bwa = Channel.value(file(params.genomes[params.assembly].bwa))
reference.hc_dbs = Channel.value(file(params.genomes[params.assembly].hc_dbs))
reference.dbsnp = Channel.value(file(params.genomes[params.assembly].dbsnp))
reference.gridss = Channel.value(file(params.genomes[params.assembly].gridss))
reference.pcgrbase = Channel.value(file(params.genomes[params.assembly].pcgr))
reference.refflat = Channel.value(file(params.genomes[params.assembly].refflat))

//if seqlevel is exome, there is a dir per exome already parsed according to exomeTag
reference.seqlevel = Channel.value(file(params.genomes[params.assembly].wgs))

//setting of refBed
refbed = Channel.fromPath("${params.bedInp}")

process intlister {

  label 'low_mem'

  input:
  file(bed) from refbed
  file(dict) from reference.dict

  output:
  file('bed.interval_list') into ( intlist_call1, intlist_call2, intlist_call3 )

  script:
  """
  ##always make interval list so we are in line with fasta
  picard BedToIntervalList I=${bed} O=bed.interval_list1 SD=${dict}

  ##BedToIntervalList (reason unknown) makes 1bp interval to 0bp interval, replace with original
  perl -ane 'if(\$F[0]=~m/^@/){print \$_;next;} if(\$F[1] == \$F[2]){\$f=\$F[1]; \$f--; \$F[1]=\$f; print join("\\t", @F[0..\$#F]) . "\\n";} else{print \$_;}' bed.interval_list1 > bed.interval_list
  """
}

/*
================================================================================
                          -0. PREPROCESS INPUT SAMPLE FILE
================================================================================
*/
/* 0.00: Input using sample.csv, bam.csv, sample_cat
*/

Channel.fromPath("${params.bamCsv}")
       .splitCsv( header: true )
       .map { row -> [file(row.tumourbam), file(row.germlinebam)] }
       .set { which_bam }

process bam_input {

  label 'low_mem'

  input:
  tuple file(tumourbam), file(germlinebam) from which_bam

  output:
  tuple val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai) into ( mutect2somaticing, mutect2_contam, mantastrelka2ing, lanceting )
  tuple val(sampleID), val(germlineID) into sampleID_germID

  script:
  sampleID = "${tumourbam}".split("\\.")[0]
  germlineID = "${germlinebam}".split("\\.")[0]
  tumourbai = "${tumourbam}".replaceAll(".bam", ".bam.bai")
  germlinebai = "${germlinebam}".replaceAll(".bam", ".bam.bai")
  """
  #! bash
  samtools index ${tumourbam}
  samtools index ${germlinebam}
  """
}

// 1.4: scatter-gather implementation for mutect2, lancet
process scat_gath {

  label 'low_mem'

  input:
  file(intlist) from intlist_call1.getVal()

  output:
  file('lancet.scatgath.*.bed') into lancet_bedding
  file('mutect2.scatgath.*.bed.interval_list') into mutect2_bedding
  file('hc.scatgath.*.bed.interval_list') into hc_bedding

  script:
  def sgcount = 2
  """
  ##strip out all but chromosomes in the interval_list (no decoys etc)
  CHRS=\$(grep -v "@" ${intlist} | cut -f 1 | uniq)
  for CHR in \$CHRS; do
    grep "SN:\$CHR\\s" ${intlist} >> used.interval_list
  done
  grep -v "@" ${intlist} >> used.interval_list

  ##generate scatters
  picard IntervalListTools \
    I=used.interval_list \
    SCATTER_COUNT=${sgcount} \
    O=./

  ##rename scatters and parse into appropriate format for tools
  ls temp*/* | while read FILE; do
    COUNTN=\$(dirname \$FILE | perl -ane '@s=split(/\\_/); print \$s[1];');
    mv \$FILE mutect2.scatgath.\${COUNTN}.bed.interval_list;
    cp mutect2.scatgath.\${COUNTN}.bed.interval_list hc.scatgath.\${COUNTN}.bed.interval_list
    grep -v @ mutect2.scatgath.\${COUNTN}.bed.interval_list | \
      cut -f 1,2,3,5 > lancet.scatgath.\${COUNTN}.bed
  done
  """
}

mutect2bedding = mutect2_bedding.flatten()
mutect2somaticing
  .map { it -> [it[0],it[1],it[2],it[3],it[4],it[5]]}
  .combine(mutect2bedding)
  .set { mutect2somaticbedding }

// 2.7.1: MuTect2
// NB --germline-resource dollar-sign{dbsnp} removed as no AF causing error

process mutct2_sg {

  label 'med_mem'

  input:
  tuple val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai), file(intlist) from mutect2somaticbedding
  file(fasta) from reference.fa
  file(fai) from reference.fai
  file(dict) from reference.dict

  output:
  tuple val(sampleID), file('*sort.mutect2.vcf') into mutect2_gt
  tuple val(sampleID), file('*.vcf.stats') into mutect2_st
  tuple val(sampleID), file('*mutect2.f1r2.tar.gz') into mutect2_f1r2

  script:
  def taskmem = task.memory == null ? "" : "--java-options \"-Xmx" + javaTaskmem("${task.memory}") + "\""
  """
  SCATGATHN=\$(echo ${intlist} | perl -ane '@s=split(/\\./);print\$s[2];')
  gatk ${taskmem} \
    Mutect2 \
    --native-pair-hmm-threads ${task.cpus} \
    --reference ${fasta} \
    --input ${germlinebam} \
    --input ${tumourbam} \
    --normal-sample ${germlineID} \
    --tumor-sample ${sampleID} \
    --output ${sampleID}"."\${SCATGATHN}".mutect2.vcf" \
    --disable-sequence-dictionary-validation true \
    --f1r2-tar-gz \${SCATGATHN}".mutect2.f1r2.tar.gz" \
    -L ${intlist}

  picard SortVcf \
    I=${sampleID}"."\${SCATGATHN}".mutect2.vcf" \
    O=${sampleID}"."\${SCATGATHN}".sort.mutect2.vcf" \
    SD=${dict}
  """
}

// 2.7.2: MuTect2_merge
mutect2_gt
  .groupTuple()
  .map { it -> tuple(it[0], it[1][0..-1].flatten()) }
  .set { mutect2_fm }

process mutct2_concat {

  label 'med_mem'

  input:
  tuple val(sampleID), file(rawvcfs) from mutect2_fm

  output:
  tuple val(sampleID), file('*mutect2.merge.vcf') into mutect2_merge

  script:
  """
  ls *.sort.mutect2.vcf > vcf.list
  picard MergeVcfs I=vcf.list O=${sampleID}".mutect2.merge.vcf"
  """
}

mutect2_st
  .groupTuple()
  .map { it -> tuple(it[0], it[1][0..-1].flatten()) }
  .set { mutect2_sm }

// 2.7.3: MuTect2 Concatenate VCFs
process mutct2_concstat {

  label 'med_mem'

  input:
  tuple val(sampleID), file(stats) from mutect2_sm

  output:
  tuple val(sampleID), file('*mutect2.merge.vcf.stats') into mutect2_stats

  script:
  """
  STATS=\$(ls *stats | perl -ane 'foreach \$k (@F){print "--stats \$k ";}')
  gatk MergeMutectStats --output ${sampleID}".mutect2.merge.vcf.stats" \$STATS
  """
}

// 2.7.4: MuTect2 Concatenate VCFs
mutect2_f1r2.groupTuple()
            .map { it -> [it[0], it[1..-1].flatten()] }
            .set { mutect2_f1r2_set }

process mutct2_f1r2_comb {

  label 'med_mem'

  input:
  tuple val(sampleID), file(mutect2_ro) from mutect2_f1r2_set

  output:
  tuple val(sampleID), file("${sampleID}.mutect2.f1r2.tar.gz") into mutect2_f1r2_comb

  script:
  """
  ALL_F1R2_INPUT=\$(for x in *.mutect2.f1r2.tar.gz; do echo -n "-I \$x "; done)
  gatk LearnReadOrientationModel \$ALL_F1R2_INPUT -O ${sampleID}.mutect2.f1r2.tar.gz
  """
}

// 2.7.5: MuTect2 Contamination
mutect2_contam
  .join(mutect2_merge)
  .join(mutect2_stats)
  .join(mutect2_f1r2_comb)
  .groupTuple()
  .map { it -> [it[0], it[1..5].flatten(), it[6], it[7], it[8]].flatten() }
  .set { mutect2_contam_merge }

process mutct2_contam_filter {

  label 'med_mem'

  publishDir path: "${params.outDir}/samples/${sampleID}/mutect2", mode: "copy", overwrite: true

  input:
  tuple val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai), file(mergevcf), file(statsvcf), file(readorient) from mutect2_contam_merge
  file(fasta) from reference.fa
  file(fai) from reference.fai
  file(dict) from reference.dict
  file(gps_files) from reference.seqlevel
  file(intlist) from intlist_call3.getVal()

  output:
  tuple val(sampleID), file('*snv_indel.pass.vcf') into mutect2_pass
  tuple val(sampleID), file('*.raw.vcf') into mutect2_raw
  file('*') into completedmutect2call

  script:
  def taskmem = task.memory == null ? "" : "--java-options \"-Xmx" + javaTaskmem("${task.memory}") + "\""
  hg = params.assembly == "GRCh37" ? "hg19" : "hg38"
  gpsgz = params.seqlevel == "exome" ? "${gps_files}/${params.exomeTag}/af-only-gnomad.${params.exomeTag}.${hg}.noChr.vcf.gz" : "${gps_files}/af-only-gnomad.wgs.${hg}.noChr.vcf.gz"
  """
  gatk ${taskmem} \
    GetPileupSummaries \
    -I ${tumourbam} \
    -V ${gpsgz} \
    -O ${sampleID}".getpileupsummaries.table" \
    -L ${intlist}

  gatk CalculateContamination \
    -I ${sampleID}".getpileupsummaries.table" \
    -O ${sampleID}".calculatecontamination.table"

  CONTAM=\$(tail -n+2 ${sampleID}.calculatecontamination.table | cut -f 2 | cut -d "." -f 1)
  if [[ \$CONTAM != 0 ]]; then
    touch ${sampleID}".CONTAMINATION.WARNING.txt"
  fi

  gatk IndexFeatureFile \
    --input ${mergevcf}

  gatk ${taskmem} \
    FilterMutectCalls \
    --reference ${fasta} \
    --contamination-table ${sampleID}".calculatecontamination.table" \
    --interval-padding 5 \
    --output ${sampleID}".mutect2.FilterMutectCalls.vcf" \
    --unique-alt-read-count 3 \
    --variant ${mergevcf} \
    --stats ${statsvcf} \
    --disable-sequence-dictionary-validation true \
    --ob-priors ${readorient} \
    -L ${intlist}

  perl ${workflow.projectDir}/bin/filter_Lancet_Mutect2_Manta-Strelka2_Format.pl \
    ID=${sampleID} \
    DP=14 \
    MD=2 \
    VCF=${sampleID}".mutect2.FilterMutectCalls.vcf"
  """
}

// 2.9: Manta output is a pre-req for Strelka2, so call both here
process mntstr {

  label 'high_mem'

  publishDir path: "${params.outDir}/samples/${sampleID}/manta-strelka2", mode: "copy"

  input:
  tuple val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai) from mantastrelka2ing
  file(fasta) from reference.fa
  file(fai) from reference.fai
  file(dict) from reference.dict
  file(bed_files) from reference.seqlevel

  output:
  tuple val(sampleID), file("${sampleID}.strelka2.snv_indel.pass.vcf") into strelka2_pass
  tuple val(sampleID), file('*.raw.vcf') into strelka2_raw
  file('*.txt') into log_mantastrelka

  script:
  def bedgz = params.seqlevel == "wgs" ? "${bed_files}/wgs.bed.gz" : "${bed_files}/${params.exomeTag}/${params.exomeTag}.bed.gz"
  def callRegions = params.seqlevel == "exome" ? "--exome --callRegions ${bedgz}" : "--callRegions ${bedgz}"
  """
  {

    samtools view -H ${tumourbam} | grep SN | cut -f 2 | sed 's/SN://g' | grep -v M > r.fa
    samtools faidx ${fasta} -r r.fa >> new.fasta

    configManta.py ${callRegions} --referenceFasta=new.fasta --normalBam=${germlinebam} --tumourBam=${tumourbam} --runDir=manta

    manta/runWorkflow.py -m local -j ${task.cpus}

    configureStrelkaSomaticWorkflow.py ${callRegions} --referenceFasta=new.fasta --indelCandidates=manta/results/variants/candidateSmallIndels.vcf.gz --normalBam=${germlinebam} --tumorBam=${tumourbam} --runDir=strelka2

    strelka2/runWorkflow.py -m local -j ${task.cpus}

    ##merge into raw snv_indel
    gatk MergeVcfs -I strelka2/results/variants/somatic.snvs.vcf.gz -I strelka2/results/variants/somatic.indels.vcf.gz -O tmp.strelka2.snv_indel.vcf

    ${workflow.projectDir}/bin/manta_strelka2_rename_filter.sh  tmp.strelka2.snv_indel.vcf tmp2.strelka2.snv_indel.vcf ${sampleID} ${germlineID}

    perl ${workflow.projectDir}/bin/filter_Lancet_Mutect2_Manta-Strelka2_Format.pl \
        ID=${sampleID} \
        DP=14 \
        MD=2 \
        VCF=tmp2.strelka2.snv_indel.vcf

  } 2>&1 | tee > ${sampleID}.manta-strelka2.log.txt
  """
}

// 2.10.1: Lancet
lancetbedding = lancet_bedding.flatten()
lanceting
  .map { it -> [it[0],it[1],it[2],it[3],it[4],it[5]]}
  .combine(lancetbedding)
  .set { lancetsbedding }

process lancet_sg {

  label 'med_mem'

  input:
  tuple val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai), file(bed) from lancetsbedding
  file(fasta) from reference.fa
  file(fai) from reference.fai
  file(dict) from reference.dict

  output:
  tuple val(sampleID), file('*.sort.lancet.vcf') into lancet_gt

  when:
  params.seqlevel == "exome"

  script:
  scatgathn = "${bed}".split("\\.")[2]
  """
  lancet \
    --num-threads ${task.cpus} \
    --ref ${fasta} \
    --bed ${bed} \
    --tumor ${tumourbam} \
    --normal ${germlinebam} | \
    perl -ane 'if(\$F[0]=~m/^\\#CHROM/){
      \$_=~s/TUMOR/${sampleID}/;
      \$_=~s/NORMAL/${germlineID}/;
      print \$_;}
    else{print \$_;}' > ${sampleID}"."${scatgathn}".lancet.vcf"

  picard SortVcf \
    I=${sampleID}"."${scatgathn}".lancet.vcf" \
    O=${sampleID}"."${scatgathn}".sort.lancet.vcf" \
    SD=${dict}
  """
}

// 2.10.2: Lancet Merge
lancet_gt
  .groupTuple()
  .map { it -> tuple(it[0], it[1][0..-1].flatten()) }
  .set { lancet_fm }

process lancet_concat {

  label 'med_mem'

  input:
  tuple val(sampleID), file(rawvcf) from lancet_fm

  output:
  tuple val(sampleID), file('*lancet.merge.vcf') into lancet_merge

  script:
  """
  ls *.sort.lancet.vcf > vcf.list
  picard MergeVcfs I=vcf.list O=${sampleID}".lancet.merge.vcf"
  """
}

/* 2.10.2: Lancet Filter
*/
process lancet_filter {

  label 'med_mem'

  publishDir path: "${params.outDir}/samples/${sampleID}/lancet", mode: "copy"

  input:
  tuple val(sampleID), file(mergevcf) from lancet_merge

  output:
  tuple val(sampleID), file('*.pass.vcf') into lancet_pass
  tuple val(sampleID), file('*.raw.vcf') into lancet_raw
  file('*') into completedlancetcall

  script:
  """
  perl ${workflow.projectDir}/bin/filter_Lancet_Mutect2_Manta-Strelka2_Format.pl \
    ID=${sampleID} \
    DP=14 \
    MD=2 \
    VCF=${mergevcf}
  """
}

/*
================================================================================
                          3.  ANNOTATION AND REPORTING
================================================================================
*/

// 3.02: Annotate Vcfs
lancet_pass
  .join(strelka2_pass)
  .join(mutect2_pass)
  .set { run_vepann }

process vepann {

  label 'med_mem'

  publishDir path: "${params.outDir}/samples/${sampleID}/${caller}", mode: "copy", pattern: "${sampleID}.${caller}.snv_indel.pass.vep.vcf"

  input:
  tuple val(sampleID), file(lan_pass), file(str_pass), file(mut_pass) from run_vepann
  file(fasta) from reference.fa
  file(fai) from reference.fai
  file(dict) from reference.dict
  file(grchver) from reference.grchvers
  file(pcgrbase) from reference.pcgrbase

  output:
  tuple val(sampleID), file('*.vep.vcf') into runGRanges

  script:
  def grch_vers = "${grchver}".split("\\/")[-1]
  """
  for VCF in *vcf; do
    VCFO=\$(echo \$VCF | sed 's/vcf/vep.vcf/')
    vep --dir_cache ${pcgrbase}/data/${grch_vers}/.vep \
      --offline \
      --assembly ${params.assembly} \
      --vcf_info_field ANN \
      --symbol \
      --species homo_sapiens \
      --check_existing \
      --cache \
      --fork ${task.cpus} \
      --af_1kg \
      --af_gnomad \
      --vcf \
      --input_file \$VCF \
      --output_file \$VCFO \
      --format "vcf" \
      --fasta ${fasta} \
      --hgvs \
      --canonical \
      --ccds \
      --mane \
      --force_overwrite \
      --verbose
  done
  """
}

/* 3.1 RData GRanges from processed VCFs
* take publishDir and check for number of files therein
* each sample has 6 associated (raw, pass per caller)
* NB increment if adding callers!
*/

//separate out impacts processing as WGS cango above walltime
impacts = ["HIGH,MODERATE,MODIFIER,LOW"]

runGRanges
  .join(sampleID_germID)
  .join(lancet_raw)
  .join(strelka2_raw)
  .join(mutect2_raw)
  .set { runGRanges_in }

process vcfGRa {

  label 'max_mem'

  publishDir "${params.outDir}/combined/variant_consensus", mode: "copy"

  input:
  tuple val(sampleID), file(vcfs), val(germlineID), file(lan_raw), file(str_raw), file(mut_raw) from runGRanges_in
  each impact from impacts

  output:
  file('*impacts.pcgr.tsv.vcf') into vcfs_pcgr
  file('*') into completedvcfGRangesConsensus
  file('*.pdf') into sendmail_vcfGRa

  script:
  def inc_ord = params.incOrder ? params.incOrder : "noord"
  def which_genome =
  """
  Rscript -e "somenone::variant_consensus(germline_id = \\"${germlineID}\\", vep_vcf_pattern = \\"snv_indel.pass.vep.vcf\\", raw_vcf_pattern = \\"raw.vcf\\", tag = \\"${params.runID}\\", which_genome = \\"hg38\\", included_order = \\"${inc_ord}\\", impacts = \\"${impact}\\")"
  """
}

// 3.2 Create VCF for PCGR from consensus (reheader based on vcf 4.2 standards)
vcfs_pcgrd = vcfs_pcgr
              .collect()
              .flatten()

process pcgr_vcf {

  label 'low_mem'
  publishDir path: "${params.outDir}/samples/${sampleID}/VEP/"

  input:
  file(vcf) from vcfs_pcgrd

  output:
  tuple val(sampleID), file(ovcf) into snvpass_pcgr

  when:
  vcf =~ "HMML_impacts.pcgr.tsv.vcf"

  script:
  sampleID = "${vcf}".split("\\.")[0]
  ovcf = "${vcf}".replace("pcgr.tsv.vcf", "snv_indel.pass.vep.vcf")
  """
    vep --dir_cache ${pcgrbase}/data/${grch_vers}/.vep \
      --offline \
      --assembly ${params.assembly} \
      --vcf_info_field ANN \
      --symbol \
      --species homo_sapiens \
      --check_existing \
      --cache \
      --fork ${task.cpus} \
      --af_1kg \
      --af_gnomad \
      --vcf \
      --input_file ${vcf} \
      --output_file ${ovcf} \
      --format "vcf" \
      --fasta ${fasta} \
      --hgvs \
      --canonical \
      --ccds \
      --mane \
      --force_overwrite \
      --verbose
  """
}
