class: Workflow
description: "immune_repertoire_pipeline"
cwlVersion: cwl:draft-3

inputs:
  - id: r1_fqs
    type:
      type: array
      items: File
  - id: r2_fqs
    type:
      type: array
      items: File
  - id: type
    type:
      type: array
      items: string
  - id: sampleName
    type:
      type: array
      items: string
  - id: Diff
    type:
      type: array
      items: string
  - id: fold
    type: string
  - id: Thred
    type: string
  - id: Project_name
    type: string
  - id: Project_key
    type: string
  - id: Species
    type: string
  - id: Genome_version
    type: string
  - id: refData
    type: string
  - id: resolution
    type: string
  - id: minexp
    type: string
  - id: mincell
    type: string
  - id: minUMI
    type: string
  - id: minGene
    type: string
  - id: maxGene
    type: string
  - id: maxpct
    type: string
  - id: SCT
    type: string

steps: 
########################################################
  - id: creat_detail_cfg
    run:
      class: CommandLineTool
      description: "..."
      hints:
        - class: DockerRequirement
          dockerImageId: immune_repertoire_pipeline
        - class: ResourceRequirement
          coresMax: 1
          ramMax: 0
          baseRam: 403
          modulus: 0 
      inputs:
        - id: Project_name
          type: string
          inputBinding:
            prefix: --Project_name
        - id: Project_key
          type: string
          inputBinding:
            prefix: --Project_key
        - id: Diff
          type:
            type: array
            items: string
          inputBinding:
            prefix: --Diff
        - id: Thred
          type: string
          inputBinding:
            prefix: --Thred
        - id: fold
          type: string
          inputBinding:
            prefix: --fold
        - id: Species
          type: string
          inputBinding:
            prefix: --Species
        - id: Genome_version
          type: string
          inputBinding:
            prefix: --Genome_version
        - id: refData
          type: string
          inputBinding:
            prefix: --refData
        - id: minexp
          type: string
          inputBinding:
            prefix: --minexp
        - id: mincell
          type: string
          inputBinding:
            prefix: --mincell
        - id: minUMI
          type: string
          inputBinging:
            prefix: --minUMI
        - id: minGene
          type: string
          inputBinding:
            prefix: --minGene
        - id: maxGene
          type: string
          inputBinding:
            prefix: --maxGene
        - id: maxpct
          type: string
          inputBinding:
            prefix: --maxpct
        - id: resolution
          type: string
          inputBinding:
            prefix: --resolution
        - id: SCT
          type: string
          inputBinding: 
            prefix: --SCT
      outputs:
        - id: detail_cfg
          type: File
          outputBinding:
            glob: ^detail\.cfg$
        - id: deg_cfg
          type: File
          outputBinding:
            glob: ^deg_cfg\.cfg$
        - id: group_file
          type: 
            type: array
            items: File
          outputBinding:
            glob: ^*\.txt$
        - id: group_name_stat
          type: File
          outputBinding:
            glob: ^Group_info\.stat$
      baseCommand: creat_detail.pl
    inputs:
      - { id: Project_name, source: "Project_name" }
      - { id: Project_key, source: "Project_key" }
      - { id: Diff, source: "Diff", nonSymbol: true }
      - { id: Thred, source: "Thred" }
      - { id: fold, source: "fold" }
      - { id: Species, source: "Species" }
      - { id: Genome_version, source: "Genome_version" }
      - { id: minexp, source: "minexp" }
      - { id: mincell, source: "mincell" }
      - { id: minUMI, source: "minUMI" }
      - { id: minGene, source: "minGene" }
      - { id: maxGene, source: "maxGene", nonSymbol: true }
      - { id: maxpct, source: "maxpct" }
      - { id: resolution, source: "resolution" }
      - { id: SCT, source: "SCT" }
      - { id: refData, source: "refData", nonSymbol: true }
    outputs:
      - id: detail_cfg
        type: File
      - id: deg_cfg
        type: File
      - id: group_file
        type: 
          type: array
          items: File
        nonSymbol: True
      - id: group_name_stat
        type: File
        nonSymbol: True
##############################################################
  - id: Phred_Change
    run:
      class: CommandLineTool
      description: "..."
      hints:
        - class: DockerRequirement
          dockerImageId: immune_repertoire_pipeline
        - class: ResourceRequirement
          coresMax: 1
          ramMax: 0
          baseRam: 100
          modulus: 0
      inputs:
        - id: r1_fqs
          type:
            type: array
            items: File
          inputBinding:
            prefix: --r1_fqs
        - id: r2_fqs
          type:
            type: array
            items: File
          inputBinding:
            prefix: --r2_fqs
        - id: sample
          type:
            type: array
            items: string
          inputBinding:
            prefix: --sample
        - id: type
          type:
            type: array
            items: string
          inputBinding:
            prefix: --type
      outputs:
        - id: data_cfg
          type: File
          outputBinding:
            glob: ^*\.cfg$
        - id: samples_name
          type: File
          outputBinding:
            glob: ^samples_name\.xls$
        - id: samples_fq_dir
          type:
            type: array
            items: directory
          outputBinding:
            glob: ^*\.fastq$
      baseCommand: Phred_Change.pl
    inputs:
      - { id: r1_fqs, source: "r1_fqs" }
      - { id: r2_fqs, source: "r2_fqs" }
      - { id: sample, source: "sampleName" }
      - { id: type, source: "type" }
    outputs:
      - id: data_cfg
        type: File
      - id: samples_name
        type: File
        nonSymbol: true
      - id: samples_fq_dir
        type:
          type: array
          items: directory
##############################################################
################################################################
  - id: data_access_fastqc
    run:
      class: CommandLineTool
      description: "..."
      hints:
        - class: DockerRequirement
          dockerImageId: immune_repertoire_pipeline
        - class: ResourceRequirement
          coresMax: 1
          ramMax: 0
          baseRam: 6000
          modulus: 0.02
      inputs:
        - id: sample_fq_dir
          type: directory
          inputBinding:
            prefix: --sample_fq_dir
      outputs:
        - id: fastqc_result
          type:
            type: array
            items: directory
          outputBinding:
            glob: ^*\_fastqc$
      baseCommand: fastqc_data_access.pl
    inputs:
      - { id: sample_fq_dir, source: "Phred_Change.samples_fq_dir" }
    outputs:
      - id: fastqc_result
        type: 
          type: array
          items: directory
##############################################################
  - id: data_access
    run:
      class: CommandLineTool
      description: "..."
      hints:
        - class: DockerRequirement
          dockerImageId: immune_repertoire_pipeline
        - class: ResourceRequirement
          coresMax: 1
          ramMax: 0
          baseRam: 6000
          modulus: 0.02
      inputs:
        - id: single_sample_data
          type: directory
          inputBinding:
            prefix: --sample_dir
      outputs:
        - id: fq1_file
          type: 
            type: array
            items: File
          outputBinding:
            glob: ^*_R1_001\.fastq$
        - id: fq2_file
          type:
            type: array
            items: File
          outputBinding:
            glob: ^*_R2_001\.fastq$    
        - id: sample_stat
          type:
            type: array
            items: File
          outputBinding:
            glob: ^*\.stat$
      baseCommand: data_access.pl
    inputs:
      - { id: single_sample_data, source: "Phred_Change.samples_fq_dir" }
    outputs:
      - id: fq1_file
        type:
          type: array
          items: File
      - id: fq2_file
        type:
          type: array
          items: File
      - id: sample_stat
        type:
          type: array
          items: File
################################################################
  - id: data_access_stat
    run:
      class: CommandLineTool
      description: "..."
      hints:
        - class: DockerRequirement
          dockerImageId: immune_repertoire_pipeline
        - class: ResourceRequirement
          coresMax: 1
          ramMax: 0
          baseRam: 100
          modulus: 0.02
      inputs:
        - id: sample_stats
          type: 
            type: array
            items: File
          inputBinding:
            prefix: --sample_stats
      outputs:
        - id: All_sample_stat
          type: File
          outputBinding:
            glob: ^AllSample_GC_Q\.stat$
      baseCommand: data_access_stat.pl
    inputs:
      - { id: sample_stats, source: "data_access.sample_stat" }
    outputs:
      - id: All_sample_stat
        type: File
################################################################
  - id: cellranger
    run: 
      class: CommandLineTool
      description: "..."
      hints: 
        - class: DockerRequirement
          dockerImageId: immune_repertoire_pipeline
        - class: ResourceRequirement
          coresMax: 8
          ramMax: 0
          baseRam: 65000
          modulus: 1
      inputs:
        - id: sample_fq_dir
          type: directory
          inputBinding: 
            prefix: --fqdir
        - id: detail_cfg
          type: File
          inputBinding:
            prefix: --detail
      outputs:
        - id: samplesResult
          type:
            type: array
            items: directory
          outputBinding:
            glob: ^*\.Result$
        - id: summary_file
          type:
            type: array
            items: File
          outputBinding:
            glob: ^*\.metrics_summary\.csv$
        - id: Analysis_Result
          type: 
            type: array
            items: directory
          outputBinding:
            glob: ^*\.analysis$
        - id: origin_results
          type: 
            type: array
            items: directory
          outputBinding:
            glob: ^*\.origin_results$
      baseCommand: cellranger.pl
    inputs:
      - { id: sample_fq_dir, source: "Phred_Change.samples_fq_dir" }
      - { id: detail_cfg, source: "creat_detail_cfg.detail_cfg" }
    outputs:
      - id: samplesResult
        type:
          type: array
          items: directory
      - id: Analysis_Result
        type:
          type: array
          items: directory
      - id: summary_file
        type:
          type: array
          items: File
      - id: origin_results
        type:
          type: array
          items: directory

#################################################################
  - id: cellranger_stat
    run:
      class: CommandLineTool
      description: "..."
      hints:
        - class: DockerRequirement
          dockerImageId: immune_repertoire_pipeline
        - class: ResourceRequirement
          coresMax: 1
          ramMax: 0
          baseRam: 100
          modulus: 0
      inputs:
        - id: summary_file
          type:
            type: array
            items: File
          inputBinding:
            prefix: --summary_file
      outputs:
        - id: statistic
          type:
            type: array
            items: directory
          outputBinding:
            glob: ^*statistic$
      baseCommand: cellranger_stat.pl
    inputs:
      - { id: summary_file, source: "cellranger.summary_file" }
    outputs:
      - id: statistic
        type:
          type: array
          items: directory
#################################################################
  - id: sc_basic_analysis
    run:
      class: CommandLineTool
      description: "..."
      hints:
        - class: DockerRequirement
          dockerImageId: immune_repertoire_pipeline
        - class: ResourceRequirement
          coresMax: 13
          ramMax: 0
          baseRam: 35000
          modulus: 0
      inputs:
        - id: config
          type: File
          inputBinding:
            prefix: --config
        - id: sample_cellranger_result
          type:
            type: array
            items: directory
          inputBinding:
            prefix: --sample_cellranger_result
      outputs:
        - id: filtered
          type: directory
          outputBinding: 
           glob: ^filtered$
        - id: analysed
          type: directory
          outputBinding: 
           glob: ^analysed$
        - id: analysed_integrated
          type: directory
          outputBinding: 
           glob: ^analysed_integrated$
        - id: cluster_diff_integrated
          type: directory
          outputBinding: 
           glob: ^cluster_diff_integrated$
        - id: sample_diff_integrated
          type: directory
          outputBinding: 
           glob: ^sample_diff_integrated$
        - id: sample_statistic
          type:
           type: array
           items: directory
          outputBinding:
           glob: ^*\.statistic$
        - id: rds
          type:
           type: array
           items: File
          outputBinding:
           glob: ^*\.Rds$
        - id: inte_rds
          type: File
          outputBinding:
           glob: ^analysed_integrated\.single_seruat\.Rds$
        - id: rda
          type:
           type: array
           items: File
          outputBinding:
           glob: ^*\.Rda$
        - id: avg
          type:
           type: array
           items: File
          outputBinding:
           glob: ^*All_cluster_Markergene_avgExp\.xls$
        - id: id_list
          type: File
          outputBinding:
           glob: ^symbol.list$
      baseCommand: sc_basic_analysis.py
    inputs:
      - { id: config, source: "creat_detail_cfg.detail_cfg"}
      - { id: sample_cellranger_result, source: "cellranger.samplesResult"}
    outputs:
      - id: filtered
        type: directory
      - id: analysed
        type: directory
      - id: analysed_integrated
        type: directory
        nonSymbol: True
      - id: cluster_diff_integrated
        type: directory
        nonSymbol: True
      - id: sample_diff_integrated
        type: directory
        nonSymbol: True
      - id: sample_statistic
        type:
         type: array
         items: directory
      - id: rds
        type:
         type: array
         items: File
      - id: rda
        type:
         type: array
         items: File
      - id: inte_rds
        type: File
        nonSymbol: True
      - id: id_list
        type: File
      - id: avg
        type:
         type: array
         items: File
#################################################################
  - id: sc_gene_functional_analysis
    run:
      class: CommandLineTool
      description: "..."
      hints:
        - class: DockerRequirement
          dockerImageId: immune_repertoire_pipeline
        - class: ResourceRequirement
          coresMax: 5
          ramMax: 0
          baseRam: 6000
          modulus: 0
      inputs:
        - id: config
          type: File
          inputBinding:
            prefix: --config
        - id: id_list
          type: File
          inputBinding:
            prefix: --id_list
        - id: statistic_file
          type:
            type: array
            items: directory
          inputBinding:
            prefix: --statistic_file
        - id: all_cluster_marker_avg
          type:
            type: array
            items: File
          inputBinding:
            prefix: --all_cluster_marker_avg
      outputs:
        - id: Anno_enrichment
          type:
           type: array
           items: directory
          outputBinding: 
           glob: ^*\.Anno_enrichment$
        - id: ppi_result
          type:
           type: array
           items: directory
          outputBinding:
           glob: ^*\.ppi_result$
        - id: TFBS_Analysis
          type:
           type: array
           items: directory
          outputBinding:
           glob: ^*\.TFBS_Analysis$
      baseCommand: sc_gene_functional_analysis.py
    inputs:
      - { id: config, source: "creat_detail_cfg.detail_cfg"}
      - { id: id_list, source: "sc_basic_analysis.id_list"}
      - { id: statistic_file, source: "sc_basic_analysis.sample_statistic"}
      - { id: all_cluster_marker_avg, source: "sc_basic_analysis.avg"}
    outputs:
      - id: Anno_enrichment
        type:
         type: array
         items: directory
      - id: ppi_result
        type:
         type: array
         items: directory
      - id: TFBS_Analysis
        type:
         type: array
         items: directory
#################################################################
  - id: sc_character_analysis
    run:
      class: CommandLineTool
      description: "..."
      hints:
        - class: DockerRequirement
          dockerImageId: immune_repertoire_pipeline
        - class: ResourceRequirement
          coresMax: 4
          ramMax: 0
          baseRam: 35000
          modulus: 0
      inputs:
        - id: rds
          type:
            type: array
            items: File
          inputBinding:
            prefix: --rds
        - id: config
          type: File
          inputBinding:
            prefix: --config
      outputs:
        - id: cell_cycle_results
          type:
           type: array
           items: directory
          outputBinding: 
           glob: ^*\.cell_cycle$
        - id: cell_typeAnno_results
          type:
           type: array
           items: directory
          outputBinding: 
           glob: ^*\.cell_typeAnno$
        - id: cell_trace_results
          type:
           type: array
           items: directory
          outputBinding: 
           glob: ^*\.cell_trace$
      baseCommand: sc_character_analysis.py
    inputs:
      - { id: config, source: "creat_detail_cfg.detail_cfg"}
      - { id: rds, source: "sc_basic_analysis.rds"}
    outputs:
      - id: cell_cycle_results
        type:
          type: array
          items: directory
      - id: cell_typeAnno_results
        type:
          type: array
          items: directory
      - id: cell_trace_results
        type:
          type: array
          items: directory
        nonSymbol: True
#################################################################
  - id: vdj_analysis
    run:
      class: CommandLineTool
      description: "..."
      hints:
        - class: DockerRequirement
          dockerImageId: immune_repertoire_pipeline
        - class: ResourceRequirement
          coresMax: 2
          ramMax: 0
          baseRam: 10000
          modulus: 0
      inputs:
        - id: rda_files
          type:
            type: array
            items: File
          inputBinding:
            prefix: --rda_files
        - id: sample_cellranger_result
          type:
            type: array
            items: directory
          inputBinding:
            prefix: --sample_cellranger_result
        - id: config
          type: File
          inputBinding:
            prefix: --config
      outputs:
        - id: b_results
          type: directory
          outputBinding: 
           glob: ^*b_results$
        - id: t_results
          type: directory
          outputBinding: 
           glob: ^*t_results$
      baseCommand: vdj_analysis.py
    inputs:
      - { id: config, source: "creat_detail_cfg.detail_cfg"}
      - { id: rda_files, source: "sc_basic_analysis.rda"}
      - { id: sample_cellranger_result, source: "cellranger.samplesResult"}
    outputs:
      - id: b_results
        type: directory
        nonSymbol: True
      - id: t_results
        type: directory
        nonSymbol: True
#################################################################
  - id: web_report
    run:
      class: CommandLineTool
      description: "..."
      hints:
        - class: DockerRequirement
          dockerImageId: immune_repertoire_pipeline
        - class: ResourceRequirement
          coresMax: 1
          ramMax: 0
          baseRam: 6000
          modulus: 0
      inputs:
        - id: gc
          type: File
          inputBinding:
            prefix: --gc
        - id: qc
          type:
           type: array
           items: directory
          inputBinding:
            prefix: --qc
        - id: sample_filter
          type: directory
          inputBinding:
            prefix: --sample_filter
        - id: sample_analysis
          type: directory
          inputBinding:
            prefix: --sample_analysis
        - id: single_enrichment
          type:
           type: array
           items: directory
          inputBinding:
            prefix: --single_enrichment
        - id: single_sample_ppi
          type:
           type: array
           items: directory
          inputBinding:
            prefix: --single_sample_ppi
        - id: single_tf_analysis
          type:
           type: array
           items: directory
          inputBinding:
            prefix: --single_tf_analysis
        - id: single_typeanno
          type:
           type: array
           items: directory
          inputBinding:
            prefix: --single_typeanno
        - id: integrated_result
          type: directory
          inputBinding:
            prefix: --integrated_result
        - id: allcluster_statistic
          type: directory
          inputBinding:
            prefix: --allcluster_statistic
        - id: groupdiff_statistic
          type: directory
          inputBinding:
            prefix: --groupdiff_statistic
        - id: cell_cycle
          type:
           type: array
           items: directory
          inputBinding:
            prefix: --cell_cycle
        - id: allsample_trace
          type:
           type: array
           items: directory
          inputBinding:
            prefix: --allsample_trace
        - id: cfg1
          type: File
          inputBinding:
            prefix: --cfg1
        - id: cfg2
          type: File
          inputBinding:
            prefix: --cfg2
        - id: cellranger_stat
          type:
           type: array
           items: directory
          inputBinding:
            prefix: --cellranger_stat
        - id: cellranger_results
          type:
           type: array
           items: directory
          inputBinding:
            prefix: --cellranger
        - id: t_results
          type: directory
          inputBinding:
            prefix: --t_results
        - id: b_results
          type: directory
          inputBinding:
            prefix: --b_results
      outputs:
        - id: Web_Report
          type: directory
          outputBinding:
           glob: ^Web_Report$
        - id: Needed_Data
          type: directory
          outputBinding:
           glob: ^Needed_Data$
      baseCommand: web_report.py
    inputs:
      - { id: gc, source: "data_access_stat.All_sample_stat" }
      - { id: qc, source: "data_access_fastqc.fastqc_result" }
      - { id: sample_filter, source: "sc_basic_analysis.filtered" }
      - { id: sample_analysis, source: "sc_basic_analysis.analysed" }
      - { id: single_enrichment, source: "sc_gene_functional_analysis.Anno_enrichment" }
      - { id: single_sample_ppi, source: "sc_gene_functional_analysis.ppi_result" }
      - { id: single_tf_analysis, source: "sc_gene_functional_analysis.TFBS_Analysis" }
      - { id: single_typeanno, source: "sc_character_analysis.cell_typeAnno_results" }
      - { id: allcluster_statistic, source: "sc_basic_analysis.cluster_diff_integrated", nonSymbol: true }
      - { id: groupdiff_statistic, source: "sc_basic_analysis.sample_diff_integrated", nonSymbol: true }
      - { id: allsample_trace, source: "sc_character_analysis.cell_trace_results", nonSymbol: true}
      - { id: cfg1, source: "Phred_Change.data_cfg" }
      - { id: cfg2, source: "creat_detail_cfg.detail_cfg" }
      - { id: cellranger_stat, source: "cellranger_stat.statistic" }
      - { id: t_results, source: "vdj_analysis.t_results", nonSymbol: true }
      - { id: b_results, source: "vdj_analysis.b_results", nonSymbol: true }
      - { id: integrated_result, source: "sc_basic_analysis.analysed_integrated", nonSymbol: true}
      - { id: cellranger_results, source: "cellranger.origin_results" }
      - { id: cell_cycle, source: "sc_character_analysis.cell_cycle_results" }
    outputs:
      - id: Web_Report
        type: directory
      - id: Needed_Data
        type: directory