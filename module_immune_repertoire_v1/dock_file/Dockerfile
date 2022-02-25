FROM hub.bmkcloud.lan/base2:v1.0

##set software
        ENV PATH /share/nas2/genome/biosoft/STAR/current/bin/:/share/nas2/genome/biosoft/samtools/1.9/bin/:/share/nas2/genome/biosoft/cellranger/6.0.0/bin:/share/nas2/genome/biosoft/fastqc/current:/share/nas2/genome/bmksoft/tool/svg2xxx/v1.0/:$PATH
        ENV PATH /share/nas2/genome/biosoft/R/4.0.3/bin/:$PATH
        ENV PATH /share/nas2/genome/bmksoft/tool/fastq_qc_stat/v1.1:$PATH
        ENV LD_LIBRARY_PATH /share/nas2/genome/biosoft/gcc/5.5.0/lib/:/share/nas2/genome/biosoft/gcc/5.5.0/lib64/:/share/nas2/genome/biosoft/gmp/4.3.2/lib/:/share/nas2/genome/biosoft/mpc/1.0.1/lib/:/share/nas2/genome/biosoft/mpfr/2.4.2/lib/:/share/nas2/genome/biosoft/isl/0.14/lib/:/share/nas2/genome/biosoft/curl/7.54.1/lib/:/share/nas2/genome/biosoft/curl/7.47.1/lib/:/share/nas2/genome/biosoft/pcre/lib/:/share/nas2/genome/biosoft/pcre2/10.34/lib/:/share/nas2/genome/biosoft/Anaconda3/2019.03/lib/:/share/nas2/genome/biosoft/libicuuc:/share/nas2/genome/biosoft/gsl/2.6/lib/:$LD_LIBRARY_PATH
        ENV PATH /share/nas2/genome/biosoft/Python/3.6.5/bin:$PATH

##set script 为了防止脚本重名，一定要在本地将下面的脚本软连接到一个目录下，检查
        ENV PATH /share/nas2/genome/cloud_soft/developer_platform/module_script/module_immune_repertoire/trunk/:$PATH
        ENV PATH /share/nas2/genome/cloud_soft/developer_platform/module_script/module_immune_repertoire/trunk/DataAssess:$PATH
        ENV PATH /share/nas2/genome/cloud_soft/developer_platform/module_script/module_immune_repertoire/trunk/ScBasicAnalysis:$PATH
        ENV PATH /share/nas2/genome/cloud_soft/developer_platform/module_script/module_immune_repertoire/trunk/VDJ:$PATH
        ENV PATH /share/nas2/genome/cloud_soft/developer_platform/module_script/module_immune_repertoire/trunk/cell_cycle:$PATH
        ENV PATH /share/nas2/genome/cloud_soft/developer_platform/module_script/module_immune_repertoire/trunk/module_enrich:$PATH
        ENV PATH /share/nas2/genome/cloud_soft/developer_platform/module_script/module_immune_repertoire/trunk/PPI:$PATH
        ENV PATH /share/nas2/genome/cloud_soft/developer_platform/module_script/module_immune_repertoire/trunk/TF_analysis:$PATH
        ENV PATH /share/nas2/genome/cloud_soft/developer_platform/module_script/module_immune_repertoire/trunk/WebReport:$PATH
        ENV PATH /share/nas2/genome/cloud_soft/developer_platform/module_script/module_immune_repertoire/trunk/cellranger:$PATH

