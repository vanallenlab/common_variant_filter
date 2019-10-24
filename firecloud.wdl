workflow CommonVariantFilter {
    String sampleId
    File maf

    Int? min_exac_ac = 10
    Int? min_filter_depth = 0
    Boolean? disable_whitelist = false
    Boolean? disable_comment = false
    Boolean? filter_noncoding = false

    String? disable_whitelist_str = if disable_whitelist == true then 'true' else 'false'
    String? disable_comment_str = if disable_comment == true then 'true' else 'false'
    String? filter_noncoding_str = if filter_noncoding == true then 'true' else 'false'

    Int? RAM = 4
    Int? SSD = 25
    Int? preemptible = 3

    String? docker_tag = "1.0.3"

    meta {
        author: "Brendan Reardon"
        email: "breardon@broadinstitute.org"
        laboratory: "Van Allen Lab"
        institution: "Dana-Farber Cancer Institute, Broad Institute of MIT & Harvard"
        github: "https://github.com/vanallenlab/common_variant_filter"
        license: "MIT License"
    }

    call commonfilterTask {
        input: sampleId=sampleId,
            maf=maf,
            min_exac_ac=min_exac_ac,
            min_filter_depth=min_filter_depth,
            disable_whitelist=disable_whitelist_str,
            filter_noncoding=filter_noncoding_str,
            RAM=RAM,
            SSD=SSD,
            preemptible=preemptible,
            docker_tag=docker_tag
    }

    output  {
        File commonfilter_annotated_maf = commonfilterTask.annotatedMAF
        File commonfilter_passed_maf = commonfilterTask.passedMAF
        File commonfilter_rejected_maf = commonfilterTask.rejectedMAF
        Int commonfilter_considered_count = commonfilterTask.consideredCount
        Int commonfilter_pass_count = commonfilterTask.passCount
        Int commonfilter_reject_count = commonfilterTask.rejectCount
    }
}

task commonfilterTask {
    String sampleId
    File maf

    Int? min_exac_ac
    Int? min_filter_depth
    String? disable_whitelist
    String? disable_comment
    String? filter_noncoding

    Int? RAM
    Int? SSD
    Int? preemptible

    String? docker_tag

    command {
        args="--min_exac_ac "${min_exac_ac}" "
        args+="--min_filter_depth "${min_filter_depth}" "
        if [ ${filter_noncoding} == 'true' ];
            then args+="--filter_noncoding "; fi
        if [ ${disable_whitelist} == 'true'];
            then args+="--disable_wl "; fi
        if [ ${disable_comment} == 'true'];
            then args+="--disable_comment "; fi

        sed '/^#/ d' < ${maf} > uncommented_maf.maf

        python /common_variant_filter.py --id ${sampleId} --maf uncommented_maf.maf $args
    }

    output  {
        File passedMAF="${sampleId}.common_variant_filter.pass.maf"
        File rejectedMAF="${sampleId}.common_variant_filter.reject.maf"
        String consideredCount=read_string("considered.txt")
		String passCount=read_string("passed.txt")
		String rejectCount=read_string("rejected.txt")
    }

    runtime {
        disks: "local-disk " + SSD + " SSD"
        docker: "vanallenlab/common_variant_filter:" + docker_tag
        memory: RAM + " GB"
        preemptible: preemptible
    }
}
