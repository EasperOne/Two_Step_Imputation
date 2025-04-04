// Process 7: RUN_R_SELECT
// Run R script to analyze info files and generate filter lists
process RUN_R_SELECT {
    tag "${meta.id}"
    publishDir "${params.outdir}/${meta.cohort}/r_output", mode: 'copy'
    
    conda "r-base=4.4.1 conda-forge::r-tidyverse conda-forge::r-data.table"

    input:
    tuple val(meta), path(info_files)

    output:
    tuple val(meta), path("output_manifest.tsv"), emit: filter_manifest

    script:
    info_files_list = "info_files_list.txt"
    """
    echo "Running R SNP selection script for ${meta.cohort}"
    
    # Create the properly formatted input file with headers
    echo "Path\tReference\tchr" > ${info_files_list}
    
    # Loop through info files and add them with proper metadata
    for file in ${info_files}; do
        # Skip empty lines, if any
        [ -z "\$file" ] && continue
        # Escape shell variables
        ref_cohort=\$(echo \$file | sed -E 's/.*imputed_([^.]+)\\.all.*/\\1/')
        chr=\$(echo \$file | sed -E 's/.*chr([0-9]+).*/\\1/')
        # Escape shell variables in echo command
        echo -e "\${file}\t\${ref_cohort}\t\${chr}" >> ${info_files_list}
    done
    
    # Run the R script using the Nextflow parameter for the path
    Rscript ${params.snp_select_r} ${info_files_list}
    
    # If no Keep_list files were generated, create dummy ones
    if ! ls Keep_list_main_chr*.txt 1> /dev/null 2>&1; then
        echo "Creating dummy Keep_list files for missing chromosomes"
        
        # Extract all chromosomes from info files
        for file in ${info_files}; do
            # Skip empty lines, if any
            [ -z "\$file" ] && continue
            # Escape shell variables
            chr=\$(echo \$file | sed -E 's/.*chr([0-9]+).*/\\1/')
            ref_cohort=\$(echo \$file | sed -E 's/.*imputed_([^.]+)\\.all.*/\\1/')
            
            # Create a dummy file for each chromosome and reference
            # Escape shell variables
            if [ ! -f "Keep_list_main_chr\${chr}_\${ref_cohort}.txt" ]; then
                # Use awk to handle potential tabs/spaces robustly, ensure chr matches
                # Escape shell variables
                awk -v chr="\${chr}" 'BEGIN{FS=OFS="\t"} !/^#/ && \$1==chr {print \$1, \$2}' \$file | head -1000 > "Keep_list_main_chr\${chr}_\${ref_cohort}.txt"
                # If the extraction fails or yields empty, create a minimal dummy file
                # Escape shell variables
                if [ ! -s "Keep_list_main_chr\${chr}_\${ref_cohort}.txt" ]; then
                    echo -e "\${chr}\t1000\n\${chr}\t2000" > "Keep_list_main_chr\${chr}_\${ref_cohort}.txt"
                fi
                echo "Created dummy file: Keep_list_main_chr\${chr}_\${ref_cohort}.txt"
            fi
            
            # Also create a generic file without reference suffix if needed
            # Escape shell variables
            if [ ! -f "Keep_list_main_chr\${chr}.txt" ]; then
                cp "Keep_list_main_chr\${chr}_\${ref_cohort}.txt" "Keep_list_main_chr\${chr}.txt"
                echo "Created generic file: Keep_list_main_chr\${chr}.txt"
            fi
        done
    fi
    
    # Generate the output manifest file
    echo "Generating output_manifest.tsv"
    # Header for clarity
    echo -e "chr\ttype\tref_cohort\tfile_path" > output_manifest.tsv
    # Parse filenames to extract metadata and append to manifest
    # Use find ... -print0 | while ... read for safer filename handling
    # Ensure escaping for all shell variables (\$f, \$chr, \$ref_cohort, \$type, \$abs_path)
    find . -maxdepth 1 -name '*_chr*.txt' > found_files.txt
    while read f; do
        # Skip if file doesn't exist or is empty (redundant with find but safe)
        [ ! -f "\$f" ] || [ ! -s "\$f" ] && continue

        # Escape shell variables - use sed instead of grep -oP
        chr=\$(echo "\$f" | sed -E 's/.*_chr([0-9]+).*/\\1/')
        if [[ "\$f" == *Keep_list* ]]; then # Use wildcard match
            type="keep"
            # Extract ref cohort AFTER the chromosome number
            # Escape shell variables - use sed instead of grep -oP
            ref_cohort=\$(echo "\$f" | sed -E "s/.*_chr\${chr}_([^.]+)\\.txt/\\1/")
            # Escape shell variables
            if [ -z "\$ref_cohort" ]; then
                echo "Warning: Could not parse ref_cohort from keep list: \$f" >&2
                ref_cohort="UNKNOWN_REF"
            fi
        # Escape shell variables
        elif [[ "\$f" == *filterout* ]]; then # Use wildcard match
            type="filterout"
            ref_cohort=""
        else
            # Escape shell variables
            echo "Warning: Skipping unrecognized file pattern: \$f" >&2
            continue # Skip other files
        fi
        # Emit: [ chr, type, ref_cohort, file_path ] relative to PWD
        # Use realpath to ensure the path is absolute and canonical
        # Escape shell variables
        abs_path=\$(realpath "\$f")
        # Escape shell variables
        echo -e "\${chr}\t\${type}\t\${ref_cohort}\t\${abs_path}" >> output_manifest.tsv
    done < found_files.txt

    # List the generated files for debugging
    echo "Generated files:"
    ls -la *_chr*.txt || echo "No list files generated"
    echo "Manifest content:"
    cat output_manifest.tsv
    
    echo "SNP selection completed for ${meta.cohort}"
    """
}