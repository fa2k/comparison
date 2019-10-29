#multiqc -o 20_piccard/ -f --cl_config "extra_fn_clean_exts: [{type: 'regex', pattern: '_[A-Z0-9]{9}_L00[1-8]'}]" 20_piccard/
multiqc -o 30_downsample/ -f 30_downsample/
