#!/usr/bin/env nextflow

nextflow.preview.output = true

include { GZIP_DECOMPRESS as GZIP_DECOMPRESS_FASTA } from './modules/gzip'
include { RESFINDER                                } from './modules/cgetools/resfinder'
include { RESFINDER     as RESFINDER_LONGREAD      } from './modules/cgetools/resfinder'
include { RESFINDER     as RESFINDER_SHORTREAD     } from './modules/cgetools/resfinder'
include { PLASMIDFINDER                            } from './modules/cgetools/plasmidfinder'
include { PLASMIDFINDER as PLASMIDFINDER_LONGREAD  } from './modules/cgetools/plasmidfinder'
include { PLASMIDFINDER as PLASMIDFINDER_SHORTREAD } from './modules/cgetools/plasmidfinder'
include { ORGFINDER_DETECT                         } from './modules/orgfinder/detect'
include { SPECIATOR                                } from './modules/speciator'
include { AMRFINDERPLUS_UPDATE                     } from './modules/amrfinderplus/update'
include { AMRFINDERPLUS_RUN                        } from './modules/amrfinderplus/run'
include { PROKKA_RUN                               } from './modules/tseemann/prokka'
include { MLST                                     } from './modules/tseemann/mlst'
include { MLST as CGEMLST                          } from './modules/cgetools/mlst'
include { MOBTYPER_RUN                             } from './modules/mobsuite/mobtyper'
include { SAMTOOLS_FAIDX                           } from './modules/samtools/faidx'
include { SAMTOOLS_FASTQ                           } from './modules/samtools/fastq'
include { TO_JSON                                  } from './modules/tojson'
include { MULTIREPORT       } from './subworkflows/multireport'
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'



def get_samplesheet() {
	ss = [
		asm_ch: Channel.empty(),
		lr_ch: Channel.empty(),
		sr_ch: Channel.empty()
	]
	if (params.samplesheet) {
		SS = Channel.fromList(samplesheetToList(params.samplesheet, "assets/schema_samplesheet.json"))
			.multiMap({x ->
				asm_ch: [x[0].subMap('sample_id','assembly_name'),x[0].assembly_fasta]
				lr_ch: [x[0].subMap('sample_id','assembly_name'),x[0].long_reads]
				sr_ch: [x[0].subMap('sample_id','assembly_name'),[x[0].short_reads_1,x[0].short_reads_2]]
			})
		ss.lr_ch = SS.lr_ch
		ss.sr_ch = SS.sr_ch
		ss.asm_ch = SS.asm_ch
	} else {
		if (params.long_reads) {
			ss.lr_ch = Channel.fromPath(params.long_reads)
				.map({x -> tuple([sample_id:x.name.replaceAll(/\.(fastq\.gz|fq\.gz|bam|cram)$/,'')],x)})
		}
		if (params.short_reads) {
			ss.sr_ch = Channel
				.fromFilePairs(params.short_reads,size:-1) { file -> file.name.replaceAll(/_(R?[12])(_001)?\.(fq|fastq)\.gz$/, '') }
				.map({id,x -> [[sample_id:id],x]})
		}
		if (params.assembly_fasta) {
			ss.asm_ch = Channel.fromPath(params.assembly_fasta)
				.map({x -> tuple([sample_id:x.name.replaceAll(/\.(fasta|fa|fna)(\.gz)?$/,'')],x)})
		}
	}
	// Filter missing values
	ss.asm_ch = ss.asm_ch.filter({x,y -> y})
	ss.sr_ch = ss.sr_ch.map({x,y -> [x,y.findAll({v->v})]}).filter({x,y -> y})
	ss.lr_ch = ss.lr_ch.filter({x,y -> y})
	return ss
}


// Get name of the organism to use for given sample
def org_name(meta) {
	if (params.org_name!=null) return params.org_name
    if (meta.containsKey('org_name')) return meta['org_name']
	return null
}

// Retreive arguments to use for the given tool on the given sample
def tool_args(tool_name,meta,org_name=null) {
	def key = tool_name + "_args"
	def default_args_key = 'default_' + key
	if (params.containsKey(key)) return params[key]
    if (meta.containsKey(key)) return meta[key]
    if (org_name==null) return params[default_args_key]
    def org_args = params.organisms.containsKey(org_name)?params.organisms[org_name] : [:]
    if (org_args.containsKey(key)) return org_args[key]
    return params[default_args_key]
}



workflow ANNOT_ORG {
	take:
		fa_ch       // channel: [ val(meta), path(assembly_fna) ]
		detected_orgname_ch  // channel: [ val(meta), val(orgname) ]
	main:
		// Update fa_ch with appropriate org_name
		fa_org_ch = fa_ch
			.join(detected_orgname_ch,remainder:true)
			.map({meta,fa,detected_org_name -> [meta,fa,org_name(meta)?:detected_org_name]})

		// MLST typing
		cgemlst_ch = fa_org_ch
			.filter({!params.skip_cgemlst})
			.map({meta,fa,org_name -> [meta, fa, tool_args('cgemlst',meta,org_name)]})
			.filter({meta,fasta,args -> args!=null})
			| CGEMLST
		MLST_ch = fa_org_ch
			.filter({!params.skip_MLST})
			.map({meta,fa,org_name -> [meta, fa, tool_args('MLST',meta,org_name)]})
			.filter({meta,fasta,args -> args!=null})
			| MLST

		// PROKKA annotations
		prokka_ch = fa_org_ch
			.filter({!params.skip_prokka})
			.map({meta,fa,org_name -> [meta, fa, tool_args('prokka',meta,org_name)]})
			.filter({meta,fasta,args -> args!=null})
			| PROKKA_RUN

	emit:
		cgemlst = cgemlst_ch
		MLST = MLST_ch
		prokka = prokka_ch
		org_name = fa_org_ch.map({meta,fa,org_name -> [meta,org_name]})
}


workflow {
	main:
		// Validate parameters and print summary of supplied ones
		validateParameters()
		log.info(paramsSummaryLog(workflow))

		// -------------------
		// Prepare SampleSheet
		// -------------------
		ss = get_samplesheet()

		// ------------------------------------------------------------------			
		// Convert input files formats
		// ------------------------------------------------------------------
		// Uncompress .fasta.gz files when needed
		ss.asm_ch = ss.asm_ch.branch({meta,f -> 
			gz: f.name =~ /\.gz$/
			fa: true
		})
		ss.asm_ch = ss.asm_ch.fa.mix(GZIP_DECOMPRESS_FASTA(ss.asm_ch.gz))

		// Convert BAM to FASTQ.GZ when long_reads are given in BAM format
		ss.lr_ch = ss.lr_ch.branch({meta,f -> 
			bam: f.name =~ /\.(bam|cram)$/
			fq: true
		})
		ss.lr_ch = ss.lr_ch.fq.mix(SAMTOOLS_FASTQ(ss.lr_ch.bam))
		
		// -------------------
		// Reads processing
		// -------------------
		PLASMIDFINDER_LONGREAD(ss.lr_ch.filter({!params.skip_plasmidfinder_longread}))
		RESFINDER_LONGREAD(ss.lr_ch.filter({!params.skip_resfinder_longread}),'nanopore')
		PLASMIDFINDER_SHORTREAD(ss.sr_ch.filter({!params.skip_plasmidfinder_shortread}))
		RESFINDER_SHORTREAD(ss.lr_ch.filter({!params.skip_resfinder_shortread}),'illumina')
			
		// ---------------------------------------------------------------------
		// Tools that can run directly on a FASTA witout specifying an organism
		// ---------------------------------------------------------------------
		fai_ch = SAMTOOLS_FAIDX(ss.asm_ch)

		// CGE - RESFINDER
		resfinder_ch = RESFINDER(ss.asm_ch.filter({!params.skip_resfinder}),'fasta')

		// Plasmid typing
		plasmidfinder_ch = ss.asm_ch.filter({!params.skip_plasmidfinder}) | PLASMIDFINDER
		
		// NCBI AMRfinder+
		if (params.skip_amrfinderplus) {
				amrfinderplus_ch = Channel.empty()
		} else {
				amrfinderplus_db = AMRFINDERPLUS_UPDATE()
				amrfinderplus_ch = AMRFINDERPLUS_RUN(
						ss.asm_ch
							.map({meta,fasta -> [meta,fasta,tool_args('amrfinderplus',meta)]})
							.filter({meta,fasta,args -> args!=null}),
						amrfinderplus_db
				)
		}
		
		// MOBsuite - MOBtyper
		mobtyper_ch = ss.asm_ch
			.filter({!params.skip_mobtyper})
			.map({meta,fasta -> [meta,fasta,tool_args('mobtyper',meta)]})
			.filter({meta,fasta,args -> args!=null})
        | MOBTYPER_RUN


    	// Run orgfinder to auto detect organism
		orgfinder_ch = ss.asm_ch.filter({!params.skip_orgfinder}) | ORGFINDER_DETECT
			
		// Speciator
		speciator_ch = ss.asm_ch.filter({!params.skip_speciator}) | SPECIATOR

			
		// ---------------------------------------------------------------------
		// Organism specific tools
		// ---------------------------------------------------------------------
		ann_ch = ANNOT_ORG(ss.asm_ch,orgfinder_ch.org_name)


		// Results aggregation and reporting
		MULTIREPORT(
			ss.asm_ch,
			fai_ch,
			ann_ch.org_name,
			orgfinder_ch.orgfinder,
			amrfinderplus_ch,
			resfinder_ch,
			mobtyper_ch,
			plasmidfinder_ch,
			ann_ch.cgemlst,
			ann_ch.MLST,
			ann_ch.prokka,
			RESFINDER_LONGREAD.out,
			PLASMIDFINDER_LONGREAD.out,
			RESFINDER_SHORTREAD.out,
			PLASMIDFINDER_SHORTREAD.out
		)

	publish:
		// Input assembly
	    multireport_folder = MULTIREPORT.out.folder
    	html_report        = MULTIREPORT.out.html
    	xlsx_report        = MULTIREPORT.out.xlsx
}


output {	
	multireport_folder {
		path { x -> x[1] >> '.' }
		mode 'copy'
	}
	html_report {
		path { x -> x[1] >> "${x[0]}" }
		mode 'copy'
	}
	xlsx_report {
		path { x -> x[1] >> "${x[0]}" }
		mode 'copy'
	}
}



