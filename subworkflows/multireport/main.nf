
include { RMD_RENDER        } from '../../modules/rmd/render'
include { ORGANIZE_FILES    } from '../../modules/organize_files'

process MULTITABLE {
	container 'registry.gitlab.unige.ch/amr-genomics/rscript:main'
    memory '8 GB'
    cpus 2
    time '30 min'
    input:
		tuple val(meta),path("db")
		each path("lib_typing.R")
    output:
        tuple val(meta), path('multitable.xlsx'), emit:'xlsx'
    script:
        """
		#!/usr/bin/env Rscript
		source("lib_typing.R")
		db <- db_load("db/samples")
		tbl <- list(
			assemblies = summarise_assembly(db),
			contigs = summarise_contigs(db),
			resistances = summarise_resistances(db),
			plasmidfinder_hits = summarise_plasmidfinder_hits(db)
		)
		tbl\$contigs <- left_join(tbl\$contigs,select(tbl\$assemblies,sample_id,assembly_name,species_name,mlst_type),by=c('sample_id','assembly_name'))
		openxlsx::write.xlsx(tbl,file="multitable.xlsx")
        """
}


workflow MULTIREPORT {
	take:
		fa
		fai
		org_name
		orgfinder
		amrfinderplus
		resfinder
		mobtyper
		plasmidfinder
		cgemlst
		MLST
		prokka
		long_resfinder
		long_plasmidfinder
		short_resfinder
		short_plasmidfinder
		
	main:
		ORGANIZE_FILES(
			Channel.empty().mix(
				fa.map(           {meta,file -> [file,"samples/${meta.sample_id}/assemblies/${meta.assembly_name}/assembly.fasta"]}),
				fai.map(          {meta,file -> [file,"samples/${meta.sample_id}/assemblies/${meta.assembly_name}/assembly.fasta.fai"]}),
				org_name.map(     {meta,val  -> [val,"samples/${meta.sample_id}/assemblies/${meta.assembly_name}/org_name.txt"]}),
				orgfinder.map(    {meta,file -> [file,"samples/${meta.sample_id}/assemblies/${meta.assembly_name}/${file.name}"]}),
				amrfinderplus.map({meta,file -> [file,"samples/${meta.sample_id}/assemblies/${meta.assembly_name}/${file.name}"]}),
				resfinder.map(    {meta,file -> [file,"samples/${meta.sample_id}/assemblies/${meta.assembly_name}/${file.name}"]}),
				plasmidfinder.map({meta,file -> [file,"samples/${meta.sample_id}/assemblies/${meta.assembly_name}/${file.name}"]}),
				cgemlst.map(      {meta,file -> [file,"samples/${meta.sample_id}/assemblies/${meta.assembly_name}/${file.name}"]}),
				mobtyper.map(     {meta,file -> [file,"samples/${meta.sample_id}/assemblies/${meta.assembly_name}/${file.name}"]}),
				MLST.map(         {meta,file -> [file,"samples/${meta.sample_id}/assemblies/${meta.assembly_name}/${file.name}"]}),
				prokka.map(       {meta,file -> [file,"samples/${meta.sample_id}/assemblies/${meta.assembly_name}/${file.name}"]}),
				long_resfinder.map(     {meta,file -> [file,"samples/${meta.sample_id}/long_reads/${file.name}"]}),
				long_plasmidfinder.map( {meta,file -> [file,"samples/${meta.sample_id}/long_reads/${file.name}"]}),
				short_resfinder.map(    {meta,file -> [file,"samples/${meta.sample_id}/short_reads/${file.name}"]}),
				short_plasmidfinder.map({meta,file -> [file,"samples/${meta.sample_id}/short_reads/${file.name}"]})
			)
			.collect({[it]})
			.map({["aggregated_folder",it]})
		)

		RMD_RENDER(
			ORGANIZE_FILES.out.map({m,x -> ["multireport.html",x,"indir='${x}'"]}),
			file("${moduleDir}/assets/multireport.Rmd"),
			file("${moduleDir}/assets/lib_typing.R")
		)

		MULTITABLE(
			ORGANIZE_FILES.out.map({m,x -> ["multitable.xlsx",x]}),
			file("${moduleDir}/assets/lib_typing.R")
		)

	emit:
		folder = ORGANIZE_FILES.out
		html = Channel.empty()
		//html   = RMD_RENDER.out  // channel: [ meta, path(html) ]
		xlsx   = MULTITABLE.out  // channel: [ meta, path(xlsx) ]			
}
