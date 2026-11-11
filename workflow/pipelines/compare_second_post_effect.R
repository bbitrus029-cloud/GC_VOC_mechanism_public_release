# Comparison is intentionally separate from calculation and opens references last.
args<-commandArgs(TRUE)
if(length(args)!=3L || !args[1]%in%c('B','D','E'))stop('Usage: compare_second_post_effect.R B|D|E COMPLETED_COMPUTATION_DIRECTORY REPORT_DIRECTORY')
source('workflow/R/common/public_paths.R')
out<-resolve_output_directory(args[3],normalizePath(getwd()))
run<-normalizePath(args[2],mustWork=TRUE);stage<-args[1]
# Fixed writer contract: independent of discovery in the supplied directory.
roster<-list(
 B=c('results/B_diagnostic.tsv','results/B_donor_values.tsv','diagnostics/B_omissions.tsv','diagnostics/B_bootstrap_draws.tsv','diagnostics/B_g_av.json'),
 D=c('results/D_results.tsv','results/D_patient_deltas.tsv','diagnostics/D_omissions.tsv','diagnostics/D_component_correlations.tsv','diagnostics/D_model_diagnostics.json','diagnostics/D_frozen_reference_moments.tsv','diagnostics/D_influence_summary.tsv'),
 E=c('results/E_results.tsv','results/E_paired_values.tsv','diagnostics/E_bootstrap_draws.tsv','diagnostics/E_omissions.tsv','diagnostics/E_influence_summary.tsv','diagnostics/E_join_population.tsv','diagnostics/E_ratio_identity.tsv','diagnostics/E_aliquot_RNA_bridge.tsv'))
required<-sort(roster[[stage]],method='radix')
qc_required<-c(paste0('qc/',stage,'_run.json'),'qc/verified_input_identities.tsv')
for(d in c('results','diagnostics','qc'))stopifnot(dir.exists(file.path(run,d)),!nzchar(Sys.readlink(file.path(run,d))))
found<-sort(unlist(lapply(c('results','diagnostics'),function(d)paste0(d,'/',list.files(file.path(run,d),recursive=TRUE,all.files=TRUE,no..=TRUE)))),method='radix')
if(!identical(found,required))stop('Incomplete or unexpected scientific output roster for ',stage,': missing=',paste(setdiff(required,found),collapse=','),'; unexpected=',paste(setdiff(found,required),collapse=','))
qc_found<-sort(paste0('qc/',list.files(file.path(run,'qc'),recursive=TRUE,all.files=TRUE,no..=TRUE)),method='radix')
stopifnot(identical(qc_found,sort(qc_required,method='radix')))
for(f in c(required,qc_required)) {
 p<-file.path(run,f);stopifnot(file.exists(p),!dir.exists(p),!nzchar(Sys.readlink(p)),file.info(p)$size>0)
 con<-file(p,'rb');seek(con,where=file.info(p)$size-1,origin='start');last<-readBin(con,'raw',n=1L);close(con)
 stopifnot(identical(last,as.raw(10L)))
}
status<-jsonlite::read_json(file.path(run,'qc',paste0(stage,'_run.json')))
start<-c(B='normalized_member_values',D='patient_component_and_H_deltas',E='paired_member_values')
stopifnot(identical(status$stage,stage),identical(status$status,'COMPUTATION_COMPLETED_NOT_YET_COMPARED'),
 identical(status$new_inference,FALSE),identical(status$starting_point,unname(start[stage])),
 identical(status$runtime$R,'4.6.1'),identical(status$runtime$profile,if(stage=='B')'epithelial' else 'supplementary'),
 length(status$runtime$packages)>0L,length(status$session)>0L)
read<-function(p) {
 if(grepl('\\.json$',p)) {z<-jsonlite::read_json(p,simplifyVector=TRUE);stopifnot(is.list(z),length(z)>0L);return(z)}
 lines<-readLines(p,warn=TRUE);stopifnot(length(lines)>1L,all(nzchar(lines)))
 widths<-count.fields(p,sep='\t',quote='',comment.char='',blank.lines.skip=FALSE)
 stopifnot(length(unique(widths))==1L,widths[1]>1L)
 z<-read.delim(p,check.names=FALSE,quote='',comment.char='',fill=FALSE,stringsAsFactors=FALSE)
 stopifnot(nrow(z)>0L,!anyDuplicated(names(z)));z
}
inputs<-read(file.path(run,'qc','verified_input_identities.tsv'))
stopifnot(identical(names(inputs),c('path','sha256')),!anyNA(inputs),!anyDuplicated(inputs$path),all(grepl('^[a-f0-9]{64}$',inputs$sha256)),
 !any(grepl('/(B_diagnostic|D_results|E_results|endpoint_ledger)\\.tsv$',inputs$path)))
files<-file.path(run,required)
# Read the complete persisted calculation set before opening accepted references.
calculated<-setNames(lapply(files,read),files)
manifest<-read.delim('PUBLIC_MANIFEST.tsv',quote='',colClasses='character');checked<-list()
# Scientific input identity is portable across relocation. Historical catalog
# hashes remain QC provenance, not a requirement to restore old controls.
input_roster<-list(B=c('results/supplementary/second_post_effect/B_member_values.tsv','results/diagnostics/supplementary/second_post_effect/B_reference.tsv','results/diagnostics/supplementary/second_post_effect/B_normalization.tsv'),D=c('results/supplementary/second_post_effect/D_patient_deltas.tsv','results/diagnostics/supplementary/second_post_effect/D_frozen_reference_moments.tsv'),E=c('results/supplementary/second_post_effect/E_paired_values.tsv','results/diagnostics/supplementary/second_post_effect/E_join_population.tsv','results/diagnostics/supplementary/second_post_effect/E_ratio_identity.tsv','results/diagnostics/supplementary/second_post_effect/E_aliquot_RNA_bridge.tsv'))
for(key in c(input_roster[[stage]],'workflow/config/second_post_effect_endpoints.tsv','workflow/config/implementation/second_post_effect_seed_map.tsv')) {
 i<-which(endsWith(inputs$path,paste0('/',key))|inputs$path==key);j<-match(key,manifest$path)
 stopifnot(length(i)==1L,!is.na(j),inputs$sha256[i]==manifest$sha256[j])
}
for(p in files) {
 ref<-file.path(if(basename(dirname(p))=='results')'results/supplementary/second_post_effect' else 'results/diagnostics/supplementary/second_post_effect',basename(p))
 j<-match(ref,manifest$path);stopifnot(!is.na(j),digest::digest(file=ref,algo='sha256')==manifest$sha256[j])
 a<-calculated[[p]];b<-read(ref)
 stopifnot(identical(names(a),names(b)))
 if(is.data.frame(a)) {
  stopifnot(nrow(a)==nrow(b))
  for(k in names(a)) {
   stopifnot(identical(is.na(a[[k]]),is.na(b[[k]])))
   if(is.numeric(a[[k]])&&is.numeric(b[[k]])) {
    stopifnot(isTRUE(all.equal(unname(a[[k]]),unname(b[[k]]),tolerance=1e-12,check.attributes=FALSE)))
    if(k%in%c('p_value','q_value')) {
     good<-is.finite(b[[k]]);stopifnot(all((a[[k]][good]>0)==(b[[k]][good]>0)))
     positive<-good & b[[k]]>0;stopifnot(all(abs(a[[k]][positive]/b[[k]][positive]-1)<1e-12))
    }
   } else stopifnot(identical(a[[k]],b[[k]]))
  }
 } else stopifnot(isTRUE(all.equal(a,b,tolerance=1e-12,check.attributes=FALSE)))
 checked[[basename(p)]]<-list(relative_path=required[match(p,files)],rows=if(is.data.frame(a))nrow(a) else NA_integer_,reference_sha256=manifest$sha256[j],calculated_sha256=digest::digest(file=p,algo='sha256'))
}
checked_paths<-unname(vapply(checked,function(x)x$relative_path,''))
stopifnot(identical(checked_paths,required),length(checked)==length(required))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
jsonlite::write_json(list(status='PASS',stage=stage,expected_outputs=required,checked_outputs=checked_paths,expected_count=length(required),checked_count=length(checked),required_QC=qc_required,completion_state=status$status,complete_output_coverage=TRUE,tolerance=1e-12,positive_tiny_P_q_relative_tolerance=1e-12,order_identities_missingness='exact',references_opened_after_persisted_calculations=TRUE,files=checked),file.path(out,'comparison.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA)
cat('PASS',stage,'persisted replay comparisons:',length(files),'files\n')
