# Extracted maintained functions. Stable scientific identifiers and numerical order are retained.
# Public adaptation: explicit inputs and outputs; original scientific rules retained.

original_integration_module <- function(context) {
  historical_namespace <- 's09_reproducibility/s08_evidence_integration'
  stopifnot(is.list(context), is.function(context$input), is.function(context$output))
  resolve_path <- function(path) {
    if (startsWith(path, context$output_root)) return(path)
    prefix <- paste0('data/processed/', historical_namespace, '/')
    if (startsWith(path,prefix)) return(context$output('processed',substring(path,nchar(prefix)+1L)))
    prefix <- paste0('results/', historical_namespace, '/')
    if (startsWith(path,prefix)) return(context$output('results',substring(path,nchar(prefix)+1L)))
    prefix <- paste0('reports/qc/', historical_namespace, '/')
    if (startsWith(path,prefix)) return(context$output('qc',substring(path,nchar(prefix)+1L)))
    context$input(path)
  }
  read_rds <- function(path,...) base::readRDS(context$input(resolve_path(path)),...)
  read_yaml <- function(path) yaml::read_yaml(context$input(resolve_path(path)))
  read_lines <- function(path,...) readLines(context$input(resolve_path(path)),...)
  write_rds <- function(x,path) {context$write_rds(x,resolve_path(path));path}
  utc <- function() format(Sys.time(),'%Y-%m-%dT%H:%M:%OS6Z',tz='UTC')
  s08_ns <- historical_namespace
  s08_sessions <- paste0('s09_reproducibility/',c('s05_bulk_effects','s06_single_cell_effects','s07_protein_bridge'))
  s08_result <- function(name) file.path('results',s08_ns,name)
  s08_guard <- function() context$guard()
  s08_read <- function(path) utils::read.delim(context$input(resolve_path(path)),check.names=FALSE,colClasses='character',na.strings=NULL,quote='',comment.char='',stringsAsFactors=FALSE)
  s08_hash <- function(path) s01_sha256(resolve_path(path))
  s08_rds <- write_rds
  route_output_strings <- function(x) {
    if(is.character(x)) {
      for(session in c('s05_bulk_effects','s06_single_cell_effects','s07_protein_bridge','s08_evidence_integration'))
        for(role in c('data/processed','results')) x<-gsub(paste0(role,'/',session,'/'),paste0(role,'/s09_reproducibility/',session,'/'),x,fixed=TRUE)
      return(x)
    }
    if(is.data.frame(x)){for(k in names(x))x[[k]]<-route_output_strings(x[[k]]);return(x)}
    if(is.list(x)){for(k in seq_along(x))x[[k]]<-route_output_strings(x[[k]]);return(x)}
    x
  }
  s08_write <- function(x,path,keys=NULL) {
    x<-route_output_strings(x)
    if(!is.null(keys))x<-s08_sort(x,keys)
    y<-x;for(k in names(y))if(is.numeric(y[[k]]))y[[k]]<-ifelse(is.na(y[[k]]),'NA',sprintf('%.17g',y[[k]]))
    context$write_tsv(y,resolve_path(path),na.strings=NULL)
    actual<-s08_read(path);expected<-as.data.frame(lapply(y,as.character),stringsAsFactors=FALSE);expected[is.na(expected)]<-'NA';rownames(expected)<-NULL
    s08_assert(identical(actual,expected),'Integration persisted strings changed')
    path
  }
s08_num <- function(x) { y <- rep(NA_real_,length(x)); ok <- !is.na(x) & !grepl('^NA|^NO_FORMAL',x); y[ok] <- as.numeric(x[ok]); y }

s08_assert <- function(ok,msg) if(!isTRUE(ok)) stop(msg,call.=FALSE)

s08_sort <- function(x,keys) {x<-x[do.call(order,c(unname(x[keys]),list(method='radix'))),,drop=FALSE];rownames(x)<-NULL;x}

s08_bind <- function(xs) {cols<-unique(unlist(lapply(xs,names))); do.call(rbind,lapply(xs,function(x){for(n in setdiff(cols,names(x)))x[[n]]<-'NA';x[cols]}))}

s08_coverage_check <- function(x,cw) {
 s08_assert(!anyDuplicated(x$hypothesis_id) && !anyDuplicated(x$endpoint_id),'Duplicated hypothesis join')
 s08_assert(setequal(x$hypothesis_id,cw$hypothesis_id),'Omitted or extra registered hypotheses')
 z<-cw[match(x$hypothesis_id,cw$hypothesis_id),]
 for(n in intersect(names(cw),names(x))) s08_assert(identical(x[[n]],z[[n]]),paste('Crosswalk semantic mismatch',n))
 sc<-x[x$multiplicity_family %in% c('SC_G_MEMBERS','SC_H_MEMBERS'),]
 s08_assert(all(sc$standard_error_status=='NA_UNAVAILABLE') && all(sc$confidence_interval_status=='NA_UNAVAILABLE') && all(is.na(s08_num(sc$standard_error))) && all(is.na(s08_num(sc$ci_lower))) && all(is.na(s08_num(sc$ci_upper))),'Fabricated single-cell uncertainty')
 s08_assert(all(sc$estimability=='ESTIMABLE') && all(is.finite(s08_num(sc$statistic))) && all(s08_num(sc$numerator_df)>0) && all(s08_num(sc$denominator_df)>0),'Native SC member test lost')
 TRUE
}

s08_projection_check <- function(x) {
 for(p in unique(x$source_table)) {
  a<-s08_read(p);b<-x[x$source_table==p,];a<-a[match(b$source_row_id,a$endpoint_id),,drop=FALSE]
  s08_assert(!anyNA(a$endpoint_id),'Missing source reference')
  for(n in names(a)) s08_assert(identical(a[[n]],b[[n]]),paste('Altered accepted value',p,n))
 }
 TRUE
}

s08_ref <- function(path,key,value) paste0(path,'#',key,'=',value)

s08_assemble <- function(verified) {
 s08_guard();cw<-s08_read('config/hypothesis_estimand_crosswalk.tsv')
 tabs<-list();for(s in s08_sessions) {
  ns<-c('all_program_estimates','all_member_estimates',if(s==s08_sessions[2])'captured_composition_estimates',if(s==s08_sessions[3])'cptac_bridge_estimates')
  for(n in ns){p<-file.path('results',s,paste0(n,'.tsv'));x<-s08_read(p);x$source_table<-p;x$source_row_id<-x$endpoint_id;x$evidence_identifier<-s08_ref(p,'endpoint_id',x$endpoint_id);tabs[[paste(s,n)]]<-x}
 }
 all<-s08_bind(tabs);s08_coverage_check(all,cw);s08_projection_check(all)
 files<-character()
 for(n in c('all_program_estimates','all_member_estimates')) files<-c(files,s08_write(s08_bind(tabs[endsWith(names(tabs),n)]),s08_result(paste0(n,'.tsv')),'endpoint_id'))
 for(n in c('bulk_reproducibility_states','single_cell_state_estimates','captured_composition_estimates','protein_estimates','cptac_bridge_estimates')) {
  s<-if(n=='bulk_reproducibility_states')s08_sessions[1] else if(n %in% c('single_cell_state_estimates','captured_composition_estimates'))s08_sessions[2] else s08_sessions[3]
  x<-s08_read(file.path('results',s,paste0(n,'.tsv')));files<-c(files,s08_write(x,s08_result(paste0(n,'.tsv')),if(n=='bulk_reproducibility_states')'construct' else 'endpoint_id'))
 }
 files<-c(files,s08_write(all,s08_result('hypothesis_coverage.tsv'),'hypothesis_id'))
 # Diagnostics remain separate and are indexed at their original row identities.
 ds<-list();for(i in seq_along(s08_sessions)) {
  p<-file.path('results',s08_sessions[i],c('robustness_diagnostics.tsv','nominal_program_diagnostics.tsv','nominal_protein_diagnostics.tsv')[i]);x<-s08_read(p)
  id<-if(i==1)paste(x$endpoint_id,x$diagnostic_id,sep='|') else x$analysis_id
  ds[[i]]<-data.frame(session=s08_sessions[i],source_table=p,source_row_id=id,evidence_identifier=s08_ref(p,if(i==1)'endpoint_id|diagnostic_id' else 'analysis_id',id),role='same_data_robustness',category=x$category,endpoint_id=x$endpoint_id,terminal_state=x$diagnostic_state,stringsAsFactors=FALSE)
 }
 bp<-'results/s07_protein_bridge/cptac_bridge_estimates.tsv';br<-s08_read(bp)
 ds[[4]]<-data.frame(session=s08_sessions[3],source_table=bp,source_row_id=paste0(br$endpoint_id,'|BCa_sensitivity'),evidence_identifier=s08_ref(bp,'endpoint_id',br$endpoint_id),role='same_data_robustness',category='prespecified_BCa_interval_sensitivity_NO_FORMAL_FAMILY',endpoint_id=br$endpoint_id,terminal_state=br$bca_state,stringsAsFactors=FALSE)
 files<-c(files,s08_write(do.call(rbind,ds),s08_result('diagnostic_source_inventory.tsv'),c('session','source_row_id')))
 a5<-list(formal=context$read_tsv('results/s05_bulk_effects/all_program_estimates.tsv'));a7<-list(protein=context$read_tsv('results/s07_protein_bridge/protein_estimates.tsv'));bs<-context$read_tsv('results/s05_bulk_effects/bulk_reproducibility_states.tsv')
 states<-lapply(c('G','H'),function(axis){b<-bs[bs$construct==axis,];v<-a5$formal$effect[match(c(b$active_source_1_endpoint_id,b$active_source_2_endpoint_id),a5$formal$endpoint_id)];dir<-if(all(is.finite(v)) && v[1]*v[2]>0)sign(v[1]) else NA_real_;p<-a7$protein[a7$protein$endpoint_id==paste0('PDC000614_protein_',axis),];data.frame(construct=axis,bulk_dimension_a=b$dimension_a_state,bulk_dimension_b=b$dimension_b_state,bulk_consensus_direction=dir,protein_effect=p$effect,protein_q=p$q_value,cross_layer_state=s01_cross_layer_state(b$dimension_a_state,dir,p$estimability=='ESTIMABLE',p$q_value,p$effect))})
 states<-do.call(rbind,states);comp<-s01_compensation_compatible(states$bulk_dimension_a[1],states$bulk_consensus_direction[1],states$bulk_dimension_a[2],states$bulk_consensus_direction[2]);states$compensation_compatible<-comp
 files<-c(files,s08_write(states,s08_result('cross_layer_classification_inputs.tsv'),'construct'))
 save<-list(all=all,states=states,bulk=bs,diagnostics=do.call(rbind,ds),files=files)
 p<-file.path('data/processed',s08_ns,'integration.rds');s08_rds(save,p)
 c(files,p)
}

s08_fmt <- function(x) {n<-s08_num(x);ifelse(is.na(n),as.character(x),formatC(n,digits=4,format='g'))}

s08_estimate_text <- function(x) {
 paste0(x$endpoint_id,': n=',x$biological_n,'; ',x$estimability,'; ',x$estimand,'=',s08_fmt(x$effect),'; 95% interval [',s08_fmt(x$ci_lower),', ',s08_fmt(x$ci_upper),']; P=',s08_fmt(x$p_value),'; q=',s08_fmt(x$q_value),'; ',x$multiplicity_family,' (frozen ',x$frozen_family_size,'); ',x$missing_reason)
}

s08_synthesis <- function(assembled) {
 a<-read_rds(file.path('data/processed',s08_ns,'integration.rds'));x<-a$all;st<-a$states
 matrix<-list();k<-0L
 limits<-c(core_bulk='Independent paired tissue cohorts; NAT is adjacent non-tumor, not healthy mucosa. Rank position is not calibrated across sources.',TCGA='Cross-platform sensitivity only; non-estimability is an evidence gap, not zero.',single_cell_state='Donor-aware within-compartment state; compartment-specific members and backgrounds; EPCAM depletion; not dominant cell source.',captured_composition='Captured EPCAM-depleted release composition, separate from state; not unbiased whole-tissue proportions and cannot explain away bulk effects.',ordinary_protein='Protein abundance only; G has six measured members with ALOX12 absent; fixed common pairs; not enzyme activity.',same_patient_bridge='Exploratory same-cohort RNA/protein coupling; overlaps ordinary protein patients; wide intervals do not establish absence of coupling.',context='Fixed context programs cannot replace or rescue G/H; sources and patients overlap core endpoints.',robustness='Same-data diagnostics cannot modify primary labels or create independent replication.',external_anchor='External literature defines the question, not a measured project outcome; competing oral, microbial, systemic, ambient and material sources remain unresolved.')
 for(axis in c('G','H')) for(layer in names(limits)) {
   y<-switch(layer,core_bulk=x[x$endpoint_id %in% unlist(a$bulk[a$bulk$construct==axis,c('active_source_1_endpoint_id','active_source_2_endpoint_id')]),],TCGA=x[x$endpoint_id==paste0('TCGA_STAD_',axis),],single_cell_state=x[x$contract=='C2' & x$construct==paste0(axis,'_core') & x$member=='NA_NOT_APPLICABLE',],captured_composition=x[x$multiplicity_family %in% c('SC_COMPOSITION_PRIMARY','SC_COMPOSITION_STRICT200','SC_COMPOSITION_CLR','SC_COMPOSITION_QUASIBINOMIAL'),],ordinary_protein=x[x$contract=='C3' & (grepl(paste0('^',axis,'_core'),x$construct) | grepl(paste0('^',axis,'__'),x$family_hypothesis)),],same_patient_bridge=x[x$endpoint_id==paste0('CPTAC_bridge_',axis),],context=x[x$construct %in% if(axis=='G')'G_buffer' else c('H_expanded','H_reductive','H_conjugation'),],NULL)
   if(layer=='robustness') {d<-a$diagnostics;refs<-d$evidence_identifier;txt<-paste0('Shared complete diagnostic inventory (not all rows concern this axis); bulk ',axis,': LOO direction stable=',a$bulk$loo_direction_stable_both_core[a$bulk$construct==axis],'; LOO replication sensitive=',a$bulk$loo_replication_label_sensitive[a$bulk$construct==axis],'; score-method sensitive=',a$bulk$score_method_sensitive[a$bulk$construct==axis],'. All source-specific qualifications retained.')} else if(layer=='external_anchor') {refs<-c('docs/plans/APPROVED_ANALYSIS_PLAN.md#4','docs/plans/APPROVED_ANALYSIS_PLAN.md#Appendix_D','reports/qc/pre_amendment_audit/construct/construct_gene_matrix.tsv#ALL_ROWS');txt<-'C4-C10 straight-chain saturated aldehydes are an external question anchor; construct directness does not measure VOCs.'} else {s08_assert(nrow(y)>0,paste('Empty layer',axis,layer));refs<-y$evidence_identifier;txt<-paste(s08_estimate_text(y),collapse=' || ')}
   k<-k+1L;matrix[[k]]<-data.frame(construct=axis,layer=layer,evidence_identifier=paste(refs,collapse=';;'),cross_layer_state=st$cross_layer_state[st$construct==axis],compensation_compatible=st$compensation_compatible[st$construct==axis],same_cohort_dependence_flag=layer %in% c('same_patient_bridge','ordinary_protein','context','robustness','captured_composition','single_cell_state'),principal_limitation=limits[[layer]],claim_ceiling='Abundance-level molecular configuration and bounded biochemical compatibility; requires direct chemical and functional measurements.',provenance_id='docs/plans/APPROVED_ANALYSIS_PLAN.md#18.1-18.3',comparison_scope='Registered bulk versus ordinary protein only; not agreement of all cellular states, members or assays',evidence_summary=txt,evidence_role=switch(layer,core_bulk='independent_population_first_estimate;independent_population_replication',TCGA='cross_platform_sensitivity',single_cell_state='donor_aware_cellular_context',captured_composition='donor_aware_cellular_context',ordinary_protein='paired_human_protein_orthogonal_evidence',same_patient_bridge='prespecified_exploratory_same_cohort_rna_protein_bridge',context='source-specific fixed context; no additional population',robustness='same_data_robustness',external_anchor='external_question_anchor'),dependence_group=switch(layer,core_bulk='GSE29272;GSE27342',TCGA='TCGA_STAD',single_cell_state='GSE206785',captured_composition='GSE206785',ordinary_protein='CPTAC_GC',same_patient_bridge='CPTAC_GC',context='Shares parent source populations',robustness='Shares parent source populations',external_anchor='External literature and construct audit'),classification_timing='Prespecified rule applied after released effects',synthesis_timing='POST_EFFECT_SYNTHESIS',stringsAsFactors=FALSE)
 }
 m<-do.call(rbind,matrix);m$bulk_dimension_a<-st$bulk_dimension_a[match(m$construct,st$construct)];m$bulk_dimension_b<-st$bulk_dimension_b[match(m$construct,st$construct)];m$bulk_consensus_direction<-st$bulk_consensus_direction[match(m$construct,st$construct)];bridge_fields<-c('bridge_state','bootstrap_ci_lower','bootstrap_ci_upper','bca_state','bca_reason','bca_ci_lower','bca_ci_upper','bridge_q_ci_inconsistent','bridge_nominal_p_ci_inconsistent','biological_n')
 for(n in bridge_fields)m[[paste0('native_',n)]]<-'NA_NOT_APPLICABLE'
 for(i in which(m$layer=='same_patient_bridge')) {
  b<-x[x$endpoint_id==paste0('CPTAC_bridge_',m$construct[i]),]
  for(n in bridge_fields)m[[paste0('native_',n)]][i]<-b[[n]]
  m$evidence_summary[i]<-paste0(m$evidence_summary[i],'; bridge_state=',b$bridge_state,'; primary percentile CI [',s08_fmt(b$bootstrap_ci_lower),', ',s08_fmt(b$bootstrap_ci_upper),']; BCa ',b$bca_state,' [',s08_fmt(b$bca_ci_lower),', ',s08_fmt(b$bca_ci_upper),']; q/CI inconsistent=',b$bridge_q_ci_inconsistent,'; nominal-P/CI inconsistent=',b$bridge_nominal_p_ci_inconsistent)
 }
 out<-s08_write(m,s08_result('cross_layer_evidence_matrix.tsv'),c('construct','layer'))
 # Every canonical estimate has an exact result-linked boundary, including unavailable entries.
 claims<-data.frame(claim_id=paste0('CLAIM_',x$hypothesis_id),evidence_identifier=x$evidence_identifier,permitted_claim=s08_estimate_text(x),prohibited_inference='VOC detection or identity; tumor origin; production or clearance capacity; activity; flux; causality; functional validation; tissue-to-breath transfer; explanation of a sensor or clinical classifier.',uncertainty_statement=paste0(x$claim_ceiling,'; ',x$uncertainty_reason,'; ',x$missing_reason,'; non-significance does not establish absence or equivalence. NAT is adjacent non-tumor tissue.'),provenance_id=paste0(x$provenance_id,';;docs/plans/APPROVED_ANALYSIS_PLAN.md#23'),wording_use='Results or caption: report the stated estimand and source-specific scale; no causal interpretation',timing='POST_EFFECT_SYNTHESIS',stringsAsFactors=FALSE)
 extra<-data.frame(claim_id=paste0('CLAIM_SYNTHESIS_',m$construct,'_',m$layer),evidence_identifier=m$evidence_identifier,permitted_claim=ifelse(m$layer=='ordinary_protein',paste0(m$construct,' registered bulk-ordinary-protein comparison: ',m$cross_layer_state,'.'),m$evidence_summary),prohibited_inference=claims$prohibited_inference[1],uncertainty_statement=m$principal_limitation,provenance_id=m$provenance_id,wording_use='Discussion or caption; retain source-specific uncertainty and dependence',timing='POST_EFFECT_SYNTHESIS')
 out<-c(out,s08_write(rbind(claims,extra),s08_result('claim_boundary_ledger.tsv'),'claim_id'))
 out
}



s08_metadata_check <- function(d) {
 loo<-d[grepl('^leave_one_gene_out',d$category),]
 s08_assert(nrow(loo)==18,'S06 LOO inventory incomplete')
 for(i in seq_len(nrow(loo))) {
   expected<-paste(setdiff(strsplit(loo$parent_members[i],',',fixed=TRUE)[[1]],loo$deleted_member[i]),collapse=',')
   s08_assert(loo$members[i]==expected,'Corrected S06 deleted/parent membership lost')
 }
 sub<-d[grepl('fixed_diagnostic__(oxidation_support|substrate_assembly)',d$category),]
 s08_assert(nrow(sub)==2,'S06 submodule inventory incomplete')
 s08_assert(all(sub$parent_state_id=='GSE206785_myeloid_G') && all(sub$population_reference_state_id=='GSE206785_myeloid_H'),'G submodule misattributed to H')
 ctx<-d[grepl('fixed_diagnostic__(G_buffer|H_)',d$category),]
 s08_assert(nrow(ctx)==16,'S06 context inventory incomplete')
 s08_assert(all(ctx$parent_state_id=='NA_NOT_APPLICABLE') && all(ctx$endpoint_id!=ctx$population_reference_state_id),'Context conflated with population reference')
 TRUE
}
# Deterministic member/context projection from the accepted S08 supplement.
# The member-selection expressions are preserved verbatim; inputs are explicit.
s08_endpoint_context <- function(integration, config, state_roles) {
 x<-integration$all;cfg<-config;sc<-state_roles
 members<-vapply(seq_len(nrow(x)),function(i) {
   z<-x[i,];if(z$member!='NA_NOT_APPLICABLE' && z$member!='NA')return(z$member)
   if(z$contract=='C2' && z$endpoint_id %in% sc$endpoint_id)return(sc$available_members[match(z$endpoint_id,sc$endpoint_id)])
   if(z$members!='NA')return(z$members)
   c<-cfg$programs[[z$construct]]$members;if(is.null(c))c<-cfg$programs$context_programs[[z$construct]]
   if(is.null(c))'NA_NOT_APPLICABLE' else paste(c,collapse=';')
 },character(1))
 info<-data.frame(endpoint_id=x$endpoint_id,hypothesis_id=x$hypothesis_id,construct=x$construct,declared_or_available_members=members,member_semantics=ifelse(x$contract=='C2' & x$endpoint_id %in% sc$endpoint_id,'Frozen detected members; no score for non-estimable state','Exact declared member set; estimability and coverage remain controlling'),reference_tissue=ifelse(x$contract=='C4','Paired NAT within each layer; bridge relates two tumor-minus-NAT deltas','NAT: adjacent non-tumor tissue'),scale=x$estimand,biological_unit=x$biological_unit,biological_n=x$biological_n,evidence_identifier=x$evidence_identifier,membership_provenance=ifelse(x$contract=='C2','data/manifests/s03_single_cell_ingress/state_role_freeze.tsv','config/analysis_config.yml; accepted source result'),stringsAsFactors=FALSE)
 s08_assert(nrow(info)==146L && !anyDuplicated(info$endpoint_id) &&
   !anyDuplicated(info$hypothesis_id),'Complete endpoint measurement context required')
 info
}

# Final S09 reconciliation of the immutable eight-column S08 source catalog.
# Logical evidence identities remain historical; physical locations are separate.
s08_source_coverage <- function() {
 p <- context$output('qc','current_input_bindings.tsv')
 context$write_tsv(context$verify_inputs(),p)
 p
}


  run <- function() {
    context$guard()
    assembled<-s08_assemble(NULL);synthesis<-s08_synthesis(NULL)
    endpoint_context<-s08_endpoint_context(read_rds(file.path('data/processed',s08_ns,'integration.rds')),
      read_yaml('config/analysis_config.yml'),s08_read('data/manifests/s03_single_cell_ingress/state_role_freeze.tsv'))
    metadata<-c(s08_write(endpoint_context,s08_result('endpoint_measurement_context.tsv'),'endpoint_id'),
      s08_source_coverage())
    paths<-c(assembled,synthesis,metadata);tables<-paths[endsWith(paths,'.tsv')]
    s08_metadata_check(s08_read('results/s06_single_cell_effects/nominal_program_diagnostics.tsv'))
    list(scope='original_integration',tables=vapply(tables,resolve_path,''),objects=resolve_path(file.path('data/processed',s08_ns,'integration.rds')),
         validation=list(status='PASS',hypotheses=146L,kind='deterministic_projection_from_explicit_rebuilt_dependencies',rendered=FALSE))
  }
  list(run=run,project=s08_assemble,synthesis=s08_synthesis,coverage_check=s08_coverage_check,endpoint_context=s08_endpoint_context,source_coverage=s08_source_coverage,resolve=resolve_path)
}
