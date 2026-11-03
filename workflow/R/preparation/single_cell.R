# Raw UMI to donor-level preparations. Source access remains subject to its terms.
# All keys below must be mapped in the explicit input manifest; no acquisition occurs.
single_cell_preparation <- function(context) {
 s03_path<-context$input
 s03_qc<-function(p)context$output('qc',p)
 s03_out<-function(p)context$output('processed',p)
 s03_manifest<-function(p)context$output('metadata',p)
 s03_input<-function(id)context$input(paste0('inputs/GSE206785/',id))
 s03_settings<-function()yaml::read_yaml(context$input('config/implementation/s03_settings.yml'))
 s03_write<-context$write_tsv;s03_json<-context$write_json;s03_rds<-context$write_rds
 s03_firewall<-function(...)context$guard()
 s03_resource<-function(label,heavy=FALSE)context$resource_check(if(heavy)8 else 6,label)
 readRDS<-context$read_rds
s03_reconstruct <- function(dependencies) {
  s03_firewall('reconstruct');s03_resource('raw_metrics_start',TRUE)
  context$guard()
  spec <- gse206785_spec()
  md <- gse206785_reconcile_metadata_sources(read.csv(s03_input('fixed_commit_metadata'),check.names=FALSE),
    read.delim(gzfile(s03_input('geo_processed_metadata_corroboration')),check.names=FALSE))$metadata
  identity <- gse206785_validate_metadata(md,list(rows=111140L,samples=48L,patients=24L,
    tissues=spec$tissues,platform_donors=spec$expected$chemistry_donors))
  soft <- gse206785_parse_soft_samples(s03_input('retained_family_soft'))
  bridge <- gse206785_reconcile_soft_metadata(soft$samples,identity$sample_ledger,48L)
  s01_require(setequal(soft$series_supplementary_files,c('GSE206785_metadata.txt.gz','GSE206785_scgex.txt.gz')),'GEO listing mismatch')
  raw <- s03_input('combined_raw_umi_local_derivative')
  parsed <- gse206785_parse_raw_umi_gzip(raw,25287L,111140L)
  parsed$structure$physical_payload_count <- 9L
  parsed$structure$documented_chain <- 'nine_fixed_commit_byte_fragments_one_logical_gzip_stream'
  parsed$structure$combined_bytes <- as.numeric(file.info(raw)$size)
  parsed$structure$combined_sha256 <- s01_sha256(raw)
  release <- gse206785_reconcile_release(parsed$structure,md,25287L,111140L,9L)
  universe <- gse206785_apply_cell_universe(parsed$metrics,md,expected_counts=spec$expected$fixed_universe)
  obj <- list(metadata=md,parsed=parsed,universe=universe,identity=identity,release=release,soft_bridge=bridge,
    cache_identity=list(input_manifest_sha256=context$manifest_sha256))
  out <- s03_rds(obj,s03_out('reconstructed_identity_qc.rds'))
  s03_resource('raw_metrics_end');out
}

s03_key <- function(md) paste(md$donor,md$tissue,md$compartment,sep='|')

s03_aggregate <- function(reconstruction) {
  s03_firewall('aggregate');s03_resource('raw_aggregation_start',TRUE)
  x <- readRDS(reconstruction);md<-x$metadata;flags<-x$universe$flags
  donors <- sort(unique(md$Patient),method='radix');comp<-gse206785_spec()$compartments
  grid <- expand.grid(tissue=c('Normal','Tumor'),donor=donors,compartment=comp,stringsAsFactors=FALSE)
  grid <- grid[,c('donor','tissue','compartment')];grid$unit_id<-s03_key(grid)
  cell <- data.frame(donor=md$Patient,tissue=md$Tissue,compartment=flags$compartment)
  j <- match(s03_key(cell),grid$unit_id)-1L;j[!flags$seven_compartment_pass]<- -1L
  counts <- gse206785_aggregate_raw_umi(s03_input('combined_raw_umi_local_derivative'),j,grid$unit_id,25287L,111140L)
  s01_require(identical(rownames(counts),x$parsed$structure$gene_names),'aggregated gene order')
  s01_validate_numeric_matrix(counts,'S03 pseudobulk',TRUE)
  s01_require(all(counts<=2^53-1),'pseudobulk exceeds exact integer range')
  grid$cell_count <- tabulate(j[j>=0]+1L,nbins=nrow(grid))
  grid$total_umi <- unname(colSums(counts));grid$eligible<-gse206785_apply_unit_eligibility(grid)
  original <- gse206785_aggregate_units(md,x$parsed$metrics,flags)
  oi <- match(s03_key(original),grid$unit_id)
  s01_require(!anyNA(oi)&&identical(as.numeric(original$total_umi),grid$total_umi[oi])&&
    identical(original$cell_count,grid$cell_count[oi]),'independent cell metrics aggregation mismatch')
  members <- unlist(gse206785_spec()$programs,use.names=FALSE)
  s01_require(identical(unname(as.matrix(original[,members])),unname(t(counts[members,oi]))),'core-count aggregation mismatch')
  grid$sample_id <- x$identity$sample_ledger$sample[match(paste(grid$donor,grid$tissue),paste(x$identity$sample_ledger$donor,x$identity$sample_ledger$tissue))]
  grid$chemistry <- x$identity$sample_ledger$chemistry[match(grid$sample_id,x$identity$sample_ledger$sample)]
  grid$unit_status <- ifelse(grid$cell_count==0,'ZERO_CAPTURE_IN_OBSERVED_SAMPLE',ifelse(grid$eligible,'ELIGIBLE','LOW_CAPTURE_OR_UMI'))
  out <- s03_rds(list(counts=counts,units=grid),s03_out('raw_pseudobulks.rds'))
  s03_resource('raw_aggregation_end');out
}
# Encode the independently persisted fibroblast matrix using the frozen byte contract.

s03_prepare <- function(reconstruction,pseudobulks,checkpoint) {
  s03_firewall('prepare');s03_resource('preparation_start')
  x<-readRDS(reconstruction);b<-readRDS(pseudobulks);spec<-gse206785_spec()
  units<-b$units;members<-unlist(spec$programs,use.names=FALSE)
  units[,members]<-t(b$counts[members,,drop=FALSE])
  # Historical structural constructors are used only as pure result-blind functions.
  detection<-gse206785_detect_members(units)
  detection<-detection[order(match(detection$compartment,spec$compartments),match(detection$program,c('G','H')),match(detection$member,members)),];rownames(detection)<-NULL
  structural<-gse206785_build_eligibility_matrix(gse206785_unit_counts(units),detection,24L)
  donor_ledger<-gse206785_donor_compartment_ledger(units,sort(unique(units$donor),method='radix'))
  composition_input<-gse206785_composition_input(x$metadata,x$universe$flags,unique(units$donor))
  composition<-gse206785_build_composition_eligibility(composition_input)
  reconcile<-list(gse206785_core_member_detection.tsv=detection,
    gse206785_compartment_program_eligibility.tsv=structural,
    gse206785_donor_compartment_eligibility.tsv=donor_ledger,
    gse206785_composition_eligibility.tsv=composition)
  # Independent detection reconciliation below remains; no private historical-table oracle is required.
  cfg<-yaml::read_yaml(s03_path('config/analysis_config.yml'))$single_cell_member_model
  prepared<-setNames(lapply(spec$compartments,function(comp) {
    j<-which(b$units$compartment==comp)
    s03_prepare_compartment(b$counts[,j,drop=FALSE],b$units[j,,drop=FALSE],cfg)
  }),spec$compartments)
  for(comp in spec$compartments) {
    d<-prepared[[comp]]$detection;z<-detection[detection$compartment==comp,]
    i<-match(z$member,d$gene_id)
    s01_require(identical(z$available,d$available[i])&&all(z$nonzero_unit_fraction==d$nonzero_fraction[i])&&
      all(z$nonzero_both_tissues_paired_donors==d$paired_both_nonzero[i]),'independent detection reconciliation')
  }
  out<-c(s03_rds(prepared,s03_out('preparation_objects.rds')),
    s03_write(b$units,s03_manifest('pseudobulk_units.tsv')),
    s03_write(detection,s03_manifest('core_member_detection.tsv')),
    s03_write(donor_ledger,s03_manifest('donor_compartment_eligibility.tsv')),
    s03_write(composition,s03_manifest('composition_donor_sets.tsv')),
    s03_write(composition_input,s03_manifest('composition_denominators.tsv')),
    s03_write(gse206785_cell_universe_table(x$universe$counts),s03_manifest('cell_universe.tsv')),
    s03_write(x$identity$sample_ledger,s03_manifest('sample_identity.tsv')))
  structural$historical_marker_status<-structural$marker_audit_status
  structural$marker_audit_status<-NULL
  out<-c(out,s03_write(structural,s03_manifest('structural_states.tsv')))
  backgrounds<-do.call(rbind,lapply(names(prepared),function(comp) {
    p<-prepared[[comp]];data.frame(compartment=comp,order=seq_along(p$background),gene_id=p$background,
      member_model_retained=p$background %in% p$modeling_background)
  }))
  pair_ids<-do.call(rbind,lapply(names(prepared),function(comp) data.frame(compartment=comp,donor_id=prepared[[comp]]$complete_pair_donors)))
  all_detection<-do.call(rbind,lapply(names(prepared),function(comp)cbind(compartment=comp,prepared[[comp]]$detection)))
  out<-c(out,s03_write(backgrounds,s03_manifest('ordered_state_member_backgrounds.tsv')),
    s03_write(pair_ids,s03_manifest('complete_pair_identities.tsv')),
    s03_write(all_detection,s03_manifest('all_gene_detection.tsv')))
  out
}

s03_pool_audit <- function(counts,units,compartments) {
  s01_assert_exact_order(units$unit_id,colnames(counts),'audit unit')
  keys<-paste(units$donor,units$compartment,sep='|')
  donors<-sort(unique(units$donor),method='radix')
  audit<-expand.grid(donor_id=donors,author_label=compartments,stringsAsFactors=FALSE)
  audit$unit_id<-paste(audit$donor_id,audit$author_label,sep='|')
  pooled<-matrix(0,nrow(counts),nrow(audit),dimnames=list(rownames(counts),audit$unit_id))
  audit$cell_count<-0L
  for(i in seq_len(nrow(audit))) {
    j<-which(keys==audit$unit_id[i]);pooled[,i]<-rowSums(counts[,j,drop=FALSE])
    audit$cell_count[i]<-sum(units$cell_count[j])
  }
  audit$total_umi<-unname(colSums(pooled))
  audit$eligible<-audit$cell_count>=20 & audit$total_umi>=10000
  list(counts=pooled,units=audit)
}

s03_panel_scores <- function(pooled,panels,background) {
  counts<-pooled$counts;units<-pooled$units
  s01_assert_exact_order(background,rownames(counts),'audit background')
  s01_assert_exact_order(units$unit_id,colnames(counts),'pooled audit unit')
  s01_validate_numeric_matrix(counts,'pooled audit counts',TRUE)
  genes<-unlist(panels,use.names=FALSE)
  mapping<-data.frame(panel=rep(names(panels),lengths(panels)),marker=genes,
    matching_rows=vapply(genes,function(g)sum(rownames(counts)==g),integer(1)))
  complete<-all(lengths(panels)==4L)&&all(mapping$matching_rows==1L)&&length(background)>=2L
  units$marker_mapping_complete<-complete
  for(panel in names(panels)) units[[panel]]<-NA_real_
  if(complete) for(j in which(units$eligible)) {
    logcpm<-log2(1+1e6*counts[,j]/units$total_umi[j])
    ranks<-(rank(logcpm,ties.method='average')-1)/(length(background)-1)
    for(panel in names(panels)) units[[panel]][j]<-mean(ranks[match(panels[[panel]],background)])
  }
  list(scores=units,mapping=mapping)
}

s03_marker <- function(pseudobulks,preparation) {
  s03_firewall('marker_audit');s03_resource('marker_audit_start')
  b<-readRDS(pseudobulks);settings<-s03_settings();panels<-settings$audit$panels
  pooled<-s03_pool_audit(b$counts,b$units,names(panels));background<-rownames(pooled$counts)
  # Persist the exact ordered background and settings before first real panel scoring.
  out<-s03_write(data.frame(order=seq_along(background),gene_id=background),s03_manifest('ordered_audit_background.tsv'))
  scored<-s03_panel_scores(pooled,panels,background)
  decisions<-lapply(names(panels),function(label)s01_marker_audit(scored$scores,label))
  names(decisions)<-names(panels)
  summary<-do.call(rbind,lapply(names(decisions),function(label) {
    d<-decisions[[label]]
    reason<-if(d$state=='NA_NONESTIMABLE') {
      if(d$eligible_donor_units<8L)'FEWER_THAN_EIGHT_ELIGIBLE_DONOR_UNITS' else 'MISSING_OR_AMBIGUOUS_MARKER_MAPPING'
    } else if(d$state=='PASS') 'EXPECTED_UNIQUE_SELECTION_AND_STRICT_MEDIAN_PASS' else
      paste(c(if(d$expected_selection_fraction<.6)'EXPECTED_SELECTION_BELOW_60_PERCENT',if(!d$median_rule_pass)'EXPECTED_MEDIAN_NOT_STRICTLY_GREATER'),collapse=';')
    data.frame(author_label=label,state=d$state,eligible_donor_units=d$eligible_donor_units,
      expected_unique_selections=sum(d$unit_adjudication==label),ambiguous_units=sum(d$unit_adjudication=='AMBIGUOUS'),
      expected_selection_fraction=d$expected_selection_fraction,median_rule_pass=d$median_rule_pass,reason=reason)
  }))
  scored$scores$selection<-'INELIGIBLE'
  for(label in names(decisions)) {
    j<-which(scored$scores$author_label==label & scored$scores$eligible)
    if(length(decisions[[label]]$unit_adjudication)) scored$scores$selection[j]<-decisions[[label]]$unit_adjudication else
      scored$scores$selection[j]<-'NA_NONESTIMABLE'
  }
  medians<-do.call(rbind,lapply(names(decisions),function(label) data.frame(author_label=label,panel=names(panels),
    median_unrounded_score=if(is.null(decisions[[label]]$donor_balanced_median))rep(NA_real_,7) else decisions[[label]]$donor_balanced_median)))
  confusion<-as.data.frame(table(factor(scored$scores$author_label,levels=names(panels)),
    factor(scored$scores$selection,levels=c(names(panels),'AMBIGUOUS','INELIGIBLE','NA_NONESTIMABLE'))))
  names(confusion)<-c('author_label','selected_panel','donor_units')
  out<-c(out,s03_rds(list(pooled=pooled,scored=scored,decisions=decisions,summary=summary),s03_out('marker_audit.rds')),
    s03_write(scored$scores,s03_qc('marker_donor_scores.tsv')),
    s03_write(scored$mapping,s03_manifest('marker_mapping.tsv')),
    s03_write(summary,s03_qc('marker_decisions.tsv')),
    s03_write(medians,s03_qc('marker_panel_medians.tsv')),
    s03_write(confusion,s03_qc('marker_confusion.tsv')))
  s03_firewall('marker_audit');out
}

s03_freeze <- function(preparation,marker) {
  s03_firewall('freeze')
  audit<-readRDS(s03_out('marker_audit.rds'));summary<-audit$summary
  states<-setNames(summary$state,summary$author_label)
  propagation<-s01_propagate_marker_states(states);composition<-s01_composition_gate(states)
  crosswalk<-s01_read_tsv(s03_path('config/hypothesis_estimand_crosswalk.tsv'))
  crosswalk<-crosswalk[crosswalk$source=='GSE206785',];rownames(crosswalk)<-NULL
  structural<-s01_read_tsv(s03_manifest('structural_states.tsv'))
  structural$endpoint_id<-paste('GSE206785',structural$compartment,structural$program,sep='_')
  structural$marker_state<-unname(states[structural$compartment])
  structural$endpoint_eligibility<-'NA_NONESTIMABLE'
  structural$reason<-structural$gate_role_reason
  structural$gate_role<-structural$gate_role_reason<-NULL
  for(i in seq_len(nrow(structural))) {
    id<-paste(structural$compartment[i],structural$program[i],sep='_')
    structural$reason[i]<-paste(c(if(!structural$minimum_eight_paired_donors_pass[i])'FEWER_THAN_EIGHT_COMPLETE_PAIRS',
      if(!structural$program_coverage_pass[i])'FROZEN_PROGRAM_COVERAGE_FAILED',
      if(id %in% names(propagation) && propagation[id]!='ELIGIBLE')paste0('LINKED_MARKER_',structural$marker_state[i])),collapse=';')
    if(structural$state_estimability[i]=='structurally_estimable' && id %in% names(propagation) && propagation[id]=='ELIGIBLE') {
      structural$endpoint_eligibility[i]<-'ELIGIBLE';structural$reason[i]<-'STRUCTURE_AND_LOCAL_MARKER_PASS'
    }
  }
  crosswalk$endpoint_eligibility<-'NA_NONESTIMABLE';crosswalk$reason<-'NOT_EVALUATED';crosswalk$structural_n<-NA_integer_
  crosswalk$marker_state<-'NOT_APPLICABLE';crosswalk$estimates_status<-'NOT_RUN_S03';crosswalk$owner_acceptance<-'READER_RECONSTRUCTION_NOT_ORIGINAL_FREEZE'
  p<-readRDS(s03_out('preparation_objects.rds'))
  sets<-s01_read_tsv(s03_manifest('composition_donor_sets.tsv'))
  for(i in seq_len(nrow(crosswalk))) {
    id<-crosswalk$endpoint_id[i];family<-crosswalk$multiplicity_family[i]
    if(id %in% structural$endpoint_id) {
      r<-structural[match(id,structural$endpoint_id),]
      crosswalk$endpoint_eligibility[i]<-r$endpoint_eligibility;crosswalk$reason[i]<-r$reason
      crosswalk$structural_n[i]<-r$complete_paired_donors;crosswalk$marker_state[i]<-r$marker_state
    } else if(family %in% c('SC_G_MEMBERS','SC_H_MEMBERS')) {
      tokens<-strsplit(crosswalk$family_hypothesis[i],'__',fixed=TRUE)[[1]];comp<-tokens[1];gene<-tokens[2]
      z<-p[[comp]];d<-z$detection[match(gene,z$detection$gene_id),]
      okay<-length(z$complete_pair_donors)>=8 && isTRUE(d$available) && states[comp]=='PASS' && gene %in% z$modeling_background
      crosswalk$endpoint_eligibility[i]<-if(okay)'ELIGIBLE' else 'NA_NONESTIMABLE'
      crosswalk$reason[i]<-if(okay)'MEMBER_PREPARATION_AND_LOCAL_MARKER_PASS' else paste(c(
        if(states[comp]!='PASS')paste0('LINKED_MARKER_',states[comp]),
        if(length(z$complete_pair_donors)<8)'FEWER_THAN_EIGHT_COMPLETE_PAIRS',
        if(!isTRUE(d$available))'MEMBER_DETECTION_FAILED',if(!gene %in% z$modeling_background)'MEMBER_ABSENT_MODELING_BACKGROUND'),collapse=';')
      crosswalk$structural_n[i]<-length(z$complete_pair_donors);crosswalk$marker_state[i]<-states[comp]
    } else if(startsWith(family,'SC_COMPOSITION_')) {
      crosswalk$endpoint_eligibility[i]<-composition$state
      crosswalk$reason[i]<-if(composition$estimable)'ALL_SEVEN_MARKERS_PASS' else paste0('ALL_SEVEN_MARKERS_REQUIRED;NONPASS=',paste(names(states)[states!='PASS'],collapse=','))
      crosswalk$structural_n[i]<-switch(family,SC_COMPOSITION_PRIMARY=sum(sets$primary_100),SC_COMPOSITION_CLR=sum(sets$primary_100),SC_COMPOSITION_STRICT200=sum(sets$strict_200),SC_COMPOSITION_QUASIBINOMIAL=sum(sets$count_aware))
      crosswalk$marker_state[i]<-if(composition$estimable)'ALL_PASS' else 'NOT_ALL_PASS'
    } else stop('unrecognized S03 endpoint')
  }
  s01_require(nrow(crosswalk)==56L&&!any(crosswalk$reason=='NOT_EVALUATED'),'complete crosswalk freeze')
  families<-aggregate(endpoint_id~multiplicity_family+frozen_family_size,data=crosswalk,FUN=length)
  member<-crosswalk[crosswalk$multiplicity_family %in% c('SC_G_MEMBERS','SC_H_MEMBERS'),]
  cfg<-yaml::read_yaml(s03_path('config/analysis_config.yml'))$single_cell_member_model
  member$future_standard_error_status<-cfg$standard_error_status
  member$future_confidence_interval_status<-cfg$confidence_interval_status
  member$uncertainty_reason<-cfg$uncertainty_reason
  c(s03_write(structural,s03_manifest('state_role_freeze.tsv')),
    s03_write(crosswalk,s03_manifest('endpoint_role_freeze.tsv')),
    s03_write(member,s03_manifest('member_role_freeze.tsv')),
    s03_write(families,s03_manifest('frozen_family_inventory.tsv')),
    s03_json(composition,s03_manifest('composition_role_freeze.json')))
}

list(s03_reconstruct=s03_reconstruct,s03_aggregate=s03_aggregate,s03_prepare=s03_prepare,s03_pool_audit=s03_pool_audit,s03_panel_scores=s03_panel_scores,s03_marker=s03_marker,s03_freeze=s03_freeze)
}
