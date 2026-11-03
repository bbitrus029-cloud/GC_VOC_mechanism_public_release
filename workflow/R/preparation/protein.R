# Original unscored protein/bridge preparation. No data or identities are embedded.
# All source files and gate comparison ledgers must be independently supplied.
protein_preparation <- function(context) {
  reader <- gcvoc_pdc000614(context)
  list2env(reader,envir=environment())
  s04_root <- context$root
  s04_path <- context$input; s04_read <- context$read_tsv
  s04_out <- function(p) context$output('processed',p)
  s04_manifest <- function(p) context$output('metadata',p)
  s04_qc <- function(p) context$output('qc',p)
  s04_write <- context$write_tsv; s04_rds <- context$write_rds; s04_json <- context$write_json
s04_prepare <- function(settings) {
  context$guard()
  cfg<-context$config
  source<-settings$source_directory
  metadata<-pdc_gate_read_linkedomics_identity(context$input(file.path(source,'GC_metatable.tsi')))
  clinical<-pdc_gate_read_clinical_identity(context$input(file.path(source,'GC_Jan2025_clinical_data.tsi')))
  s01_require(all(metadata$canonical_case_identifier%in%clinical),'clinical case corroboration')
  rna<-pdc_gate_read_masked_matrix(context$input(file.path(source,settings$rna_file)),metadata,'RNA')
  protein<-pdc_gate_read_masked_matrix(context$input(file.path(source,settings$protein_file)),metadata,'protein')
  s01_require(rna$sha256==settings$rna_sha256&&protein$sha256==settings$protein_sha256,'immutable matrix changed')
  s01_require(identical(dim(protein$values),c(11294L,189L))&&identical(dim(rna$values),c(62696L,172L)),'physical dimensions')
  s01_require(sum(!is.finite(rna$values))==0&&sum(!is.finite(protein$values))==254539,'finite masks drift')
  reconciliation<-pdc_gate_reconcile_official_exclusions(s04_root)
  provenance<-pdc_gate_assess_authoritative_protein_provenance(s04_root,protein,reconciliation)
  s01_require(!any(reconciliation$exclusions$authoritative_sample_identifier%in%protein$identity$original_matrix_column),'excluded column reappeared')
  protein<-pdc_gate_attach_exclusions(protein,reconciliation$exclusions$authoritative_sample_identifier,unique(reconciliation$exclusions$canonical_case_identifier))
  rna<-pdc_gate_attach_exclusions(rna,reconciliation$exclusions$authoritative_sample_identifier,unique(reconciliation$exclusions$canonical_case_identifier))
  gate_ledger<-s04_read(s04_path(file.path(settings$protein_gate,'pdc000614_tmt13_gate_pair_and_quadruplet_ledger.tsv')))
  eq<-s04_read(s04_path(file.path(settings$protein_gate,'pdc000614_tmt13_gate_endpoint_equivalence_assessment.tsv')))
  ordinary<-list();bridge<-list();backgrounds<-list();coverage<-list();ledgers<-list()
  add_background<-function(endpoint,layer,genes,prior=NULL) {
    payload<-pdc_gate_background_payload(genes)
    if(!is.null(prior)) {
      s01_require(identical(genes,jsonlite::fromJSON(prior$current_background_features_json))&&
        payload$sha256==prior$current_background_sha256,'exact ordered Gate background mismatch')
    }
    backgrounds[[length(backgrounds)+1L]]<<-data.frame(endpoint_id=endpoint,layer=layer,position=seq_along(genes),gene_id=genes)
  }
  for(axis in c('G','H')) {
    members<-cfg$programs[[paste0(axis,'_core')]]$protein$members;id<-paste0('protein_',axis,'_program')
    z<-pdc_gate_assess_pair_endpoint(protein,members,id)
    expected<-gate_ledger[gate_ledger$endpoint_id==id,];expected<-expected[order(expected$patient_id,method='radix'),];rownames(expected)<-NULL
    s01_require(identical(z$ledger,expected),'ordinary exact ledger mismatch')
    s01_require(z$complete_n==cfg$coverage_contracts$protein_core_common_pairs[[axis]],'axis common pair count')
    add_background(id,'protein',z$background,eq[eq$endpoint_id==id,])
    ordinary[[id]]<-z;ledgers[[id]]<-z$ledger
    coverage[[id]]<-data.frame(endpoint_id=id,members=paste(members,collapse=';'),missing_members='',complete_n=z$complete_n,
      structural_eligibility='ELIGIBLE',reason='EXACT_GATE_COMMON_COMPLETE_PAIRS',model_validity='DEFERRED_UNTIL_AUTHORIZED_INFERENCE')
    id<-paste0('bridge_',axis);b<-pdc_gate_assess_bridge_endpoint(rna,protein,members,id)
    expected<-gate_ledger[gate_ledger$endpoint_id==id,];expected<-expected[order(expected$patient_id,method='radix'),];rownames(expected)<-NULL
    s01_require(identical(b$ledger,expected),'bridge exact ledger mismatch')
    s01_require(b$complete_n==cfg$sources$CPTAC_bridge[[paste0(axis,'_complete_quadruplets')]],'bridge n')
    add_background(id,'protein',b$protein_background,eq[eq$endpoint_id==id,])
    prior<-eq[eq$endpoint_id==id,]
    s01_require(grepl(paste0(';rna_background_feature_n_',length(b$rna_background),
      ';rna_background_sha256_',pdc_gate_background_payload(b$rna_background)$sha256),prior$reason,fixed=TRUE),
      'bridge RNA ordered background hash mismatch')
    add_background(id,'RNA',b$rna_background)
    bridge[[id]]<-b;ledgers[[id]]<-b$ledger
    coverage[[id]]<-data.frame(endpoint_id=id,members=paste(members,collapse=';'),missing_members='',complete_n=b$complete_n,
      structural_eligibility='ELIGIBLE',reason='EXACT_GATE_COMPLETE_QUADRUPLETS',model_validity='DEFERRED_UNTIL_AUTHORIZED_INFERENCE')
  }
  # Materialize the approved deterministic context rule now, without claiming earlier file existence.
  context_time<-format(Sys.time(),tz='UTC',usetz=TRUE)
  for(id in names(cfg$coverage_contracts$context_programs)) {
    members<-cfg$coverage_contracts$context_programs[[id]];missing<-setdiff(members,rownames(protein$values))
    if(length(missing)) {
      z<-list(endpoint_id=id,members=members,complete_n=0L,background=character(),retained_cases=character(),retained_opaque=character())
      state<-'NA_NONESTIMABLE';reason<-'MISSING_FIXED_MEMBER_NO_SUBSTITUTION'
    } else {
      z<-pdc_gate_assess_pair_endpoint(protein,members,id,minimum_pairs=0L)
      state<-if(z$complete_n>=3L)'ELIGIBLE' else 'NA_NONESTIMABLE';reason<-if(state=='ELIGIBLE')'APPROVED_PROGRAM_WIDE_FINITE_MASK_RULE' else 'FEWER_THAN_THREE_COMPLETE_PAIRS'
      ledgers[[id]]<-z$ledger;add_background(id,'protein',z$background)
    }
    ordinary[[id]]<-z
    coverage[[id]]<-data.frame(endpoint_id=id,members=paste(members,collapse=';'),missing_members=paste(missing,collapse=';'),complete_n=z$complete_n,
      structural_eligibility=state,reason=reason,model_validity='DEFERRED_UNTIL_AUTHORIZED_INFERENCE')
  }
  cov<-do.call(rbind,coverage);rownames(cov)<-NULL;bg<-do.call(rbind,backgrounds);rownames(bg)<-NULL
  ledger<-do.call(rbind,ledgers);rownames(ledger)<-NULL
  identities<-rbind(data.frame(layer='protein',protein$identity),data.frame(layer='RNA',rna$identity))
  feature_map<-rbind(data.frame(layer='protein',source_row=seq_len(nrow(protein$values)),gene_id=rownames(protein$values)),
    data.frame(layer='RNA',source_row=seq_len(nrow(rna$values)),gene_id=rownames(rna$values)))
  # Persist one untransformed matrix per layer with opaque columns; no program values are created.
  out<-c(s04_rds(list(protein=protein,rna=rna,ordinary=ordinary,bridge=bridge),s04_out('unscored_inputs.rds')),
    s04_write(cov,s04_manifest('endpoint_coverage.tsv')),s04_write(bg,s04_manifest('ordered_backgrounds.tsv')),
    s04_write(ledger,s04_manifest('pair_quadruplet_membership.tsv')),s04_write(identities,s04_manifest('sample_identity.tsv')),
    s04_write(feature_map,s04_manifest('feature_identity.tsv')),
    s04_json(provenance,s04_qc('publication_provenance_verification.json')),
    s04_write(reconciliation$ledger,s04_manifest('publication_exclusions.tsv')))
  c(out,s04_json(list(context_materialized_utc=context_time,rule='Approved Sections 5.5 and 16.2; no scientific choice added',
    ordinary_G_n=ordinary$protein_G_program$complete_n,ordinary_H_n=ordinary$protein_H_program$complete_n,
    bridge_G_n=bridge$bridge_G$complete_n,bridge_H_n=bridge$bridge_H$complete_n,
    disease_effect_access_occurred=FALSE),s04_qc('ingress_summary.json')))
}

  list(prepare=s04_prepare)
}
