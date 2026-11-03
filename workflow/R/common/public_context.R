# Trusted bundled controls, authoritative acquisitions and explicit reconstructed inputs.
public_analysis_context <- function(scope, manifest_path, output_root, root=getwd(), create=TRUE) {
  root <- normalizePath(root,mustWork=TRUE)
  if (!file.exists(file.path(root,'PUBLIC_MANIFEST.tsv'))) stop('Run from the public release root.')
  output_root <- resolve_output_directory(output_root,root)
  inside <- function(p,d) identical(p,d) || startsWith(p,paste0(d,'/'))
  clean_path <- function(p) {
    if(length(p)!=1L || is.na(p) || !nzchar(trimws(p))) stop('Invalid input path: empty identity.')
    p <- path.expand(p)
    parts <- strsplit(p,'/',fixed=TRUE)[[1L]]
    if(any(parts=='..')) stop('Invalid input path: traversal is not an input identity.')
    parts <- parts[nzchar(parts) & parts!='.']
    paste0(if(startsWith(p,'/'))'/' else '',paste(parts,collapse='/'))
  }
  key_of <- function(p) {
    p <- clean_path(p)
    q <- if(startsWith(p,'/'))p else file.path(root,p)
    if(file.exists(q)) q <- normalizePath(q,mustWork=TRUE)
    if(inside(q,root)) substring(q,nchar(root)+2L) else p
  }
  pub_path <- file.path(root,'PUBLIC_MANIFEST.tsv')
  pub_hash <- digest::digest(file=pub_path,algo='sha256')
  pub <- read.delim(pub_path,quote='',colClasses='character')
  if(anyDuplicated(pub$path)) stop('Duplicate public manifest identity.')
  checked <- new.env(parent=emptyenv()); resources <- new.env(parent=emptyenv())
  produced <- new.env(parent=emptyenv()); written <- character()
  trusted_file <- function(key) {
    j <- match(key,pub$path)
    if(is.na(j) || pub$sha256[j]=='SELF_HASH_EXCLUDED') stop('Missing bundled authority: ',key)
    p <- file.path(root,key)
    if(!file.exists(p) || dir.exists(p) || normalizePath(p,mustWork=TRUE)!=p ||
       digest::digest(file=p,algo='sha256')!=pub$sha256[j]) stop('Bundled input integrity mismatch: ',key)
    assign(p,pub$sha256[j],checked);p
  }
  aliases <- read.delim(trusted_file('workflow/config/bundled_input_aliases.tsv'),quote='',colClasses='character')
  catalog <- read.delim(trusted_file('workflow/config/reconstruction_inputs.tsv'),quote='"',colClasses='character',check.names=FALSE)
  acquisitions <- read.delim(trusted_file('workflow/config/external_input_identities.tsv'),quote='',colClasses='character')
  if(anyDuplicated(aliases$logical_input) || any(!aliases$bundled_path %in% pub$path)) stop('Invalid bundled input aliases.')
  bundled_key <- function(key) {
    j <- match(key,aliases$logical_input)
    dest <- if(!is.na(j))aliases$bundled_path[j] else key
    if(dest %in% pub$path)dest else NA_character_
  }
  authoritative_hash <- function(keys) {
    exact <- catalog$binding_rule=='EXACT_ACQUISITION'
    h <- catalog$authority_sha256[exact & catalog$logical_input %in% keys]
    j <- acquisitions$binding_rule=='EXACT_ACQUISITION' & acquisitions$object %in% basename(keys)
    h <- unique(c(h,acquisitions$sha256[j]));h <- h[nzchar(h)]
    if(length(h)>1L) stop('Conflicting authoritative acquisition identities.')
    h
  }
  if (!file.exists(manifest_path) || dir.exists(manifest_path)) stop('Missing input manifest: ',manifest_path,'; see workflow/config/reconstruction_inputs.tsv. No fallback or download.')
  manifest_path <- normalizePath(manifest_path,mustWork=TRUE)
  m <- read.delim(manifest_path,quote='',comment.char='',colClasses='character',check.names=FALSE)
  if (!identical(names(m),c('logical_input','local_path','sha256')) || anyNA(m) || anyDuplicated(m$logical_input) || any(!nzchar(m$logical_input))) stop('Input manifest must have unique logical_input, local_path, sha256 columns.')
  if (nrow(m) && any(!grepl('^[a-f0-9]{64}$',m$sha256))) stop('Input manifest needs exact SHA-256 values.')
  m$key <- vapply(m$logical_input,key_of,'')
  if(anyDuplicated(m$key)) stop('Duplicate normalized input identity.')
  for(i in seq_len(nrow(m))) {
    if(!is.na(bundled_key(m$key[i]))) stop('Bundled input cannot be shadowed by a reader mapping: ',m$logical_input[i])
    p <- path.expand(m$local_path[i]);if(!startsWith(p,'/'))p <- file.path(dirname(manifest_path),p)
    if(!file.exists(p) || dir.exists(p)) stop('Missing declared input: ',m$logical_input[i],'; supply the file explicitly; no fallback or download.')
    m$local_path[i] <- normalizePath(p,mustWork=TRUE)
    if(inside(m$local_path[i],output_root)) stop('Input cannot be inside the output directory.')
    if(digest::digest(file=p,algo='sha256')!=m$sha256[i]) stop('Declared input integrity mismatch: ',m$logical_input[i])
  }
  manifest_sha256 <- digest::digest(file=manifest_path,algo='sha256')
  authority_guard <- function() {
    if(digest::digest(file=pub_path,algo='sha256')!=pub_hash) stop('Public manifest changed during execution.')
    for(key in c('workflow/config/bundled_input_aliases.tsv','workflow/config/reconstruction_inputs.tsv','workflow/config/external_input_identities.tsv')) trusted_file(key)
    if(digest::digest(file=manifest_path,algo='sha256')!=manifest_sha256) stop('Input manifest changed during execution.')
    invisible(TRUE)
  }
  output_path <- function(p) {
    p <- clean_path(p)
    if(!startsWith(p,'/')) p <- file.path(output_root,p)
    if(!startsWith(p,paste0(output_root,'/'))) stop('Invalid output path: outside declared directory.')
    ancestor <- p
    repeat {
      link <- Sys.readlink(ancestor)
      if(!is.na(link) && nzchar(link)) stop('Invalid output path: symbolic link.')
      if(ancestor==output_root)break
      ancestor <- dirname(ancestor)
    }
    if(dir.exists(output_root) && normalizePath(output_root,mustWork=TRUE)!=output_root) stop('Invalid output path: output root changed.')
    p
  }
  input_one <- function(p) {
    authority_guard()
    p <- clean_path(p)
    absolute <- if(startsWith(p,'/'))p else file.path(root,p)
    canonical <- if(file.exists(absolute))normalizePath(absolute,mustWork=TRUE) else absolute
    if(inside(absolute,output_root) || inside(canonical,output_root)) {
      q <- output_path(absolute)
      if(!file.exists(q) || dir.exists(q) || !inside(canonical,output_root) || q!=canonical ||
         !exists(canonical,produced,inherits=FALSE)) stop('Unregistered or escaped run-produced input: ',p)
      if(digest::digest(file=canonical,algo='sha256')!=get(canonical,produced,inherits=FALSE)) stop('Run-produced input integrity mismatch: ',p)
      return(canonical)
    }
    if(identical(canonical,manifest_path)) return(manifest_path)
    key <- key_of(p);dest <- bundled_key(key)
    if(!is.na(dest)) return(trusted_file(dest))
    i <- match(key,m$key)
    if(is.na(i)) i <- match(canonical,m$local_path)
    if(is.na(i)) stop('Missing declared input: ',key,'; see REPRODUCIBILITY.md and reconstruction_inputs.tsv. No private fallback or download.')
    q <- m$local_path[i]
    if(!file.exists(q) || dir.exists(q) || normalizePath(q,mustWork=TRUE)!=q) stop('Declared input target changed: ',key)
    actual <- digest::digest(file=q,algo='sha256')
    if(actual!=m$sha256[i]) stop('Declared input integrity mismatch: ',key)
    h <- authoritative_hash(c(key,m$key[m$local_path==q]))
    if(length(h) && actual!=h) stop('Authoritative acquisition SHA-256 mismatch: ',key)
    assign(q,actual,checked);q
  }
  input <- function(p) vapply(p,input_one,'',USE.NAMES=FALSE)
  output <- function(role,name='') file.path(output_root,role,name)
  guard <- function(...) {
    authority_guard()
    for(p in ls(checked,all.names=TRUE)) if(digest::digest(file=p,algo='sha256')!=get(p,checked)) stop('Input changed during execution: ',p)
    invisible(TRUE)
  }
  allowed <- function(p) {
    p <- output_path(p);parent <- dirname(p)
    if(file.exists(p) && !p %in% written) stop('Invalid output path: refusing existing unowned content.')
    if(!dir.exists(parent) && !dir.create(parent,recursive=TRUE,showWarnings=FALSE)) stop('Output path creation failed: ',parent)
    if(!dir.exists(parent)) stop('Output path parent missing.')
    written <<- union(written,p)
    p
  }
  is_link <- function(p) {link<-Sys.readlink(p);!is.na(link)&&nzchar(link)}
  atomic <- function(p, writer, validate) {
    p <- allowed(p);stage <- paste0(p,'.stage')
    if(file.exists(stage) || is_link(stage)) stop('Unresolved staged output: ',stage)
    writer(stage);validate(stage)
    if(!file.rename(stage,p)) stop('Atomic promotion failed: ',p)
    p <- output_path(p)
    if(normalizePath(p,mustWork=TRUE)!=p) stop('Invalid output path: changed after write.')
    assign(p,digest::digest(file=p,algo='sha256'),produced)
    p
  }
  read_tsv <- function(p, ...) utils::read.delim(input(p), stringsAsFactors = FALSE,
    check.names = FALSE, quote = '', comment.char = '', ...)
  read_rds <- function(p, ...) base::readRDS(input(p), ...)
  write_tsv <- function(x, p, na.strings = 'NA') atomic(p, function(f) {
    old <- options(digits = 17); on.exit(options(old))
    utils::write.table(x, f, sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA')
  }, function(f) {
    classes <- vapply(x, function(z) if (is.numeric(z)) 'numeric' else
      if (is.logical(z)) 'logical' else 'character', '')
    z <- utils::read.delim(f, check.names = FALSE, quote = '', comment.char = '',
      stringsAsFactors = FALSE, colClasses = classes, na.strings = na.strings)
    if (!identical(names(z), names(x)) || nrow(z) != nrow(x)) stop('TSV schema readback')
    for (k in names(x)) {
      if (!identical(is.na(x[[k]]), is.na(z[[k]]))) stop('TSV missingness readback: ', k)
      ok <- if (is.numeric(x[[k]])) isTRUE(all.equal(as.numeric(x[[k]]), z[[k]], tolerance = 1e-13))
        else identical(as.character(x[[k]]), as.character(z[[k]]))
      if (!ok) stop('TSV value readback: ', k)
    }
  })
  write_rds <- function(x, p) atomic(p, function(f) saveRDS(x, f, version = 3, compress = FALSE),
    function(f) if (!isTRUE(all.equal(x, base::readRDS(f), tolerance = 0,
      check.environment = FALSE))) stop('RDS value readback'))
  write_json <- function(x, p) atomic(p, function(f) jsonlite::write_json(x, f,
    auto_unbox = TRUE, pretty = TRUE, digits = NA, na = 'null', null = 'null'),
    function(f) jsonlite::read_json(f))
  write_text <- function(x, p) atomic(p, function(f) writeLines(unname(x), f),
    function(f) if (!identical(readLines(f), unname(x))) stop('Text readback'))
  resource_check <- function(minimum = 8, stage = 'atomic_target') {
    z <- readLines('/proc/meminfo')
    kb <- function(k) as.numeric(sub('^[^:]+:[[:space:]]*([0-9]+).*$', '\\1',
      z[startsWith(z, paste0(k, ':'))])) / 1024^2
    a <- kb('MemAvailable'); s <- kb('SwapTotal') - kb('SwapFree')
    rising <- exists('swap', resources) && s > resources$swap
    resources$swap <- s
    if (a < minimum || (s > 2 && rising)) stop('Resource headroom limit: ', stage)
    invisible(list(available_GiB = a, swap_GiB = s, rising = rising))
  }
  verify_inputs <- function() {
    guard(); p <- sort(ls(checked,all.names=TRUE),method='radix')
    data.frame(path=p,sha256=vapply(p,function(x)get(x,checked),''),stringsAsFactors=FALSE)
  }
  bindings <- function() {
    keys <- unique(c(m$logical_input,aliases$logical_input))
    q <- input(keys)
    data.frame(path=keys,bytes=as.numeric(file.info(q)$size),sha256=vapply(q,function(p)digest::digest(file=p,algo='sha256'),''))
  }
  runtime <- function(profile) {
    if(getRversion()!='4.6.1') stop('Qualified R 4.6.1 is required.')
    z <- read.delim(input(paste0('workflow/environment/',profile,'_packages.tsv')),colClasses='character')
    # Complete qualified package versions, independent of installation prefix.
    for(i in seq_len(nrow(z))) {
      actual <- tryCatch(as.character(packageVersion(z$package[i])),error=function(e)'MISSING')
      if(actual!=as.character(package_version(z$version[i]))) stop('Qualified package mismatch: ',z$package[i],' expected ',z$version[i],'; found ',actual)
    }
    list(R=as.character(getRversion()),profile=profile,packages=z)
  }
  config <- yaml::read_yaml(input('config/analysis_config.yml'))
  config$public_terminal_states <- jsonlite::read_json(input('workflow/config/reconstruction_terminal_states.json'))
  if(create) {
    if(!dir.exists(output_root) && !dir.create(output_root,recursive=TRUE,showWarnings=FALSE)) stop('Output path creation failed.')
    if(!dir.exists(output_root)) stop('Output path was not created.')
  }
  list(scope=scope,root=root,output_root=output_root,input=input,canonical=input,output=output,
    allowed=allowed,atomic=atomic,read_tsv=read_tsv,read_rds=read_rds,write_tsv=write_tsv,
    write_rds=write_rds,write_json=write_json,write_text=write_text,guard=guard,
    verify_inputs=verify_inputs,resource_check=resource_check,bindings=bindings,config=config,
    manifest_path=manifest_path,manifest_sha256=manifest_sha256,runtime=runtime)
}
