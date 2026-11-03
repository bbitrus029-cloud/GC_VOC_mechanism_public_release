# Walk components in filesystem order: resolve existing parents before handling ..,
# while keeping nonexistent components absolute without creating them.
resolve_output_directory <- function(path, release_root) {
  fail <- function(message) stop('Invalid output path: ', message, call.=FALSE)
  if (length(path) != 1L || is.na(path) || !nzchar(trimws(path)))
    fail('destination must not be empty or whitespace-only.')
  path <- path.expand(path)
  if (startsWith(path, '~')) fail('home-directory expansion is unresolved.')
  if (!startsWith(path, '/')) path <- file.path(release_root, path)
  parts <- strsplit(path, '/', fixed=TRUE)[[1L]]
  parts <- parts[nzchar(parts) & parts != '.']
  out <- '/'
  for (i in seq_along(parts)) {
    part <- parts[i]
    if (part == '..') { out <- dirname(out); next }
    candidate <- file.path(out, part)
    link <- Sys.readlink(candidate)
    if (!is.na(link) && nzchar(link)) {
      # Check the final link (including a dangling link) before resolving it.
      if (i == length(parts)) fail('destination must not be a symbolic link.')
      candidate <- tryCatch(normalizePath(candidate, mustWork=TRUE),
        error=function(e) fail('a parent symbolic link cannot be resolved.'))
      if (!dir.exists(candidate)) fail('an existing parent is not a directory.')
    } else if (file.exists(candidate)) {
      if (!dir.exists(candidate)) fail('destination or parent is an existing file.')
      candidate <- normalizePath(candidate, mustWork=TRUE)
    }
    out <- candidate
  }
  if (out == release_root || startsWith(out, paste0(release_root, '/')))
    fail('choose a directory outside the immutable release tree.')
  if (dir.exists(out)) {
    if (file.access(out, 5L) != 0L)
      fail('destination directory must be readable and searchable to verify emptiness.')
    entries <- list.files(out, all.files=TRUE, no..=FALSE)
    if (length(entries) != 2L || !setequal(entries, c('.', '..')))
      fail('destination directory must be verifiably empty; existing contents are preserved.')
  }
  out
}
