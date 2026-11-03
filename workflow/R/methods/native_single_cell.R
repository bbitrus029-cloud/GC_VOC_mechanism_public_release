# Extracted maintained functions. Stable scientific identifiers and numerical order are retained.
# Numerical methods retained from the study implementation.

s06_extract_native <- function(native,members) {
 n<-nrow(native$table);s01_require(length(native$df.test)%in%c(1L,n)&&length(native$df.total)%in%c(1L,n),'native QL degrees-of-freedom dimensions')
 j<-match(members,rownames(native$table));s01_require(!anyNA(j),'native member extraction missing')
 data.frame(member=members,effect=native$table$logFC[j],statistic=native$table$F[j],numerator_df=rep_len(native$df.test,n)[j],denominator_df=rep_len(native$df.total,n)[j],p_value=native$table$PValue[j])
}
