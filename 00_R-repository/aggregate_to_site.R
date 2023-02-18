library(metabaR)

aggregate_to_site <- function(dt_clean, stats) {
  # stats="FoD" or "logASV"
  # aggregate motus
  dt_motuAgg <- aggregate_motus(dt_clean, groups=dt_clean$motus$TAXID)
  
  # aggregate pcrs by site, output frequency of detection 
  dt_st <- dt_motuAgg
  dt_st$pcrs <- left_join(dt_motuAgg$pcrs, dt_motuAgg$samples %>% mutate(sample_id=rownames(.)), by="sample_id")
  rownames(dt_st$pcrs) <- rownames(dt_motuAgg$pcrs)
  dt_st$pcrs$sample_id <- dt_st$pcrs$site
  
  # load the new 'samples' table: sites
  sites <- read.table("./02_metabaR/surfsedi_sites.txt", header=T, sep='\t')
  rownames(sites) <- sites$site
  dt_st$samples <- sites
  
  # frequency of detection at sites
  if (stats=="FoD") {
    data <- aggregate_pcrs(dt_st, FUN=FUN_agg_pcrs_prob)
  } else {
    data <- aggregate_pcrs(dt_st, FUN=FUN_agg_pcrs_mean)
  }
  
                 
  
  return(data)
}