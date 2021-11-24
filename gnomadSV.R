gnomadSV <- function(path){
  #path <- "~/Documents/PhD project/SVEnsemble/gnomad_v2_sv.sites.vcf.gz"
  gnomad.vcf <- readVcf(path)
  gnomad.bnd.gr <- rowRanges(gnomad.vcf)
  mcols(gnomad.bnd.gr) <- cbind(mcols(gnomad.bnd.gr), info(gnomad.vcf))
  testthat::expect_equal(sum(end(gnomad.bnd.gr)-start(gnomad.bnd.gr) != 0), 0) #checking all start equals end
  testthat::expect_equal(sum(is.na(gnomad.bnd.gr$END)), 0) #checking no END is na.
  
  #gnomad.bnd.gr <- gnomad.bnd.gr %>% filter(seqnames==CHR2) #removing interchromosomal events
  gnomad.bnd.grs <- GRangesList(inter.chr = gnomad.bnd.gr %>% filter(seqnames!=CHR2),
                                   intra.chr = gnomad.bnd.gr %>% filter(seqnames==CHR2 | is.na(CHR2))) #gnomadSV v2.1 contains many 
  
  #intra_chrs:
  gnomad.bnd.grs$intra.chr <- gnomad.bnd.grs$intra.chr %>% filter(!END-end<0) #removing records where end is greater than start (>600 of these guys)
  gnomad.bnd.grs$intra.chr <- gnomad.bnd.grs$intra.chr %>% mutate(end=start, sourceId=names(.), partner=paste(names(.), 'bp2', sep = "_"))
  gnomad.bnd.grs$intra.chr_bp2 <- gnomad.bnd.grs$intra.chr %>% `end<-`(., value=.$END) %>% `start<-`(., value=.$END) %>% mutate(sourceId=partner, partner=names(.)) %>% `names<-`(., .$sourceId)
  testthat::expect_length(partner(c(gnomad.bnd.grs$intra.chr, gnomad.bnd.grs$intra.chr_bp2)), length(c(gnomad.bnd.grs$intra.chr, gnomad.bnd.grs$intra.chr_bp2))) #testing for unpaired bps
  
  #inter_chrs:
  gnomad.bnd.grs$inter.chr <- gnomad.bnd.grs$inter.chr %>% mutate(end=start, sourceId=names(.), partner=paste(names(.), 'bp2', sep = "_"))
  gnomad.bnd.grs$inter.chr_bp2 <- GRanges(seqnames = gnomad.bnd.grs$inter.chr$CHR2, ranges = IRanges(gnomad.bnd.grs$inter.chr$END, width=1, names=gnomad.bnd.grs$inter.chr$partner)) %>% 
    `mcols<-`(., value=mcols(gnomad.bnd.grs$inter.chr)) %>% mutate(sourceId=names(.), partner=gnomad.bnd.grs$inter.chr$sourceId) 
  testthat::expect_length(partner(c(gnomad.bnd.grs$inter.chr, gnomad.bnd.grs$inter.chr_bp2)), length(c(gnomad.bnd.grs$inter.chr, gnomad.bnd.grs$inter.chr_bp2)))
  
  testthat::expect_length(partner(unlist(gnomad.bnd.grs, use.names = FALSE)), length(unlist(gnomad.bnd.grs)))
  unlist(gnomad.bnd.grs, use.names = FALSE) #use.names=TRUE will append list name as prefix to gr names, thus failing the partner().
}
