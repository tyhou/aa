# Shared functions for haploid ALT-allele frequency and crossover analysis.

read_haploid_vcf <- function(vcf_file) {
  lines <- readLines(vcf_file)
  header_i <- grep("^#CHROM\\t", lines)
  if (length(header_i) != 1) stop("Cannot identify one #CHROM line in ", vcf_file)
  header <- strsplit(lines[header_i], "\t", fixed = TRUE)[[1]]
  samples <- header[10:length(header)]
  records <- lines[(header_i + 1):length(lines)]
  records <- records[nzchar(records) & !startsWith(records, "#")]
  mat <- do.call(rbind, strsplit(records, "\t", fixed = TRUE))

  fixed <- data.frame(CHROM=mat[,1], POS=as.numeric(mat[,2]), ID=mat[,3],
                      REF=mat[,4], ALT=mat[,5], stringsAsFactors=FALSE)
  gt <- matrix(NA_character_, nrow(mat), length(samples), dimnames=list(NULL,samples))
  vaf <- matrix(NA_real_, nrow(mat), length(samples), dimnames=list(NULL,samples))
  for (i in seq_len(nrow(mat))) {
    keys <- strsplit(mat[i,9], ":", fixed=TRUE)[[1]]
    gi <- match("GT", keys); ai <- match("AD", keys)
    for (j in seq_along(samples)) {
      values <- strsplit(mat[i,j+9], ":", fixed=TRUE)[[1]]
      if (!is.na(gi) && gi <= length(values)) gt[i,j] <- values[gi]
      if (!is.na(ai) && ai <= length(values)) {
        ad <- suppressWarnings(as.numeric(strsplit(values[ai], ",", fixed=TRUE)[[1]]))
        if (length(ad)>=2 && !anyNA(ad) && sum(ad)>0) vaf[i,j] <- sum(ad[-1])/sum(ad)
      }
    }
  }

  contig_lines <- lines[grep("^##contig=<", lines)]
  contigs <- do.call(rbind, lapply(contig_lines, function(x) data.frame(
    CHROM=sub(".*ID=([^,>]+).*", "\\1", x),
    Length=as.numeric(sub(".*length=([0-9]+).*", "\\1", x)),
    stringsAsFactors=FALSE)))
  if (is.null(contigs) || anyNA(contigs$Length)) {
    contigs <- aggregate(POS~CHROM, fixed, max); names(contigs)[2] <- "Length"
  }
  fixed$CHROM <- factor(fixed$CHROM, levels=contigs$CHROM)
  list(fixed=fixed, gt=gt, vaf=vaf, samples=samples, contigs=contigs)
}

haploid_gt <- function(x) {
  z <- rep(NA_integer_, length(x))
  z[x %in% c("0","0/0","0|0")] <- 0L
  z[x %in% c("1","1/1","1|1")] <- 1L
  z
}

informative_mask <- function(vcf, ref_parent="11730", alt_parent="22971") {
  if (!all(c(ref_parent,alt_parent) %in% vcf$samples)) stop("Parent sample missing")
  haploid_gt(vcf$gt[,ref_parent])==0L & haploid_gt(vcf$gt[,alt_parent])==1L
}

make_frequency_table <- function(vcf, mask, prefix, out_dir=".") {
  out <- cbind(vcf$fixed[mask,c("CHROM","POS","ID","REF","ALT")],
               as.data.frame(vcf$vaf[mask,,drop=FALSE], check.names=FALSE))
  out_file <- file.path(out_dir, paste0(prefix,"_alt_allele_frequency_1391.tsv"))
  write.table(out, out_file, sep="\t", quote=FALSE, row.names=FALSE, na="NA")
  out
}

prepare_horizontal_map <- function(freq, contigs, samples, labels) {
  missing <- setdiff(samples, names(labels))
  labels[missing] <- missing
  long <- do.call(rbind, lapply(samples, function(s) data.frame(
    CHROM=as.character(freq$CHROM), POS=freq$POS, Sample=s,
    Label=unname(labels[s]), AlleleFreq=freq[[s]], stringsAsFactors=FALSE)))
  long$ChrY <- length(contigs$CHROM)-match(long$CHROM,contigs$CHROM)+1
  bg <- do.call(rbind, lapply(samples, function(s) data.frame(
    CHROM=contigs$CHROM, Length=contigs$Length, Sample=s,
    Label=unname(labels[s]),
    ChrY=length(contigs$CHROM)-seq_along(contigs$CHROM)+1,
    stringsAsFactors=FALSE)))
  list(long=long,bg=bg)
}

horizontal_frequency_map <- function(long,bg,contigs,title_text,small=FALSE) {
  fs <- if(small) 11 else 12
  ggplot2::ggplot() +
    ggplot2::geom_rect(data=bg,
      ggplot2::aes(xmin=0,xmax=Length/1e6, ymin=ChrY-0.32,ymax=ChrY+0.32),
      fill="white",color="black",linewidth=0.4) +
    ggplot2::geom_point(data=long[is.na(long$AlleleFreq),],
      ggplot2::aes(POS/1e6,ChrY), shape=124,size=if(small)4.2 else 5.2,color="grey70") +
    ggplot2::geom_point(data=long[!is.na(long$AlleleFreq),],
      ggplot2::aes(POS/1e6,ChrY,color=AlleleFreq), shape=124,size=if(small)4.2 else 5.2) +
    ggplot2::scale_color_gradient2(low="#8B2671",mid="#FFCA28",high="#069F72",
      midpoint=0.5,limits=c(0,1),na.value="grey70",name="ALT frequency") +
    ggplot2::scale_y_continuous(
      breaks=rev(seq_along(contigs$CHROM)),labels=paste0("Chr",contigs$CHROM),
      limits=c(0.5,length(contigs$CHROM)+0.5)) +
    ggplot2::scale_x_continuous(
      limits=c(0,max(contigs$Length)/1e6),expand=c(0.01,0)) +
    ggplot2::labs(title=title_text,x="Chromosomal position (Mb)",y=NULL) +
    ggplot2::theme_bw(base_size=fs) + ggplot2::theme(
      panel.grid=ggplot2::element_blank(),
      axis.text=ggplot2::element_text(color="black"),
      axis.text.y=ggplot2::element_text(face="bold",size=fs),
      plot.title=ggplot2::element_text(hjust=0.5,face="bold",size=fs+2),
      legend.position="right",plot.margin=ggplot2::margin(8,10,8,10))
}

plot_frequency_maps <- function(freq,contigs,dataset,prefix,combined_groups,
                                sample_labels=NULL, out_dir=".") {
  all_samples <- setdiff(names(freq),c("CHROM","POS","ID","REF","ALT"))
  if(is.null(sample_labels)) sample_labels <- setNames(all_samples,all_samples)
  
  outdir <- file.path(out_dir, paste0(prefix,"_individual_allele_frequency"))
  dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
  
  all_data <- prepare_horizontal_map(freq,contigs,all_samples,sample_labels)

  pdf(file.path(out_dir, paste0(prefix,"_allele_frequency_all_samples.pdf")),
      width=8,height=5.2,onefile=TRUE)
  for(s in all_samples) {
    p <- horizontal_frequency_map(
      all_data$long[all_data$long$Sample==s,],
      all_data$bg[all_data$bg$Sample==s,],contigs,unname(sample_labels[s]),FALSE)
    print(p)
    ggplot2::ggsave(file.path(outdir,paste0(s,"_allele_frequency.png")),p,width=8,height=5.2,dpi=300)
    ggplot2::ggsave(file.path(outdir,paste0(s,"_allele_frequency.pdf")),p,width=8,height=5.2)
  }
  dev.off()

  for(group_name in names(combined_groups)) {
    requested_slots <- combined_groups[[group_name]]
    blank_slots <- which(is.na(requested_slots))
    ss <- requested_slots[!is.na(requested_slots)]
    missing <- setdiff(ss,all_samples)
    if(length(missing)) stop(group_name," samples missing: ",paste(missing,collapse=", "))
    d <- prepare_horizontal_map(freq,contigs,ss,sample_labels)

    lev <- character(length(requested_slots))
    lev[!is.na(requested_slots)] <- unname(sample_labels[ss])
    if(length(blank_slots)) {
      lev[blank_slots] <- paste0(".__BLANK_FACET_",seq_along(blank_slots),"__")
    }
    d$long$Label <- factor(d$long$Label,levels=lev)
    d$bg$Label <- factor(d$bg$Label,levels=lev)
    facet_labels <- setNames(lev,lev)
    if(length(blank_slots)) facet_labels[lev[blank_slots]] <- ""
    
    p <- horizontal_frequency_map(d$long,d$bg,contigs,paste(dataset,group_name),TRUE) +
      ggplot2::facet_wrap(~Label,ncol=3,drop=FALSE,
        labeller=ggplot2::as_labeller(facet_labels)) + 
      ggplot2::theme(
        strip.text=ggplot2::element_text(size=12,face="bold"),
        axis.title.x=ggplot2::element_text(size=12),
        legend.title=ggplot2::element_text(size=11),
        legend.text=ggplot2::element_text(size=10))
        
    nr <- ceiling(length(requested_slots)/3); h <- nr*4.0
    stem <- file.path(out_dir, paste0(prefix,"_",group_name,"_allele_frequency_combined"))
    ggplot2::ggsave(paste0(stem,".png"),p,width=16,height=h,dpi=300,limitsize=FALSE)
    ggplot2::ggsave(paste0(stem,".pdf"),p,width=16,height=h,limitsize=FALSE)
  }
}

calculate_crossovers <- function(vcf,mask,sample_meta,prefix,threshold=0.8, out_dir=".") {
  required <- c("WGA","Line","DPI","Sample_info")
  if(!all(required%in%names(sample_meta))) stop("Incomplete sample metadata")
  missing <- setdiff(sample_meta$WGA,vcf$samples)
  if(length(missing)) stop("Crossover samples missing: ",paste(missing,collapse=", "))
  fixed <- vcf$fixed[mask,c("CHROM","POS")]
  vf <- vcf$vaf[mask,sample_meta$WGA,drop=FALSE]
  events <- list(); summaries <- list()

  for(s in sample_meta$WGA) {
    ee <- list(); callable <- 0L
    for(chr in levels(fixed$CHROM)) {
      idx <- which(fixed$CHROM==chr); idx <- idx[order(fixed$POS[idx])]
      if(length(idx)<2) next
      a <- vf[idx[-length(idx)],s]
      b <- vf[idx[-1],s]
      valid <- !is.na(a)&!is.na(b)
      delta <- b-a
      hit <- valid & abs(delta)>threshold
      callable <- callable+sum(valid)
      k <- which(hit)
      if(length(k)) ee[[length(ee)+1]] <- data.frame(
        WGA=s,CHROM=chr,
        SNP1_POS=fixed$POS[idx[k]],SNP2_POS=fixed$POS[idx[k+1]],
        Event_POS=fixed$POS[idx[k]],
        Interval_bp=fixed$POS[idx[k+1]]-fixed$POS[idx[k]],
        SNP1_VAF=a[k],SNP2_VAF=b[k],Delta_VAF=delta[k],
        Abs_Delta_VAF=abs(delta[k]),stringsAsFactors=FALSE)
    }
    e <- if(length(ee))do.call(rbind,ee)else NULL
    if(!is.null(e))events[[length(events)+1]] <- e
    summaries[[length(summaries)+1]] <- data.frame(
      WGA=s,Callable_adjacent_pairs=callable,
      Crossover_count=if(is.null(e))0L else nrow(e))
  }
  
  event_df <- if(length(events)) do.call(rbind,events) else data.frame()
  summary_df <- do.call(rbind,summaries)
  
  if(nrow(event_df) > 0) {
      event_df <- merge(sample_meta,event_df,by="WGA",all.y=TRUE,sort=FALSE)
  }
  summary_df <- merge(sample_meta,summary_df,by="WGA",all.y=TRUE,sort=FALSE)
  summary_df <- summary_df[match(sample_meta$WGA,summary_df$WGA),]
  
  write.table(event_df,file.path(out_dir,paste0(prefix,"_crossover_events_VAFdiff_gt0.8_1391.tsv")),
    sep="\t",quote=FALSE,row.names=FALSE)
  write.table(summary_df,file.path(out_dir,paste0(prefix,"_crossover_summary_VAFdiff_gt0.8_1391.tsv")),
    sep="\t",quote=FALSE,row.names=FALSE)

  if(nrow(event_df) > 0) {
      ct <- vcf$contigs
      ct$Offset <- c(0,head(cumsum(ct$Length),-1))
      event_df$GenomePos <- event_df$Event_POS+ct$Offset[match(event_df$CHROM,ct$CHROM)]
      mids <- ct$Offset+ct$Length/2
      ends <- head(cumsum(ct$Length),-1)
      event_df$Sample_info <- factor(event_df$Sample_info,levels=rev(sample_meta$Sample_info))
      p <- ggplot2::ggplot(event_df,ggplot2::aes(GenomePos/1e6,Sample_info)) +
        ggplot2::geom_point(shape=124,size=4,color="#8B2671") +
        ggplot2::geom_vline(xintercept=ends/1e6,linetype="dotted",color="grey65") +
        ggplot2::scale_x_continuous(breaks=mids/1e6,labels=paste0("Chr",ct$CHROM)) +
        ggplot2::labs(x="Chromosome",y=NULL,
          title="Crossover positions: |VAF2 - VAF1| > 0.8") +
        ggplot2::theme_classic(base_size=12) +
        ggplot2::theme(axis.text.y=ggplot2::element_text(size=9))
      ggplot2::ggsave(file.path(out_dir,paste0(prefix,"_crossover_positions_VAFdiff_gt0.8_1391.png")),p,
        width=12,height=6.5,dpi=300)
      ggplot2::ggsave(file.path(out_dir,paste0(prefix,"_crossover_positions_VAFdiff_gt0.8_1391.pdf")),p,
        width=12,height=6.5)
  }
  invisible(list(events=event_df,summary=summary_df))
}