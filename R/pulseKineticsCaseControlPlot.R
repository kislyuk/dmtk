library(plyr)
library(ggplot2)

data_fn = "@inputData"
pdf("@outputFile", width=24, height=6)
xmin=@xmin; xmax=@xmax; ymin=@ymin; ymax=@ymax
ylabel = "@ylabel"
title = "@title"
ref_name = "@ref_name"
xlabel = paste("Reference position:", ref_name);
baseline = @baseline
min_target_cov=30
max_plot_width=500; max_plot_pages=99
xstep=20; ystep=1
print_sequence = TRUE; line_offset = 0.025
top_strand_label="@top_strand_label"; bottom_strand_label="@bottom_strand_label";
data <- read.csv(data_fn)
#if (is.null(data$spreadLo)) data$spreadLo = data$myVar
#if (is.null(data$spreadHi)) data$spreadHi = data$myVar
if (is.null(xmin)) xmin = floor(min(data$RefPosn))
if (is.null(xmax) || xmin >= xmax) xmax = ceiling(max(data$RefPosn))
xmax = min(xmax, xmin + max_plot_width*max_plot_pages)
# if (is.null(ymin)) ymin = floor(min(data$spreadLo, 0, na.rm=TRUE)); if (is.null(ymax)) ymax = ceiling(max(data$spreadHi, na.rm=TRUE))
# if (is.null(ymin)) ymin = floor(min(data$myVar, 0, na.rm=TRUE)); if (is.null(ymax)) ymax = ceiling(max(data$myVar, na.rm=TRUE))
if (is.null(ymin)) {
    hi_cov_data = subset(data, effectiveCoverageFwdStrand >= min_target_cov)
    if (is.null(hi_cov_data$myVarFwdStrand)) hi_cov_data = data
    ymin = floor(min(c(hi_cov_data$myVarFwdStrand, hi_cov_data$myVarRevStrand), na.rm=TRUE))
}
if (is.null(ymax) || ymin >= ymax) ymax = ceiling(max(c(data$myVarFwdStrand, data$myVarRevStrand), na.rm=TRUE))
if (ymax-ymin < 20) ystep = 1 else if (ymax-ymin < 200) ystep = 10 else if (ymax-ymin < 2000) ystep = 100 else if (ymax-ymin < 20000) ystep = 1000 else if (ymax-ymin < 200000) ystep = 10000 else ystep = 100000
ymin = floor(ymin/ystep)*ystep
ymax = ceiling(ymax/ystep)*ystep
xmin = floor(xmin/xstep)*xstep
xmax = ceiling(xmax/xstep)*xstep
mod_hilite_fwd = data.frame(cbind(x=data$RefPosn[data$ModifiedFwdStrand==1], y=baseline))
mod_hilite_rev = data.frame(cbind(x=data$RefPosn[data$ModifiedRevStrand==1], y=baseline))
if (is.null(mod_hilite_fwd$x)) mod_hilite_fwd$x = -1
if (is.null(mod_hilite_rev$x)) mod_hilite_rev$x = -1

data$myVarFwdStrand[sapply(data$myVarFwdStrand, is.na)] = baseline
data$myVarRevStrand[sapply(data$myVarRevStrand, is.na)] = baseline

for (x_plot_lo in seq(xmin, xmax, max_plot_width)) {
	x_plot_hi = x_plot_lo + max_plot_width + xstep
	if (xmax < max_plot_width) x_plot_hi = xmax
	data_this_pane = subset(data, RefPosn >= x_plot_lo & RefPosn <= x_plot_hi)
	data_low_cov_fwd = subset(data_this_pane, effectiveCoverageFwdStrand < min_target_cov)
	data_low_cov_rev = subset(data_this_pane, effectiveCoverageRevStrand < min_target_cov)
	p1 <- ggplot(data_this_pane) +
		geom_rect(data=mod_hilite_rev, aes(xmin=x-0.3, xmax=x+0.3, ymin=ymin, ymax=ymax), colour="grey80", fill="grey80") +
		geom_rect(aes(xmin=RefPosn-0.4, xmax=RefPosn+0.4, ymin=baseline, ymax=myVarRevStrand, fill=factor(ModifiedRevStrand))) +
#                        geom_errorbar(aes(x=RefPosn, ymin=spreadLoRevStrand, ymax=spreadHiRevStrand), width=0.4) +
		scale_x_continuous(breaks = seq(x_plot_lo, x_plot_hi, xstep)) +
		scale_y_continuous(breaks = seq(ymin, ymax, ystep)) +
		coord_cartesian(ylim = c(ymin, ymax), xlim=c(x_plot_lo, x_plot_hi)) +
		ylab(ylabel) + xlab("") +
		geom_text(aes(x=x_plot_lo, y=ymax, hjust=0, vjust=1, label=top_strand_label), colour="grey30") +
		theme_gray(base_size = 18) +
#		scale_fill_gradient2(low="grey", mid="red", high="blue", midpoint=500) +
#		scale_fill_hue(c=mapply(min, 100, data_this_pane$effectiveCoverageRevStrand*100/1000)) +
		opts(legend.position = "none", title = title, plot.margin = unit(c(0, 2, -.5, 0), "lines"), axis.text.x = theme_text(size=0))
	if (!empty(data_low_cov_rev)) {
	    p1 <- p1 + geom_rect(data=data_low_cov_rev, aes(xmin=RefPosn-0.4, xmax=RefPosn+0.4, ymin=baseline, ymax=myVarRevStrand), colour="grey", fill="grey", alpha=0.8)
	}
	if (print_sequence) {
        p1 <- p1 + geom_text(aes(x=RefPosn, y=ymin+(ymax-ymin)*line_offset, vjust=0, label=nucleotide), size=2, color="grey")
        p1 <- p1 + geom_text(aes(x=RefPosn, y=ymin, vjust=0, label=refNucleotide), size=2, color="grey")
        seq_labels <- data.frame(x=x_plot_lo, y=c(ymin, ymin+(ymax-ymin)*line_offset), text=c("Ref:", "Cons:"))
        p1 <- p1 + geom_text(aes(x=x, y=y, label=text, vjust=0, hjust=0), data=seq_labels, size=2)
    }
#    p1 <- p1 + geom_line(aes(x=RefPosn, y=effectiveCoverageFwdStrand))
	p2 <- ggplot(data_this_pane) +
		geom_rect(data=mod_hilite_fwd, aes(xmin=x-0.3, xmax=x+0.3, ymin=ymin, ymax=ymax), colour="grey80", fill="grey80") +
		geom_rect(aes(xmin=RefPosn-0.4, xmax=RefPosn+0.4, ymin=baseline, ymax=myVarFwdStrand, fill=factor(ModifiedFwdStrand))) +
#                        geom_errorbar(aes(x=RefPosn, ymin=spreadLoFwdStrand, ymax=spreadHiFwdStrand), width=0.4) +
		scale_x_continuous(breaks = seq(x_plot_lo, x_plot_hi, xstep)) +
		scale_y_reverse(breaks = seq(ymin, ymax, ystep)) +
		coord_cartesian(ylim = c(ymin, ymax), xlim=c(x_plot_lo, x_plot_hi)) +
		ylab(ylabel) + xlab(paste(xlabel, " (", x_plot_lo, "..", x_plot_hi, ")", sep="")) +
		geom_text(aes(x=x_plot_lo, y=ymax, hjust=0, vjust=0, label=bottom_strand_label), colour="grey30") +
		theme_gray(base_size = 18) +
		opts(legend.position = "none", plot.margin = unit(c(-.5, 2, 0, 0), "lines"))
    if (!empty(data_low_cov_fwd)) {
        p2 <- p2 + geom_rect(data=data_low_cov_fwd, aes(xmin=RefPosn-0.4, xmax=RefPosn+0.4, ymin=baseline, ymax=myVarFwdStrand), colour="grey", fill="grey", alpha=0.8)
    }
    if (print_sequence) {
        p2 <- p2 + geom_text(aes(x=RefPosn, y=ymin+(ymax-ymin)*line_offset, vjust=1, label=nucleotideRevStrand), size=2, color="grey")
        p2 <- p2 + geom_text(aes(x=RefPosn, y=ymin, vjust=1, label=refNucleotideRevStrand), size=2, color="grey")
        seq_labels <- data.frame(x=x_plot_lo, y=c(ymin, ymin+(ymax-ymin)*line_offset), text=c("Ref:", "Cons:"))
        p2 <- p2 + geom_text(aes(x=x, y=y, label=text, vjust=1, hjust=0), data=seq_labels, size=2)
    }

	grid.newpage()
	pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1)))
	print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
	print(p2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
}
dev.off()
unlink(code_fn); unlink(data_fn)
