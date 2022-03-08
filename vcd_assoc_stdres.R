library(vcd)

legend_res <- structure(function (fontsize = 10.5, fontfamily = "", x = unit(1, 
    "lines"), y = unit(0.1, "npc"), height = unit(0.8, "npc"), 
    width = unit(0.7, "lines"), digits = 2, check_overlap = TRUE, 
    text = NULL, steps = 200, ticks = 10, pvalue = TRUE, range = NULL) 
{
    if (!is.unit(x)) 
        x <- unit(x, "native")
    if (!is.unit(y)) 
        y <- unit(y, "npc")
    if (!is.unit(width)) 
        width <- unit(width, "lines")
    if (!is.unit(height)) 
        height <- unit(height, "npc")
    function(residuals, shading, autotext) {
        res <- as.vector(residuals)
        if (is.null(text)) 
            text <- autotext
        p.value <- attr(shading, "p.value")
        legend <- attr(shading, "legend")
        if (all(residuals == 0)) {
            pushViewport(viewport(x = x, y = y, just = c("left", 
                "bottom"), default.units = "native", height = height, 
                width = width))
            grid.lines(y = 0.5)
            grid.text(0, x = unit(1, "npc") + unit(0.8, "lines"), 
                y = 0.5, gp = gpar(fontsize = fontsize, fontfamily = fontfamily))
            warning("All residuals are zero.")
        }
        else {
            if (is.null(range)) 
                range <- range(res)
            if (length(range) != 2) 
                stop("Range must have length two!")
            if (is.na(range[1])) 
                range[1] <- min(res)
            if (is.na(range[2])) 
                range[2] <- max(res)
            pushViewport(viewport(x = x, y = y, just = c("left", 
                "bottom"), yscale = range, default.units = "native", 
                height = height, width = width))
            if (is.null(legend$col.bins)) {
                col.bins <- seq(range[1], range[2], length = steps)
                at <- NULL
            }
            else {
                col.bins <- sort(unique(c(legend$col.bins, range)))
                col.bins <- col.bins[col.bins <= range[2] & col.bins >= 
                  range[1]]
                at <- col.bins
            }
            y.pos <- col.bins[-length(col.bins)]
            y.height <- diff(col.bins)
            grid.rect(x = unit(rep.int(0, length(y.pos)), "npc"), 
                y = y.pos, height = y.height, default.units = "native", 
                gp = gpar(fill = shading(y.pos + 0.5 * y.height)$fill, 
                  col = 0), just = c("left", "bottom"))
            grid.rect(gp = gpar(fill = "transparent"))
            if (is.null(at)) 
                at <- seq(from = head(col.bins, 1), to = tail(col.bins, 
                  1), length = ticks)
            lab <- format(round(at, digits = digits), nsmall = digits)
            tw <- lab[which.max(nchar(lab))]
            grid.text(format(signif(at, digits = digits)), x = unit(1, 
                "npc") + unit(0.8, "lines") + unit(1, "strwidth", 
                tw), y = at, default.units = "native", just = c("right", 
                "center"), gp = gpar(fontsize = fontsize, fontfamily = fontfamily), 
                check.overlap = check_overlap)
            grid.segments(x0 = unit(1, "npc"), x1 = unit(1, "npc") + 
                unit(0.5, "lines"), y0 = at, y1 = at, default.units = "native")
        }
        popViewport(1)
        grid.text(text, x = x, y = unit(1, "npc") - y + unit(1, 
            "lines"), gp = gpar(fontsize = fontsize, fontfamily = fontfamily, 
            lineheight = 0.8), just = c("left", "bottom"))
        if (!is.null(p.value) && pvalue) {
            grid.text(p.value, x = x, y = y - unit(1, "lines"), gp = gpar(fontsize = fontsize, 
                fontfamily = fontfamily, lineheight = 0.8), just = c("left", 
                "top"))
        }
    }
}, class = "grapcon_generator")


shading_hcl2 <- structure(function (observed, residuals = NULL, expected = NULL, 
    df = NULL, h = NULL, c = NULL, l = NULL, interpolate = c(2, 
        4), lty = 1, eps = NULL, line_col = "black", p.value = NULL, 
    level = 0.95, p = NULL, ...) 
{
    if (is.null(h)) 
        h <- c(260, 0)
    if (is.null(c)) 
        c <- c(100, 20)
    if (is.null(l)) 
        l <- c(90, 50)
    my.h <- rep(h, length.out = 2)
    my.c <- rep(c, length.out = 2)
    my.l <- rep(l, length.out = 2)
    lty <- rep(lty, length.out = 2)
    if (is.null(expected) && !is.null(residuals)) 
        stop("residuals without expected values specified")
    if (!is.null(expected) && is.null(df) && is.null(p.value)) {
        warning("no default inference available without degrees of freedom")
        p.value <- NA
    }
    if (is.null(expected) && !is.null(observed)) {
        expected <- loglin(observed, 1:length(dim(observed)), 
            fit = TRUE, print = FALSE)
        df <- expected$df
        expected <- expected$fit
    }
    if (is.null(residuals) && !is.null(observed)) 
        residuals <- (observed - expected)/sqrt(expected)
    if (is.null(p.value)) 
        p.value <- function(observed, residuals, expected, df) pchisq(sum(as.vector(residuals)^2), 
            df, lower.tail = FALSE)
    if (!is.function(p.value) && is.na(p.value)) {
        max.c <- my.c[1]
        p.value <- NULL
    }
    else {
        if (is.function(p.value)) 
            p.value <- p.value(observed, residuals, expected, 
                df)
        max.c <- ifelse(p.value < (1 - level), my.c[1], my.c[2])
    }
    if (!is.function(interpolate)) {
        col.bins <- sort(interpolate)
        interpolate <- stepfun(col.bins, seq(0, 1, length = length(col.bins) + 
            1))
        col.bins <- sort(unique(c(col.bins, 0, -col.bins)))
    }
    else {
        col.bins <- NULL
    }
    legend <- NULL
    if (!is.null(col.bins)) {
        res2 <- col.bins
        res2 <- c(head(res2, 1) - 1, res2[-1] - diff(res2)/2, 
            tail(res2, 1) + 1)
        legend.col <- hcl2hex(ifelse(res2 > 0, my.h[1], my.h[2]), 
            max.c * pmax(pmin(interpolate(abs(res2)), 1), 0), 
            my.l[1] + diff(my.l) * pmax(pmin(interpolate(abs(res2)), 
                1), 0), ...)
        lty.bins <- 0
        legend.lty <- lty[2:1]
        legend <- list(col = legend.col, col.bins = col.bins, 
            lty = legend.lty, lty.bins = lty.bins)
    }
    rval <- function(x) {
        res <- as.vector(x)
        fill <- hcl2hex(ifelse(res > 0, my.h[1], my.h[2]), max.c * 
            pmax(pmin(interpolate(abs(res)), 1), 0), my.l[1] + 
            diff(my.l) * pmax(pmin(interpolate(abs(res)), 1), 
                0), ...)
        dim(fill) <- dim(x)
        col <- rep(line_col, length.out = length(res))
        if (!is.null(eps)) {
            eps <- abs(eps)
            col[res > eps] <- hcl2hex(my.h[1], max.c, my.l[2], 
                ...)
            col[res < -eps] <- hcl2hex(my.h[2], max.c, my.l[2], 
                ...)
        }
        dim(col) <- dim(x)
        ltytmp <- ifelse(x > 0, lty[1], lty[2])
        if (!is.null(eps)) 
            ltytmp[abs(x) < abs(eps)] <- lty[1]
        dim(ltytmp) <- dim(x)
        return(structure(list(col = col, fill = fill, lty = ltytmp), 
            class = "gpar"))
    }
    attr(rval, "legend") <- legend
    attr(rval, "p.value") <- paste("p-value =\n" ,format.pval(p))
    return(rval)
}, class = "grapcon_generator")

#function structplot modified
struct <- function (x, residuals = NULL, expected = NULL, condvars = NULL, 
    shade = NULL, type = c("observed", "expected"), 
    residuals_type = NULL, df = NULL, split_vertical = NULL, 
    spacing = spacing_equal, spacing_args = list(), gp = NULL, 
    gp_args = list(), labeling = labeling_border, labeling_args = list(), 
    core = struc_mosaic, core_args = list(), legend = NULL, legend_args = list(), 
    main = NULL, sub = NULL, margins = unit(3, "lines"), 
    title_margins = NULL, legend_width = NULL, main_gp = gpar(fontsize = 20), 
    sub_gp = gpar(fontsize = 15), newpage = TRUE, pop = TRUE, 
    return_grob = FALSE, keep_aspect_ratio = NULL, prefix = "", p = p,
    ...) 
{
    if (is.null(shade)) 
        shade <- !is.null(gp) || !is.null(expected)
    type <- match.arg(type)
    if (is.null(residuals)) {
        residuals_type <- if (is.null(residuals_type)) 
            "pearson"
        else match.arg(tolower(residuals_type), c("pearson", 
            "deviance", "ft"))
    }
    else {
        if (is.null(residuals_type)) 
            residuals_type <- ""
    }
    if (is.structable(x)) {
        if (is.null(split_vertical)) 
            split_vertical <- attr(x, "split_vertical")
        x <- as.table(x)
    }
    if (is.null(split_vertical)) 
        split_vertical <- FALSE
    d <- dim(x)
    dl <- length(d)
    dn <- dimnames(x)
    if (is.null(dn)) 
        dn <- dimnames(x) <- lapply(d, seq)
    dnn <- names(dimnames(x))
    if (is.null(dnn)) 
        dnn <- names(dn) <- names(dimnames(x)) <- LETTERS[1:dl]
    if (any(nas <- is.na(x))) 
        x[nas] <- 0
    if ((is.null(expected) && is.null(residuals)) || !is.numeric(expected)) {
        if (!is.null(df)) 
            warning("Using calculated degrees of freedom.")
        if (inherits(expected, "formula")) {
            fm <- loglm(expected, x, fitted = TRUE)
            expected <- fitted(fm)
            df <- fm$df
        }
        else {
            if (is.null(expected)) 
                expected <- if (is.null(condvars)) 
                  as.list(1:dl)
                else lapply((condvars + 1):dl, c, seq(condvars))
            fm <- loglin(x, expected, fit = TRUE, print = FALSE)
            expected <- fm$fit
            df <- fm$df
        }
    }
    if (is.null(residuals)) 
        residuals <- switch(residuals_type, pearson = (x - expected)/sqrt(ifelse(expected > 
            0, expected, 1)), deviance = {
            tmp <- 2 * (x * log(ifelse(x == 0, 1, x/ifelse(expected > 
                0, expected, 1))) - (x - expected))
            tmp <- sqrt(pmax(tmp, 0))
            ifelse(x > expected, tmp, -tmp)
        }, ft = sqrt(x) + sqrt(x + 1) - sqrt(4 * expected + 1))
    if (any(nas <- is.na(residuals))) 
        residuals[nas] <- 0
    if (length(split_vertical) == 1) 
        split_vertical <- rep(c(split_vertical, !split_vertical), 
            length.out = dl)
    if (is.null(keep_aspect_ratio)) 
        keep_aspect_ratio <- dl < 3
    if (is.function(spacing)) {
        if (inherits(spacing, "grapcon_generator")) 
            spacing <- do.call("spacing", spacing_args)
        spacing <- spacing(d, condvars)
    }
    if (shade) {
        if (is.null(gp)) 
            gp <- shading_hcl2(observed = x, p = p)
        if (is.function(gp)) {
            if (is.null(legend) || (is.logical(legend) && legend)) 
                legend <- legend_res
            gpfun <- if (inherits(gp, "grapcon_generator")) 
                do.call("gp", c(list(x, residuals, expected, 
                  df), as.list(gp_args)))
            else gp
            gp <- gpfun(residuals)
        }
        else if (!is.null(legend) && !(is.logical(legend) && 
            !legend)) 
            stop("gp argument must be a shading function for drawing a legend")
    }
    else {
        if (!is.null(gp)) {
            warning("gp parameter ignored since shade = FALSE")
            gp <- NULL
        }
    }
    if (is.null(gp)) 
        gp <- gpar(fill = grey(0.8))
    size <- prod(d)
    FUN <- function(par) {
        if (is.structable(par)) 
            par <- as.table(par)
        if (length(par) < size || is.null(dim(par))) 
            array(par, dim = d)
        else par
    }
    gp <- structure(lapply(gp, FUN), class = "gpar")
    if (newpage) 
        grid.newpage()
    if (keep_aspect_ratio) 
        pushViewport(viewport(width = 1, height = 1, default.units = "snpc"))
    pushViewport(vcdViewport(mar = margins, oma = title_margins, 
        legend = shade && !(is.null(legend) || is.logical(legend) && 
            !legend), main = !is.null(main), sub = !is.null(sub), 
        keep_aspect_ratio = keep_aspect_ratio, legend_width = legend_width, 
        prefix = prefix))
    if (inherits(legend, "grapcon_generator")) 
        legend <- do.call("legend", legend_args)
    if (shade && !is.null(legend) && !(is.logical(legend) && 
        !legend)) {
        seekViewport(paste(prefix, "legend", sep = ""))
        residuals_type <- switch(residuals_type, deviance = "deviance\nresiduals:", 
            ft = "Freeman-Tukey\nresiduals:", pearson = "Pearson\nresiduals:", 
            residuals_type)
        legend(residuals, gpfun, "standardized\nresiduals:")
    }
    if (!is.null(main)) {
        seekViewport(paste(prefix, "main", sep = ""))
        if (is.logical(main) && main) 
            main <- deparse(substitute(x))
        grid.text(main, gp = main_gp)
    }
    if (!is.null(sub)) {
        seekViewport(paste(prefix, "sub", sep = ""))
        if (is.logical(sub) && sub && is.null(main)) 
            sub <- deparse(substitute(x))
        grid.text(sub, gp = sub_gp)
    }
    seekViewport(paste(prefix, "plot", sep = ""))
    if (inherits(core, "grapcon_generator")) 
        core <- do.call("core", core_args)
    core(residuals = residuals, observed = if (type == "observed") 
        x
    else expected, expected = if (type == "observed") 
        expected
    else x, spacing = spacing, gp = gp, split_vertical = split_vertical, 
        prefix = prefix)
    upViewport(1)
    if (is.logical(labeling)) 
        labeling <- if (labeling) 
            labeling_border
        else NULL
    if (!is.null(labeling)) {
        if (inherits(labeling, "grapcon_generator")) 
            labeling <- do.call("labeling", c(labeling_args, 
                list(...)))
        labeling(dn, split_vertical, condvars, prefix)
    }
    seekViewport(paste(prefix, "base", sep = ""))
    if (pop) 
        popViewport(1 + keep_aspect_ratio)
    else upViewport(1 + keep_aspect_ratio)
    if (return_grob) 
        invisible(structure(structable(if (type == "observed") x else expected, 
            split_vertical = split_vertical), grob = grid.grab()))
    else invisible(structable(if (type == "observed") x else expected, 
        split_vertical = split_vertical))}
vcdViewport <- function(mar = rep.int(2.5, 4),
                        legend_width = unit(5, "lines"),
                        oma = NULL,
                        legend = FALSE, main = FALSE, sub = FALSE,
                        keep_aspect_ratio = TRUE,
                        prefix = "")


{
  
  if (is.null(legend_width))
    legend_width <- unit(5 * legend, "lines")
  if (!is.unit(legend_width))
    legend_width <- unit(legend_width, "lines")

  if (legend && !main && !sub && keep_aspect_ratio) main <- sub <- TRUE
  mar <- if (!is.unit(mar))
    unit(pexpand(mar, 4, rep.int(2.5, 4), c("top","right","bottom","left")), "lines")
  else
    rep(mar, length.out = 4)
  if (is.null(oma)) {
    space <- if (legend && keep_aspect_ratio)
      legend_width + mar[2] + mar[4] - mar[1] - mar[3]
    else unit(0, "lines")
    oma <- if (main && sub)
      max(unit(2, "lines"), 0.5 * space)
    else if (main)
      unit.c(max(unit(2, "lines"), space), unit(0, "lines"))
    else if (sub)
      unit.c(unit(0, "lines"), max(unit(2, "lines"), space))
    else
      0.5 * space
  }
  oma <- if (!is.unit(oma))
    unit(pexpand(oma, 2, rep.int(2, 2), c("top","bottom")), "lines")
  else
    rep(oma, length.out = 2)

  
  vpPlot <- vpStack(viewport(layout.pos.col = 2, layout.pos.row = 3),
                    viewport(width = 1, height = 1, name = paste(prefix, "plot", sep = ""),
                             default.units = if (keep_aspect_ratio) "snpc" else "npc"))
  vpMarginBottom <- viewport(layout.pos.col = 2, layout.pos.row = 4,
                             name = paste(prefix, "margin_bottom", sep = ""))
  vpMarginLeft <- viewport(layout.pos.col = 1, layout.pos.row = 3,
                           name = paste(prefix, "margin_left", sep = ""))
  vpMarginTop <- viewport(layout.pos.col = 2, layout.pos.row = 2,
                          name = paste(prefix, "margin_top", sep = ""))
  vpMarginRight <- viewport(layout.pos.col = 3, layout.pos.row = 3,
                            name = paste(prefix, "margin_right", sep = ""))
  vpCornerTL <- viewport(layout.pos.col = 1, layout.pos.row = 2,
                         name = paste(prefix, "corner_top_left", sep = ""))
  vpCornerTR <- viewport(layout.pos.col = 3, layout.pos.row = 2,
                         name = paste(prefix, "corner_top_right", sep = ""))
  vpCornerBL <- viewport(layout.pos.col = 1, layout.pos.row = 4,
                         name = paste(prefix, "corner_bottom_left", sep = ""))
  vpCornerBR <- viewport(layout.pos.col = 3, layout.pos.row = 4,
                         name = paste(prefix, "corner_bottom_right", sep = ""))

  vpLegend <- viewport(layout.pos.col = 4, layout.pos.row = 3,
                       name = paste(prefix, "legend", sep = ""))
  vpLegendTop <- viewport(layout.pos.col = 4, layout.pos.row = 2,
                          name = paste(prefix, "legend_top", sep = ""))
  vpLegendSub <- viewport(layout.pos.col = 4, layout.pos.row = 4,
                          name = paste(prefix, "legend_sub", sep = ""))
  vpBase <- viewport(layout = grid.layout(5, 4,
                       widths = unit.c(mar[4], unit(1, "null"), mar[2], legend_width),
                       heights = unit.c(oma[1], mar[1], unit(1, "null"), mar[3], oma[2])),
                     name = paste(prefix, "base", sep = ""))
  vpMain <- viewport(layout.pos.col = 1:4, layout.pos.row = 1,
                     name = paste(prefix, "main", sep = ""))
  vpSub <- viewport(layout.pos.col = 1:4, layout.pos.row = 5,
                    name = paste(prefix, "sub", sep = ""))

  vpTree(vpBase, vpList(vpMain, vpMarginBottom, vpMarginLeft, vpMarginTop,
                        vpMarginRight, vpLegendTop, vpLegend,
                        vpLegendSub, vpCornerTL, vpCornerTR,
                        vpCornerBL, vpCornerBR, vpPlot, vpSub))

}

#function assocplot modified
aspl <- function (x, row_vars = NULL, col_vars = NULL, compress = TRUE, 
                                  xlim = NULL, ylim = NULL, spacing = spacing_conditional(sp = 0), 
                                  spacing_args = list(), split_vertical = NULL, keep_aspect_ratio = FALSE, 
                                  xscale = 0.9, yspace = unit(0.5, "lines"), main = NULL, 
                                  sub = NULL, ..., gp_axis = gpar(lty = 3), simulate = FALSE, B = 2000) 
      {
              if (is.logical(main) && main) 
                      main <- deparse(substitute(x))
                  else if (is.logical(sub) && sub) 
                          sub <- deparse(substitute(x))
                      if (!inherits(x, "ftable")) 
                              x <- structable(x)
                           tab <- as.table(x)
                      if (simulate == TRUE){
                          t2 <- chisq.test(tab, simulate.p.value = T, B = B)
                          df <- NULL
                          }
                      else{
                          t2 <- chisq.test(tab)
                          df <- t2$parameter
                          }
                          p = t2$p.value
                          tab <- t(tab)
                         t3 <- t(t2$stdres)
                          t4 <- t(t2$expected)
                      dl <- length(dim(tab))
                      cond <- rep(TRUE, dl)
                      cond[length(attr(x, "row.vars")) + c(0, length(attr(x, 
                                                                                                                                "col.vars")))] <- FALSE
                      if (inherits(spacing, "grapcon_generator")) 
                              spacing <- do.call("spacing", spacing_args)
                          spacing <- spacing(dim(tab), condvars = which(cond))
                          if (is.null(split_vertical)) 
                                 split_vertical <- attr(x, "split_vertical")
                              
                                  struct(tab, t3, t4, df = df, spacing = spacing, split_vertical = split_vertical, 
                                                        core = struc_assoc(compress = compress, xlim = xlim, 
                                                                           ylim = ylim, yspace = yspace, xscale = xscale, gp_axis = gp_axis), 
                                                        keep_aspect_ratio = keep_aspect_ratio, 
                                                        main = main, sub = sub, p = p, ...)
                              }